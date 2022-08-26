#ifndef BASICMESSAGE_H
#define BASICMESSAGE_H

#ifdef HAVE_SCHED_YIELD
    #include "sched.h"
#endif
#include <atomic>
#include <new>

// The underlying queue below is taken from OpenPA which was originally a part of MPICH.
// The project has not seen updates in years, and the new C/C++ atomics have been used
// to duplicate the functionality in a (hopefully) future-proof way. The following
// explanation is taken directly from the OpenPA source code:

/* The shadow head pointer makes this queue implementation perform much better
   than a standard queue.  Unfortunately, it also makes it a bit non-intuitive
   when read the code.  The following is an excerpt from "Design and Evaluation
   of Nemesis,  a Scalable, Low-Latency, Message-Passing Communication
   Subsystem" by D. Buntinas, G.  Mercier, and W. Gropp that gives an
   explanation:

      A process must access both the head and tail of the queue when it is
      enqueuing an element on an empty queue or when it is dequeuing an element
      that is the last element in the queue. In these cases, if the head and
      tail were in the same cache line, only one L2 cache miss would be
      encountered. If the queue has more elements in it, however, then the
      enqueuer only needs to access the tail, and the dequeuer only needs to
      access the head. If the head and tail were in the same cache line, then
      there would be L2 misses encountered as a result of false sharing each
      time a process enqueues an element after another has been dequeued from
      the same queue, and vice versa. In this case it would be better if the
      head and tail were in separate cache lines.

      Our solution is to put the head and tail in the same cache line and have a
      shadow head pointer in a separate cache line. The shadow head is
      initialized to NULL. The dequeuer uses the shadow head in place of the
      real head except when the shadow head is NULL, meaning that the queue has
      become empty. If the shadow head is NULL when the dequeuer tries to
      dequeue, it checks the value of the real head. If the real head is not
      NULL, meaning that an element has been enqueued on the queue since the
      last time the queue became empty, the dequeuer initializes its shadow head
      to the value of the real head and sets the real head to NULL. In this way,
      only one L2 cache miss is encountered when enqueuing onto an empty queue
      or dequeuing from a queue with one element. And because the tail and
      shadow head are in separate cache lines, there are no L2 cache misses from
      false sharing.

      We found that using a shadow head pointer reduced one-way latency by about
      200 ns on a dual 2 GHz Xeon node.
*/

namespace mdl {

#ifdef __cpp_lib_hardware_interference_size
    using std::hardware_constructive_interference_size;
    using std::hardware_destructive_interference_size;
#else
    // 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned │ ...
    constexpr std::size_t hardware_constructive_interference_size = 128;
    constexpr std::size_t hardware_destructive_interference_size = 128;
#endif

class CXX_Queue_element_hdr_t {
protected:
    friend class CXX_Queue_info_t;
    std::atomic<CXX_Queue_element_hdr_t *> next = nullptr;
public:
    CXX_Queue_element_hdr_t() = default;
    virtual ~CXX_Queue_element_hdr_t() = default;
};
static_assert(std::atomic<CXX_Queue_element_hdr_t *>::is_always_lock_free,"Atomics must be lock free or performance will be unacceptable");

class CXX_Queue_info_t {
    // This is aligned this way so that head and tail are in the same cache line
    struct alignas(2 * sizeof(CXX_Queue_element_hdr_t *)) normal_ {
        std::atomic<CXX_Queue_element_hdr_t *> head = nullptr;
        std::atomic<CXX_Queue_element_hdr_t *> tail = nullptr;
    } normal;
    // The shadow head should be in a different cache line to prevent false sharing
    char pad[hardware_destructive_interference_size - sizeof(normal)];
    struct shadow_ {
        std::atomic<CXX_Queue_element_hdr_t *> head = nullptr;
    } shadow;
public:
    CXX_Queue_info_t() = default;

    // This MUST be called before dequeue below.
    bool empty() {
        if (shadow.head.load(std::memory_order_relaxed) == nullptr) {
            if (normal.head.load(std::memory_order_relaxed) == nullptr) return true;
            else {
                shadow.head.store(normal.head.load(std::memory_order_relaxed),std::memory_order_release);
                normal.head.store(nullptr,std::memory_order_relaxed);
            }
        }
        return false;
    }

    template<typename MESSAGE>
    void enqueue(MESSAGE *m) {
        m->next.store(nullptr,std::memory_order_relaxed);
        // We need to make sure that the above store happens before the exchange.
        // We use memory_order_release on the exchange for this, but on POWERPC
        // the eieio instruction is sufficient (release uses lwsync normally).
        // The idea is eieio is cheap, lwsync is ~100's of cycles and sync is ~1000's.
        // eieio: Orders I/O and store/store to memory (but separately)
        // lwsync: orders load/load, store/store and load/store to memory
        // sync: orders everything
#if defined(__GNUC__) && defined(__PPC__)
        __asm__ __volatile__  ( "eieio"  ::: "memory" );
#else
        std::atomic_thread_fence(std::memory_order_release);
#endif
        auto prev = normal.tail.exchange(m,std::memory_order_relaxed);
        if (prev==nullptr) normal.head.store(m,std::memory_order_relaxed);
        else prev->next.store(m,std::memory_order_relaxed);
    }
    template<typename MESSAGE>
    void enqueue(MESSAGE &m) { enqueue(&m); }

    // You must call empty() before this function
    template<typename MESSAGE>
    MESSAGE &dequeue() {
        auto e = shadow.head.load(std::memory_order_relaxed);
        auto n = e->next.load(std::memory_order_relaxed);
        if (n != nullptr) shadow.head.store(n,std::memory_order_relaxed);
        else {
            shadow.head.store(nullptr,std::memory_order_relaxed);
            CXX_Queue_element_hdr_t *expected = e;
            if (!normal.tail.compare_exchange_strong(expected,nullptr,std::memory_order_relaxed)) {
                while ((n=e->next.load(std::memory_order_relaxed)) == nullptr) {
#if defined(HAVE_SCHED_YIELD)
                    sched_yield();
#elif defined(_MSC_VER)
                    SwitchToThread();
#endif
                }
                shadow.head.store(n,std::memory_order_relaxed);
            }
        }
        e->next.store(nullptr,std::memory_order_relaxed);
        std::atomic_thread_fence(std::memory_order_acquire);
        return *static_cast<MESSAGE *>(e);
    }
};

// IMPORTANT: Message queues internally queue only a basicMessage and not the
// derived type. When the message is dequeued, it is dynamically cast correctly.
// This allows messages to be put on more or less specific queues as required
// and the "replyQueue" can be generic. For example,
//   ewaldMessage -> cudaMessage -> basicMessage
// An ewaldMessage is sent to the CUDA thread which accepts a cudaMessage.
// When returned, it goes onto a queue that accepts an ewaldMessage.
struct basicMessage : public CXX_Queue_element_hdr_t {
    CXX_Queue_info_t *replyQueue;
    basicMessage() : replyQueue(0) {}
    virtual ~basicMessage() {} // This class needs to be polymorphic (so dynamic_cast works)
    void sendBack() {if (replyQueue) replyQueue->enqueue(this);}
};

struct basicQueue : public CXX_Queue_info_t {
public:
    basicMessage &dequeue() { return CXX_Queue_info_t::dequeue<basicMessage>(); }

    void enqueue(basicMessage *m) { CXX_Queue_info_t::enqueue(m); }
    void enqueue(basicMessage &m) { CXX_Queue_info_t::enqueue(&m); }
    void enqueue(const basicMessage &C,basicQueue &Q) {
        // We do modify "M", but we are done before we return. Promise.
        // This allows patterns like Q.enqueueAndWait(MessageType(...))
        basicMessage &M = const_cast<basicMessage &>(C);
        M.replyQueue = &Q;
        enqueue(M);
    }
};

template<class messageType>
struct messageQueue : public basicQueue {

    messageType &wait() {
        while (empty()) {
            // This is important in the case where we have oversubscribed the CPU
#if defined(HAVE_SCHED_YIELD)
            sched_yield();
#elif defined(_MSC_VER)
            SwitchToThread();
#endif
        }
        return dequeue();
    }
    messageType &dequeue() {
        basicMessage &M = basicQueue::dequeue();
        return dynamic_cast<messageType &>(M);
    }
};
} // namespace mdl
#endif
