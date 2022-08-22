#ifndef BASICMESSAGE_H
#define BASICMESSAGE_H

#ifdef HAVE_SCHED_YIELD
    #include "sched.h"
#endif
#include <atomic>
#include <new>

namespace mdl {

#ifdef __cpp_lib_hardware_interference_size
    using std::hardware_constructive_interference_size;
    using std::hardware_destructive_interference_size;
#else
    // 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned │ ...
    constexpr std::size_t hardware_constructive_interference_size = 64;
    constexpr std::size_t hardware_destructive_interference_size = 64;
#endif

class CXX_Queue_element_hdr_t {
protected:
    friend class CXX_Queue_info_t;
    std::atomic<CXX_Queue_element_hdr_t *> next = nullptr;
public:
    CXX_Queue_element_hdr_t() = default;
    virtual ~CXX_Queue_element_hdr_t() = default;
};

class CXX_Queue_info_t {
    struct alignas(hardware_destructive_interference_size) normal_ {
        std::atomic<CXX_Queue_element_hdr_t *> head = nullptr;
        std::atomic<CXX_Queue_element_hdr_t *> tail = nullptr;
    } normal;
    struct alignas(hardware_destructive_interference_size) shadow_ {
        std::atomic<CXX_Queue_element_hdr_t *> head = nullptr;
    } shadow;
public:
    CXX_Queue_info_t() = default;

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
        auto prev = normal.tail.exchange(m);
        if (prev==nullptr) normal.head.store(m,std::memory_order_relaxed);
        else prev->next.store(m,std::memory_order_relaxed);

    }
    template<typename MESSAGE>
    void enqueue(MESSAGE &m) { enqueue(&m); }

    template<typename MESSAGE>
    MESSAGE &dequeue() {
        auto e = shadow.head.load(std::memory_order_relaxed);
        if (e->next.load(std::memory_order_relaxed) != nullptr) shadow.head.store(e->next.load(std::memory_order_relaxed),std::memory_order_relaxed);
        else {
            shadow.head.store(nullptr,std::memory_order_relaxed);
            CXX_Queue_element_hdr_t *expected = e;
            if (!normal.tail.compare_exchange_strong(expected,nullptr)) {
                while (e->next.load(std::memory_order_relaxed) == nullptr) {
#if defined(HAVE_SCHED_YIELD)
                    sched_yield();
#elif defined(_MSC_VER)
                    SwitchToThread();
#endif
                }
                shadow.head.store(e->next.load(std::memory_order_relaxed),std::memory_order_relaxed);
            }
        }
        e->next.store(nullptr,std::memory_order_release);
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
