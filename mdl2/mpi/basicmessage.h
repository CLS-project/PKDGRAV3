#ifndef BASICMESSAGE_H
#define BASICMESSAGE_H

#include "opa_queue.h"

namespace mdl {

// IMPORTANT: Message queues internally queue only a basicMessage and not the
// derived type. When the message is dequeued, it is dynamically cast correctly.
// This allows messages to be put on more or less specific queues as required
// and the "replyQueue" can be generic. For example,
//   ewaldMessage -> cudaMessage -> basicMessage
// An ewaldMessage is sent to the CUDA thread which accepts a cudaMessage.
// When returned, it goes onto a queue that accepts an ewaldMessage.

struct basicMessage {
    OPA_Queue_element_hdr_t hdr;
    struct OPA_Queue_info_t *replyQueue;
    basicMessage() : replyQueue(0) {OPA_Queue_header_init(&hdr);}
    virtual ~basicMessage() {} // This class needs to be polymorphic (so dynamic_cast works)
    void sendBack() {if (replyQueue) OPA_Queue_enqueue(replyQueue, this, basicMessage, hdr);}
    };

struct basicQueue : public OPA_Queue_info_t {
    basicQueue() { OPA_Queue_init(this); }
    bool empty() {return OPA_Queue_is_empty(this);}

    void enqueue(basicMessage *m) { OPA_Queue_enqueue(this,  m, basicMessage, hdr); }
    void enqueue(basicMessage &m) { enqueue(&m); }
    void enqueue(const basicMessage &C,basicQueue &Q) {
	// We do modify "M", but we are done before we return. Promise.
	// This allows patterns like Q.enqueueAndWait(MessageType(...))
	basicMessage &M = const_cast<basicMessage&>(C);
	M.replyQueue = &Q;
	enqueue(M);
	}

    basicMessage &dequeue() {
	basicMessage *M;
	OPA_Queue_dequeue(this, M, basicMessage, hdr);
	return *M;
	}
    };

template<class messageType>
struct messageQueue : public basicQueue {

    messageType &wait() {
	while (OPA_Queue_is_empty(this)) {
	    // This is important in the case where we have oversubscribed the CPU
#ifdef _MSC_VER
	    SwitchToThread();
#else
	    sched_yield();
#endif
	    }
	return dequeue();
	}
    messageType &dequeue() {
	basicMessage &M = basicQueue::dequeue();
	return dynamic_cast<messageType&>(M);
	}
    };
} // namespace mdl
#endif
