#ifndef BASICMESSAGE_H
#define BASICMESSAGE_H

#include "opa_queue.h"

namespace mdl {

struct basicMessage {
    OPA_Queue_element_hdr_t hdr;
    struct OPA_Queue_info_t *replyQueue;
    basicMessage() : replyQueue(0) {OPA_Queue_header_init(&hdr);}
    virtual ~basicMessage() {} // This needs to be polymorphic
    void sendBack() {if (replyQueue) OPA_Queue_enqueue(replyQueue, this, basicMessage, hdr);}
    };

template<class messageType>
struct basicQueue : public OPA_Queue_info_t {
    basicQueue() { OPA_Queue_init(this); }
    bool empty() {return OPA_Queue_is_empty(this);}

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
    void enqueue(messageType *M) { basicMessage *m = M; OPA_Queue_enqueue(this,  M, basicMessage, hdr); }
    void enqueue(messageType &M) { enqueue(&M); }
    void enqueue(const messageType &C,basicQueue<messageType> &Q, bool bWait=false) {
	// We do modify "M", but we are done before we return. Promise.
	messageType &M = const_cast<messageType&>(C);
	M.replyQueue = &Q;
	enqueue(M);
//	OPA_Queue_enqueue(this, &M, messageType, hdr);
	if (bWait) Q.wait();
	}
    void enqueueAndWait(const messageType &M){
	basicQueue<messageType> wait;
	enqueue(M,wait,true);
	}
    messageType &dequeue() {
	basicMessage *M;
	OPA_Queue_dequeue(this, M, basicMessage, hdr);
	return * dynamic_cast<messageType*>(M);
	}
    };
} // namespace mdl
#endif
