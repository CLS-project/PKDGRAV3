#ifndef MDLGPU_H
#define MDLGPU_H
#include "basicmessage.h"

namespace mdl {
namespace gpu {

class Message : public basicMessage {
public:
    virtual void finish() = 0;
    Message() = default;
    virtual ~Message() = default;
};
struct MessageQueue : public messageQueue<Message> {};

class Client {
protected:
    int nGPU;
    MessageQueue done;
public:
    Client() : nGPU(0) {}
    int flushCompleted(); // Returns how many are still outstanding
    operator basicQueue &() {++nGPU; return done;}
};


} // namespace gpu
} // namespace mdl
#endif
