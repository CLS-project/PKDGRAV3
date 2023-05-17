#ifndef MDLMETAL_H
#define MDLMETAL_H
#include "Metal.hpp"
#include "gpu/mdlgpu.h"

#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <atomic>

namespace mdl {
namespace metal {
class Device;
class Stream;

class metalMessage : public gpu::Message {
    friend class Stream; // So we can launch()
protected:
    virtual void launch(Stream &stream,MTL::CommandBuffer *cbuf) = 0;
public:
    virtual void finish() {}
    metalMessage() {}
    virtual ~metalMessage() {}
};
struct metalMessageQueue : public messageQueue<metalMessage> {};

class Stream : public basicMessage {
    friend class Device;
    friend class run_on_complete;
protected:
    class Device &device;  // Device associated with this stream
    MTL::CommandQueue *queue;
    metalMessage *message; // Message currently in progress (or NULL)
    std::unordered_map<std::string,MTL::ComputePipelineState *> pipelines;

    // This should be handled in a better way, but for now let the cheese happen.
    const int requestBufferSize = 2*1024*1024;
    const int resultsBufferSize = 2*1024*1024;
    void kernel_finished(MTL::CommandBuffer *cb);
    struct run_on_complete {
        Stream *stream;
        void operator()(MTL::CommandBuffer *cb) {
            stream->kernel_finished(cb);
        }
        run_on_complete() = delete;
        run_on_complete(Stream *stream) : stream(stream) {}
    } finished;
    void launch(metalMessage &M);
public:
    Device &getDevice() {return device;}
    explicit Stream(class Device &device);
    MTL::ComputePipelineState *getPipeline(const unsigned char *library_data, int length, const char *name);
    ~Stream();
};
class StreamQueue : public messageQueue<class Stream> {};

class Work {
public:
    Work();
};

class Device {
    friend class Stream;
    friend class METAL;
protected:
    class METAL *metal;
    static void kernel_finished( void  *userData );
    void kernel_finished( Stream *stream );
    StreamQueue free_streams;
    std::atomic_int busy_streams;
    MTL::Device *pDevice;
    int nStreams;
    MTL::CommandQueue *newCommandQueue() { return pDevice->newCommandQueue(); }
    static bool compareBusy(const Device &lhs,const Device &rhs) {return lhs.busy_streams < rhs.busy_streams; }
    std::unordered_map<const unsigned char *,MTL::Library *> libraries;
    std::unordered_map<std::string,MTL::Function *> functions;
    std::unordered_map<const void *,MTL::Buffer *> buffers;
public:
    explicit Device(METAL *metal,MTL::Device *pDevice,int nStreams);
    virtual ~Device();
    bool empty() {return free_streams.empty();}
    void launch(metalMessage &M);
    MTL::Library *getLibrary(const unsigned char *library_data,int length);
    MTL::Function *getFunction(const unsigned char *library_data, int length, const char *name);
    MTL::Buffer *getBuffer(const void *pointer, NS::UInteger length);
    auto maxThreadgroupMemoryLength() { return pDevice->maxThreadgroupMemoryLength(); }
};

class METAL : public metalMessageQueue {
protected:
    std::list<Device> devices;
    std::unordered_map<const unsigned char *,dispatch_data_t> libraries;
    int nDevices;
public:
    METAL();
    virtual ~METAL();
    void initialize(int nStreamsPerDevice=8);
    bool isActive() { return devices.size() > 0; }
    void launch();
    dispatch_data_t getLibrary(const unsigned char *library_data,int length);
};
} // namespace metal
} // namespace mdl
#endif
