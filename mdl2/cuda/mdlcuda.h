#ifndef MDLCUDA_H
#define MDLCUDA_H

#include "basicmessage.h"
#include "cuda.h"
#include "cuda_runtime_api.h"

#include <vector>
#include <list>
#include <atomic>

namespace mdl {

class Device;
class cudaMessage : public basicMessage {
    friend class Device; // So we can launch()
    int iDevice; // Device number to use or -1 for any (the normal case)
protected:
    virtual void launch(Device &device,cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) = 0;
public:
    virtual void finish() {}
    cudaMessage(int iDevice=-1) : iDevice(iDevice) {}
    virtual ~cudaMessage() {}
    int getDevice() const {return iDevice;}
};

struct cudaMessageQueue : public messageQueue<cudaMessage> {};

class Stream : public basicMessage {
    friend class Device;
protected:
    cudaStream_t stream; // CUDA execution stream
    class Device *device; // Device associated with this stream
    cudaMessage *message; // Message currently in progress (or NULL)

    // This should be handled in a better way, but for now let the cheese happen.
    const int requestBufferSize = 2*1024*1024;
    const int resultsBufferSize = 2*1024*1024;
    void *pCudaBufIn, *pCudaBufOut;

public:
    explicit Stream(class Device *device);
    ~Stream();
    cudaStream_t getStream();
};
class StreamQueue : public messageQueue<class Stream> {};

class Work {
public:
    Work();
};

class Device {
    friend class Stream;
    friend class CUDA;
    int DevAttrMaxBlocksPerMultiprocessor,DevAttrMaxThreadsPerMultiprocessor,DevAttrWarpSize;
    int DevAttrSingleToDoublePrecisionPerfRatio;
public:
    auto WarpSize() const {return DevAttrWarpSize;}
    auto MaxBlocksPerMultiprocessor() const {return DevAttrMaxBlocksPerMultiprocessor;}
    auto MaxThreadsPerMultiprocessor() const {return DevAttrMaxThreadsPerMultiprocessor;}
    auto SingleToDoublePrecisionPerfRatio() const {return DevAttrSingleToDoublePrecisionPerfRatio;}
protected:
    static void CUDART_CB kernel_finished( void  *userData );
    void kernel_finished( Stream *stream );
    StreamQueue free_streams;
    std::atomic_int busy_streams;
    int iDevice, nStreams;
    static bool compareBusy(const Device &lhs,const Device &rhs) {return lhs.busy_streams < rhs.busy_streams; }
public:
    explicit Device(int iDevice,int nStreams);
    bool empty() {return free_streams.empty();}
    void launch(cudaMessage &M);
};

class CUDA : public cudaMessageQueue {
protected:
    std::list<Device> devices;
    Device **pDevice;
    int nDevices;
public:
    CUDA();
    ~CUDA();
    void initialize(int nStreamsPerDevice=8);
    auto numDevices() const { return devices.size(); }
    bool isActive()   const { return numDevices() > 0; }
    void initiate();
};

} // namespace mdl
#endif
