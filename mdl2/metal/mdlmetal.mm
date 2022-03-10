#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION

#include "mdlmetal.h"
using namespace mdl::metal;

METAL::METAL() {

}

METAL::~METAL() {

}

void METAL::initialize(int nStreamsPerDevice) {
    //auto main_queue = dispatch_get_main_queue();
    //NS::Error *error;
    auto deviceList = MTL::CopyAllDevices();
    for(auto iDevice = 0; iDevice<deviceList->count(); ++iDevice) {
        auto pDevice = reinterpret_cast<MTL::Device*>(deviceList->object(iDevice));
        if (!pDevice->lowPower())
            devices.emplace_back(this,pDevice,nStreamsPerDevice);
    }
	deviceList->release();
}

void METAL::launch() {
    while (!empty()) { // A message is waiting. Find a stream if we can.
        assert(devices.size()>0); // This would be odd at this point. No progress could be made.
        // Find the device with the fewest busy streams (most idle streams)
        auto device = std::min_element(devices.begin(),devices.end(),Device::compareBusy);
        if (device != devices.end() && !device->empty()) {
            metalMessage &M = dequeue();
            device->launch(M); // Launch message M on the most idle device.
        }
        else break; // No free streams to launch work at this time
    }
}

dispatch_data_t METAL::getLibrary(const unsigned char *data,int length) {
    auto lib = libraries.find(data); // Has the dispatch data been cached?
    if (lib != libraries.end()) return lib->second;
    else {
        auto main_queue = dispatch_get_main_queue();
        auto dispatch_data = dispatch_data_create(data,length,main_queue,DISPATCH_DATA_DESTRUCTOR_DEFAULT);
        libraries.emplace(data,dispatch_data);
        return dispatch_data;
        // return dispatch_data_empty;
    }
}
MTL::Library *Device::getLibrary(const unsigned char *data,int length) {
    auto lib = libraries.find(data);
    if (lib != libraries.end()) return lib->second;
    else {
        NS::Error *error;
        auto dlib = pDevice->newLibrary(metal->getLibrary(data,length),&error);
        if (error->code()!=0) {
            fprintf(stderr,"ERROR: cannot create METAL library, error=%ld\n",error->code());
            assert(error->code()==0);
            abort();
        }
        libraries.emplace(data,dlib);
        return dlib;
    }
}
MTL::Function *Device::getFunction(const unsigned char *data, int length, const char *name) {
    // If we have compiled the function for this device, then return the MTL::Function* reference
    auto ii = functions.find(name);
    if (ii != functions.end()) return ii->second;
    // Otherwise we need to create it (and possibly the library as well)
    else {
        auto f = getLibrary(data,length)->newFunction(NS::String::string(name,NS::ASCIIStringEncoding));
        functions.emplace(name,f);
        return f;
    }
}

MTL::Buffer *Device::getBuffer(const void* pointer, NS::UInteger length) {
    auto ii = buffers.find(pointer);
    if (ii != buffers.end()) return ii->second;
    else {
        auto b = pDevice->newBuffer(pointer,length,MTL::ResourceStorageModeManaged,nullptr);
        buffers.emplace(pointer,b);
        return b;
    }
}

Device::Device(METAL *metal,MTL::Device *pDevice,int nStreams)
    : metal(metal), pDevice(pDevice->retain()), nStreams(nStreams)
{
    for (auto i=0; i<nStreams; ++i) {
        free_streams.enqueue(new Stream(*this));
    }
}

Device::~Device() {
    pDevice->release();
}

void Device::launch(metalMessage &M) {
    if (free_streams.empty()) abort();
    auto &stream = free_streams.dequeue();
    stream.message = &M; // Save the message (for kernel_finished)
    ++busy_streams; // This is atomic
    stream.launch(M);
}

Stream::Stream(class Device &device)
    : device(device), queue(device.newCommandQueue()), message(nullptr), finished(this)
{
}

Stream::~Stream() {
    queue->release();
}

MTL::ComputePipelineState *Stream::getPipeline(const unsigned char *data, int length, const char *name) {
    // If we have compiled the function for this device, then return the MTL::Function* reference
    auto ii = pipelines.find(name);
    if (ii != pipelines.end()) return ii->second;
    // Otherwise we need to create it (and possibly the library as well)
    else {
        NS::Error *error;
        auto f = device.getFunction(data,length,name);
	    auto p = device.pDevice->newComputePipelineState(f,&error);
        if (error->code()!=0) {
            fprintf(stderr,"ERROR: cannot create Pipeline state for %s, error=%ld\n",name,error->code());
            assert(error->code()==0);
            abort();
        }
        pipelines.emplace(name,p);
        return p;
    }
}

void Stream::launch(metalMessage &M) {
	auto cbuf = queue->commandBuffer();
    M.launch(*this,cbuf);
	cbuf->addCompletedHandler(finished);
	cbuf->commit();
}

void Stream::kernel_finished(MTL::CommandBuffer*cb) {
    message->sendBack();
    message = NULL;
    --device.busy_streams; // This is atomic
    device.free_streams.enqueue(this);
}
