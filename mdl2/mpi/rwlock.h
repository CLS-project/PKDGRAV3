#ifndef RWLOCK_HINCLUDED
#define RWLOCK_HINCLUDED
#include <atomic>
#include <thread>

namespace mdl {

// This is an implementation of a single writer, multiple reader lock.
class RWLock {
protected:
    std::atomic<int> count;
public:
    RWLock() : count(0) {}
    void lock_read() {
    	while(true) {
	    int prev = count;
	    if (prev>=0 && count.compare_exchange_weak(prev,prev+1)) return;
	    std::this_thread::yield();
	    }
	}
    void unlock_read() {
    	--count;
	}
    void lock_write() {
    	while(true) {
	    int prev = count;
	    if (prev<=0 && count.compare_exchange_weak(prev,prev-1)) return;
	    std::this_thread::yield();
	    }
	}
    void unlock_write() {
    	++count;
	}
    };

class read_lock {
protected:
    RWLock &lock;
public:
    read_lock(RWLock &lock) : lock(lock) {
    	lock.lock_read();
	}
    ~read_lock() {
    	lock.unlock_read();
	}
    };

class write_lock {
protected:
    RWLock &lock;
public:
    write_lock(RWLock &lock) : lock(lock) {
    	lock.lock_write();
	}
    ~write_lock() {
    	lock.unlock_write();
	}
    };

} // namespace mdl
#endif