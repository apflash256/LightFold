#pragma once

#include <math/vecmath.h>

#include <atomic>
#include <mutex>
#include <functional>

namespace lightfold {

    inline uint32_t FloatToBits(float f) {
        uint32_t ui;
        memcpy(&ui, &f, sizeof(float));
        return ui;
    }
    inline float BitsToFloat(uint32_t ui) {
        float f;
        memcpy(&f, &ui, sizeof(uint32_t));
        return f;
    }

    class AtomicFloat {
    public:
        // AtomicFloat Public Methods
        explicit AtomicFloat(float v = 0) { bits = FloatToBits(v); }
        operator float() const { return BitsToFloat(bits); }
        float operator=(float v) { bits = FloatToBits(v); return v; }
        void Add(float v) {
            uint32_t oldBits = bits, newBits;
            do {
                newBits = FloatToBits(BitsToFloat(oldBits) + v);
            } while (!bits.compare_exchange_weak(oldBits, newBits));
        }

    private:
        // AtomicFloat Private Data
        std::atomic<uint32_t> bits;
    };

    class Barrier {
    public:
        Barrier(int count) : count(count) { }
        ~Barrier() { }
        void Wait();

    private:
        std::mutex mutex;
        std::condition_variable cv;
        int count;
    };

    void ParallelFor(std::function<void(int64_t)> func, int64_t count,
        int chunkSize = 1);
    extern thread_local int ThreadIndex;
    void ParallelFor2D(std::function<void(Point2i)> func, const Point2i& count);
    int MaxThreadIndex();
    int NumSystemCores();

    void ParallelInit();
    void ParallelCleanup();
    void MergeWorkerThreadStats();

} // namespace lightfold