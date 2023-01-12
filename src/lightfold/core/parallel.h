#include <atomic>

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

} // namespace lightfold