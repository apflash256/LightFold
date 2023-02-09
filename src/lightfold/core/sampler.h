#pragma once

#include <core/camera.h>
#include <math/vecmath.h>

#include <rng/pcg_random.h>
#include <memory>
#include <vector>
#include <random>

namespace lightfold {

    static std::uniform_real_distribution<float> UniformFloat(0.f, 1.f);
    using UniformUInt32 = std::uniform_int_distribution<uint32_t>;

    static const float OneMinusEpsilon = 0x1.fffffep-1;

    class Sampler {
    public:
        // Sampler Public Methods
        Sampler(int64_t samplesPerPixel) : samplesPerPixel(samplesPerPixel) { }
        virtual void StartPixel(const Point2i& p);
        virtual float Get1D() = 0;
        virtual Point2f Get2D() = 0;
        CameraSample GetCameraSample(const Point2i& pRaster);
        void Request1DArray(int n);
        void Request2DArray(int n);
        virtual int RoundCount(int n) const { return n; }
        const float* Get1DArray(int n);
        const Point2f* Get2DArray(int n);
        virtual bool StartNextSample();
        virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
        virtual bool SetSampleNumber(int64_t sampleNum);

        // Sampler Public Data
        const int64_t samplesPerPixel;

    protected:
        // Sampler Protected Data
        Point2i currentPixel;
        int64_t currentPixelSampleIndex;
        std::vector<int> samples1DArraySizes, samples2DArraySizes;
        std::vector<std::vector<float>> sampleArray1D;
        std::vector<std::vector<Point2f>> sampleArray2D;

    private:
        // Sampler Private Data
        size_t array1DOffset, array2DOffset;
    };

    class PixelSampler : public Sampler {
    public:
        // PixelSampler Public Methods
        PixelSampler(int64_t samplesPerPixel, int nSampledDimensions) : Sampler(samplesPerPixel) {
            for (int i = 0; i < nSampledDimensions; ++i) {
                samples1D.push_back(std::vector<float>(samplesPerPixel));
                samples2D.push_back(std::vector<Point2f>(samplesPerPixel));
            }
        }
        bool StartNextSample();
        bool SetSampleNumber(int64_t);
        float Get1D();
        Point2f Get2D();

    protected:
        // PixelSampler Protected Data
        std::vector<std::vector<float>> samples1D;
        std::vector<std::vector<Point2f>> samples2D;
        int current1DDimension = 0, current2DDimension = 0;
        pcg32 rng;
    };

    class GlobalSampler : public Sampler {
    public:
        // GlobalSampler Public Methods
        bool StartNextSample();
        void StartPixel(const Point2i&);
        bool SetSampleNumber(int64_t sampleNum);
        float Get1D();
        Point2f Get2D();
        GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) { }
        virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
        virtual float SampleDimension(int64_t index, int dimension) const = 0;

    private:
        // GlobalSampler Private Data
        int dimension;
        int64_t intervalSampleIndex;
        static const int arrayStartDim = 5;
        int arrayEndDim;
    };

    template <typename T>
    inline void Shuffle(T* samp, int count, int nDimensions, pcg32& rng) {
        for (int i = 0; i < count; ++i) {
            int other = i + UniformUInt32(0, count - i - 1)(rng);
            for (int j = 0; j < nDimensions; ++j)
                std::swap(samp[nDimensions * i + j],
                    samp[nDimensions * other + j]);
        }
    }

    inline float RadicalInverse(int baseIndex, uint64_t a) {
        unsigned int base = Primes[baseIndex];
        uint64_t limit = ~0ull / base - base;
        float invBase = 1.f / base, invBaseM = 1;
        uint64_t reversedDigits = 0;
        while (a && reversedDigits < limit) {
            uint64_t next = a / base;
            uint64_t digit = a - next * base;
            reversedDigits = reversedDigits * base + digit;
            invBaseM *= invBase;
            a = next;
        }
        return std::min(reversedDigits * invBaseM, OneMinusEpsilon);
    }

    inline uint64_t InverseRadicalInverse(uint64_t inverse, int base, int nDigits) {
        uint64_t index = 0;
        for (int i = 0; i < nDigits; ++i) {
            uint64_t digit = inverse % base;
            inverse /= base;
            index = index * base + digit;
        }
        return index;
    }

    std::vector<uint16_t> ComputeRadicalInversePermutations(pcg32& rng);

    inline float ScrambledRadicalInverse(int baseIndex, uint64_t a, const uint16_t* perm) {
        unsigned int base = Primes[baseIndex];
        uint64_t limit = ~0ull / base - base;
        float invBase = 1.f / base, invBaseM = 1;
        uint64_t reversedDigits = 0;
        while (1 - (base - 1) * invBaseM < 1 && reversedDigits < limit) {
            uint64_t next = a / base;
            int digitValue = a - next * base;
            reversedDigits = reversedDigits * base + perm[digitValue];
            invBaseM *= invBase;
            a = next;
        }
        return std::min(invBaseM * reversedDigits, OneMinusEpsilon);
    }

    static constexpr int kMaxResolution = 128;

    class HaltonSampler : public GlobalSampler {
    public:
        // HaltonSampler Public Methods
        HaltonSampler(int samplesPerPixel, const Bounds2i& sampleBounds);
        int64_t GetIndexForSample(int64_t sampleNum) const;
        float SampleDimension(int64_t index, int dimension) const;
        std::unique_ptr<Sampler> Clone(int seed);

    private:
        // HaltonSampler Private Data
        static std::vector<uint16_t> radicalInversePermutations;
        Point2i baseScales, baseExponents;
        int sampleStride;
        int multInverse[2];
        mutable Point2i pixelForOffset = Point2i(std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
        mutable int64_t offsetForCurrentPixel;

        // HaltonSampler Private Methods
        const uint16_t* PermutationForDimension(int dim) const {
            return &radicalInversePermutations[PrimeSums[dim]];
        }
        static uint64_t multiplicativeInverse(int64_t a, int64_t n) {
            int64_t x, y;
            extendedGCD(a, n, &x, &y);
            return x % n;
        }
        static void extendedGCD(uint64_t a, uint64_t b, int64_t* x, int64_t* y) {
            if (b == 0) {
                *x = 1;
                *y = 0;
                return;
            }
            int64_t d = a / b, xp, yp;
            extendedGCD(b, a % b, &xp, &yp);
            *x = yp;
            *y = xp - (d * yp);
        }
    };

    template <typename Predicate>
    inline int FindInterval(int size, const Predicate& pred) {
        int first = 0, len = size;
        while (len > 0) {
            int half = len >> 1, middle = first + half;
            // Bisect range based on value of _pred_ at _middle_
            if (pred(middle)) {
                first = middle + 1;
                len -= half + 1;
            }
            else
                len = half;
        }
        return Clamp(first - 1, 0, size - 2);
    }

    struct Distribution1D {
        // Distribution1D Public Methods
        Distribution1D(const float* f, int n) : func(f, f + n), cdf(n + 1) {
            // Compute integral of step function at $x_i$
            cdf[0] = 0;
            for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;

            // Transform step function integral into CDF
            funcInt = cdf[n];
            if (funcInt == 0) {
                for (int i = 1; i < n + 1; ++i) cdf[i] = float(i) / float(n);
            }
            else {
                for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
            }
        }
        int Count() const { return (int)func.size(); }
        float SampleContinuous(float u, float* pdf, int* off = nullptr) const {
            // Find surrounding CDF segments and _offset_
            int offset = FindInterval((int)cdf.size(),
                [&](int index) { return cdf[index] <= u; });
            if (off)
                *off = offset;
            // Compute offset along CDF segment
            float du = u - cdf[offset];
            if ((cdf[offset + 1] - cdf[offset]) > 0)
                du /= (cdf[offset + 1] - cdf[offset]);

            // Compute PDF for sampled offset
            if (pdf)
                *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;

            // Return $x\in{}[0,1)$ corresponding to sample
            return (offset + du) / Count();
        }
        int SampleDiscrete(float u, float* pdf = nullptr,
            float* uRemapped = nullptr) const {
            // Find surrounding CDF segments and _offset_
            int offset = FindInterval((int)cdf.size(),
                [&](int index) { return cdf[index] <= u; });
            if (pdf)
                *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
            if (uRemapped)
                *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
            return offset;
        }
        float DiscretePDF(int index) const {
            return func[index] / (funcInt * Count());
        }

        // Distribution1D Public Data
        std::vector<float> func, cdf;
        float funcInt;
    };

} // namespace lightfold