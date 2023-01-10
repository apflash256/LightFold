#pragma once

#include <math/vecmath.h>
#include <objects/camera.h>
#include <rng/pcg_random.h>

#include <memory>
#include <vector>
#include <random>

namespace lightfold {

    static std::uniform_real_distribution<float> UniformFloat(0.f, 1.f);

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

} // namespace lightfold