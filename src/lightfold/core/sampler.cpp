#include <core/sampler.h>
#include <core/sobolmatrices.h>

#include <cassert>

namespace lightfold {

    CameraSample Sampler::GetCameraSample(const Point2i& pRaster) {
        CameraSample cs;
        cs.pFilm = (Point2f)(pRaster + (Vector2f)Get2D());
        cs.pLens = Get2D();
        return cs;
    }

    void Sampler::StartPixel(const Point2i& p) {
        currentPixel = p;
        currentPixelSampleIndex = 0;
        array1DOffset = array2DOffset = 0;
    }

    bool Sampler::StartNextSample() {
        array1DOffset = array2DOffset = 0;
        return ++currentPixelSampleIndex < samplesPerPixel;
    }

    bool Sampler::SetSampleNumber(int64_t sampleNum) {
        array1DOffset = array2DOffset = 0;
        currentPixelSampleIndex = sampleNum;
        return currentPixelSampleIndex < samplesPerPixel;
    }

    void Sampler::Request1DArray(int n) {
        samples1DArraySizes.push_back(n);
        sampleArray1D.push_back(std::vector<float>(n * samplesPerPixel));
    }

    void Sampler::Request2DArray(int n) {
        samples2DArraySizes.push_back(n);
        sampleArray2D.push_back(std::vector<Point2f>(n * samplesPerPixel));
    }

    const float* Sampler::Get1DArray(int n) {
        if (array1DOffset == sampleArray1D.size())
            return nullptr;
        return &sampleArray1D[array1DOffset++][currentPixelSampleIndex * n];
    }

    const Point2f* Sampler::Get2DArray(int n) {
        if (array2DOffset == sampleArray2D.size())
            return nullptr;
        return &sampleArray2D[array2DOffset++][currentPixelSampleIndex * n];
    }

    bool PixelSampler::StartNextSample() {
        current1DDimension = current2DDimension = 0;
        return Sampler::StartNextSample();
    }

    bool PixelSampler::SetSampleNumber(int64_t sampleNum) {
        current1DDimension = current2DDimension = 0;
        return Sampler::SetSampleNumber(sampleNum);
    }

    float PixelSampler::Get1D() {
        if (current1DDimension < samples1D.size())
            return samples1D[current1DDimension++][currentPixelSampleIndex];
        else
            return UniformFloat(rng);
    }

    Point2f PixelSampler::Get2D() {
        if (current2DDimension < samples2D.size())
            return samples2D[current2DDimension++][currentPixelSampleIndex];
        else
            return Point2f(UniformFloat(rng), UniformFloat(rng));
    }

    void GlobalSampler::StartPixel(const Point2i& p) {
        Sampler::StartPixel(p);
        dimension = 0;
        intervalSampleIndex = GetIndexForSample(0);
        arrayEndDim = arrayStartDim + sampleArray1D.size() + 2 * sampleArray2D.size();
        for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
            int nSamples = samples1DArraySizes[i] * samplesPerPixel;
            for (int j = 0; j < nSamples; ++j) {
                int64_t index = GetIndexForSample(j);
                sampleArray1D[i][j] =
                    SampleDimension(index, arrayStartDim + i);
            }
        }
        int dim = arrayStartDim + samples1DArraySizes.size();
        for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
            int nSamples = samples2DArraySizes[i] * samplesPerPixel;
            for (int j = 0; j < nSamples; ++j) {
                int64_t idx = GetIndexForSample(j);
                sampleArray2D[i][j].x = SampleDimension(idx, dim);
                sampleArray2D[i][j].y = SampleDimension(idx, dim + 1);
            }
            dim += 2;
        }
        assert(dim == arrayEndDim);
    }

    bool GlobalSampler::StartNextSample() {
        dimension = 0;
        intervalSampleIndex = GetIndexForSample(currentPixelSampleIndex + 1);
        return Sampler::StartNextSample();
    }

    bool GlobalSampler::SetSampleNumber(int64_t sampleNum) {
        dimension = 0;
        intervalSampleIndex = GetIndexForSample(sampleNum);
        return Sampler::SetSampleNumber(sampleNum);
    }

    float GlobalSampler::Get1D() {
        if (dimension >= arrayStartDim && dimension < arrayEndDim)
            dimension = arrayEndDim;
        return SampleDimension(intervalSampleIndex, dimension++);
    }

    Point2f GlobalSampler::Get2D() {
        if (dimension + 1 >= arrayStartDim && dimension < arrayEndDim)
            dimension = arrayEndDim;
        Point2f p(SampleDimension(intervalSampleIndex, dimension),
            SampleDimension(intervalSampleIndex, dimension + 1));
        dimension += 2;
        return p;
    }

    std::vector<uint16_t> ComputeRadicalInversePermutations(pcg32& rng) {
        std::vector<uint16_t> perms;
        int permArraySize = 0;
        for (int i = 0; i < PrimeTableSize; ++i)
            permArraySize += Primes[i];
        perms.resize(permArraySize);
        uint16_t* p = &perms[0];
        for (int i = 0; i < PrimeTableSize; ++i) {
            for (int j = 0; j < Primes[i]; ++j)
                p[j] = j;
            Shuffle(p, Primes[i], 1, rng);
            p += Primes[i];
        }
        return perms;
    }

    std::vector<uint16_t> HaltonSampler::radicalInversePermutations;

    HaltonSampler::HaltonSampler(int samplesPerPixel, const Bounds2i& sampleBounds)
        : GlobalSampler(samplesPerPixel) {
        if (radicalInversePermutations.size() == 0) {
            pcg32 rng;
            radicalInversePermutations = ComputeRadicalInversePermutations(rng);
        }
        Vector2i res = sampleBounds.pMax - sampleBounds.pMin;
        for (int i = 0; i < 2; ++i) {
            int base = (i == 0) ? 2 : 3;
            int scale = 1, exp = 0;
            while (scale < std::min(res[i], kMaxResolution)) {
                scale *= base;
                ++exp;
            }
            baseScales[i] = scale;
            baseExponents[i] = exp;
        }
        sampleStride = baseScales[0] * baseScales[1];
        multInverse[0] = multiplicativeInverse(baseScales[1], baseScales[0]);
        multInverse[1] = multiplicativeInverse(baseScales[0], baseScales[1]);
    }

    int64_t HaltonSampler::GetIndexForSample(int64_t sampleNum) const {
        if (currentPixel != pixelForOffset) {
            offsetForCurrentPixel = 0;
            if (sampleStride > 1) {
                Point2i pm(currentPixel[0] % kMaxResolution, currentPixel[1] % kMaxResolution);
                for (int i = 0; i < 2; ++i) {
                    uint64_t dimOffset = (i == 0) ?
                        InverseRadicalInverse(pm[i], 2, baseExponents[i]) :
                        InverseRadicalInverse(pm[i], 3, baseExponents[i]);
                    offsetForCurrentPixel += dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
                }
                offsetForCurrentPixel %= sampleStride;
            }
            pixelForOffset = currentPixel;
        }
        return offsetForCurrentPixel + sampleNum * sampleStride;
    }

    float HaltonSampler::SampleDimension(int64_t index, int dim) const {
        if (dim == 0)
            return RadicalInverse(dim, index >> baseExponents[0]);
        else if (dim == 1)
            return RadicalInverse(dim, index / baseScales[1]);
        else
            return ScrambledRadicalInverse(dim, index, PermutationForDimension(dim));
    }

    std::unique_ptr<Sampler> HaltonSampler::Clone(int seed) {
        return std::unique_ptr<Sampler>(new HaltonSampler(*this));
    }

    inline uint64_t SobolIntervalToIndex(const uint32_t m, uint64_t frame,
        const Point2i& p) {
        if (m == 0) return 0;

        const uint32_t m2 = m << 1;
        uint64_t index = uint64_t(frame) << m2;

        uint64_t delta = 0;
        for (int c = 0; frame; frame >>= 1, ++c)
            if (frame & 1)  // Add flipped column m + c + 1.
                delta ^= VdCSobolMatrices[m - 1][c];

        // flipped b
        uint64_t b = (((uint64_t)((uint32_t)p.x) << m) | ((uint32_t)p.y)) ^ delta;

        for (int c = 0; b; b >>= 1, ++c)
            if (b & 1)  // Add column 2 * m - c.
                index ^= VdCSobolMatricesInv[m - 1][c];

        return index;
    }

    static const float FloatOneMinusEpsilon = 0.99999994;

    inline float SobolSample(int64_t a, int dimension, uint32_t scramble = 0) {
        uint32_t v = scramble;
        for (int i = dimension * SobolMatrixSize; a != 0; a >>= 1, i++)
            if (a & 1) v ^= SobolMatrices32[i];
        return std::min(v * 2.3283064365386963e-10f /* 1/2^32 */,
            FloatOneMinusEpsilon);
    }

    int64_t SobolSampler::GetIndexForSample(int64_t sampleNum) const {
        return SobolIntervalToIndex(log2Resolution, sampleNum,
            Point2i(currentPixel - sampleBounds.pMin));
    }

    float SobolSampler::SampleDimension(int64_t index, int dim) const {
        float s = SobolSample(index, dim);
        // Remap Sobol$'$ dimensions used for pixel samples
        if (dim == 0 || dim == 1) {
            s = s * resolution + sampleBounds.pMin[dim];
            s = Clamp(s - currentPixel[dim], 0.f, OneMinusEpsilon);
        }
        return s;
    }

    std::unique_ptr<Sampler> SobolSampler::Clone(int seed) {
        return std::unique_ptr<Sampler>(new SobolSampler(*this));
    }

} // namespace lightfold