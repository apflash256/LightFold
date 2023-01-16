#pragma once

#include <core/spectra.h>
#include <core/filter.h>
#include <core/parallel.h>
#include <math/vecmath.h>

#include <vector>
#include <memory>
#include <string>
#include <mutex>

namespace lightfold {

    struct FilmTilePixel {
        Spectrum contribSum = 0.f;
        float filterWeightSum = 0.f;
    };

    struct Pixel {
        float rgb[3] = { 0, 0, 0 };
        float filterWeightSum = 0;
        AtomicFloat splatRGB[3];
    };

    class FilmTile {
    public:
        // FilmTile Public Methods
        FilmTile(const Bounds2i& pixelBounds, const GeoVector2f& filterRadius,
            const float* filterTable, int filterTableSize)
            : pixelBounds(pixelBounds), filterRadius(filterRadius),
            invFilterRadius(1 / filterRadius.x, 1 / filterRadius.y),
            filterTable(filterTable), filterTableSize(filterTableSize) {
            pixels = std::vector<FilmTilePixel>(std::max(0, pixelBounds.Area()));
        }
        void AddSample(const Point2f& pFilm, const Spectrum& L, float sampleWeight = 1.) {
            Point2f pFilmDiscrete = pFilm - GeoVector2f(0.5f, 0.5f);
            Point2i p0 = (Point2i)Ceil(pFilmDiscrete - filterRadius);
            Point2i p1 = (Point2i)(Floor(pFilmDiscrete + filterRadius) + Point2i(1, 1));
            p0 = Max(p0, pixelBounds.pMin);
            p1 = Min(p1, pixelBounds.pMax);
            int* ifx = (int*)alloca((p1.x - p0.x) * sizeof(int));
            for (int x = p0.x; x < p1.x; ++x) {
                float fx = std::abs((x - pFilmDiscrete.x) *
                    invFilterRadius.x * filterTableSize);
                ifx[x - p0.x] = std::min((int)std::floor(fx), filterTableSize - 1);
            }
            int* ify = (int*)alloca((p1.y - p0.y) * sizeof(int));
            for (int y = p0.y; y < p1.y; ++y) {
                float fy = std::abs((y - pFilmDiscrete.y) *
                    invFilterRadius.y * filterTableSize);
                ify[y - p0.y] = std::min((int)std::floor(fy), filterTableSize - 1);
            }
            for (int y = p0.y; y < p1.y; ++y) {
                for (int x = p0.x; x < p1.x; ++x) {
                    int offset = ify[y - p0.y] * filterTableSize + ifx[x - p0.x];
                    float filterWeight = filterTable[offset];
                    FilmTilePixel& pixel = GetPixel(Point2i(x, y));
                    pixel.contribSum += L * sampleWeight * filterWeight;
                    pixel.filterWeightSum += filterWeight;
                }
            }
        }
        Bounds2i GetPixelBounds() const { return pixelBounds; }
        FilmTilePixel& GetPixel(const Point2i& p) {
            int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
            int offset = (p.x - pixelBounds.pMin.x) +
                (p.y - pixelBounds.pMin.y) * width;
            return pixels[offset];
        }

    private:
        // FilmTile Private Data
        const Bounds2i pixelBounds;
        const GeoVector2f filterRadius, invFilterRadius;
        const float* filterTable;
        const int filterTableSize;
        std::vector<FilmTilePixel> pixels;
    };

    class Film {
    public:
        // Film Public Methods
        Film(const Point2i& resolution, const Bounds2f& cropWindow,
            std::unique_ptr<Filter> filter,
            float diagonal, const char* filename, float scale);
        Bounds2i GetSampleBounds() const;
        Bounds2f GetPhysicalExtent() const;
        std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i& sampleBounds);
        void MergeFilmTile(std::unique_ptr<FilmTile> tile);
        void SetImage(const Spectrum* img) const;
        void AddSplat(const Point2f& p, const Spectrum& v);
        void WriteImage(float splatScale = 1);
        void Clear();
        
        // Film Public Data
        const Point2i fullResolution;
        const float diagonal;
        std::unique_ptr<Filter> filter;
        const char* filename;
        Bounds2i croppedPixelBounds;

    private:
        // Film Private Data
        std::unique_ptr<Pixel[]> pixels;
        static constexpr int filterTableWidth = 16;
        float filterTable[filterTableWidth * filterTableWidth];
        std::mutex mutex;
        const float scale;

        // Film Private Methods
        Pixel& GetPixel(const Point2i& p) {
            int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
            int offset = (p.x - croppedPixelBounds.pMin.x) + (p.y - croppedPixelBounds.pMin.y) * width;
            return pixels[offset];
        }
    };

} // namespace lightfold