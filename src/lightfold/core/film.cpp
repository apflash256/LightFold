#include <core/film.h>
#include <core/imageio.h>

namespace lightfold {

    Film::Film(const Point2i& resolution, const Bounds2f& cropWindow,
        std::unique_ptr<Filter> filt, float diagonal,
        const char* filename, float scale)
        : fullResolution(resolution), diagonal(diagonal * .001f),
        filter(std::move(filt)), filename(filename), scale(scale) {
        croppedPixelBounds =
            Bounds2i(Point2i(std::ceil(fullResolution.x * cropWindow.pMin.x),
                             std::ceil(fullResolution.y * cropWindow.pMin.y)),
                     Point2i(std::ceil(fullResolution.x * cropWindow.pMax.x),
                             std::ceil(fullResolution.y * cropWindow.pMax.y)));
        pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
        int offset = 0;
        for (int y = 0; y < filterTableWidth; ++y) {
            for (int x = 0; x < filterTableWidth; ++x, ++offset) {
                Point2f p;
                p.x = (x + 0.5f) * filter->radius.x / filterTableWidth;
                p.y = (y + 0.5f) * filter->radius.y / filterTableWidth;
                filterTable[offset] = filter->Evaluate(p);
            }
        }
    }

    Bounds2i Film::GetSampleBounds() const {
        Bounds2f floatBounds(
            Floor(Point2f(croppedPixelBounds.pMin) + GeoVector2f(0.5f, 0.5f) -
                filter->radius),
            Ceil(Point2f(croppedPixelBounds.pMax) - GeoVector2f(0.5f, 0.5f) +
                filter->radius));
        return (Bounds2i)floatBounds;
    }

    Bounds2f Film::GetPhysicalExtent() const {
        float aspect = (float)fullResolution.y / (float)fullResolution.x;
        float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
        float y = aspect * x;
        return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
    }

    std::unique_ptr<FilmTile> Film::GetFilmTile(
        const Bounds2i& sampleBounds) {
        GeoVector2f halfPixel = GeoVector2f(0.5f, 0.5f);
        Bounds2f floatBounds = (Bounds2f)sampleBounds;
        Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
        Point2i p1 = (Point2i)(Floor(floatBounds.pMax - halfPixel + filter->radius) + Point2i(1, 1));
        Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
        return std::unique_ptr<FilmTile>(new FilmTile(tilePixelBounds, filter->radius, filterTable, filterTableWidth));
    }

    void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
        std::lock_guard<std::mutex> lock(mutex);
        for (Point2i pixel : tile->GetPixelBounds()) {
            const FilmTilePixel& tilePixel = tile->GetPixel(pixel);
            Pixel& mergePixel = GetPixel(pixel);
            RGB rgb = tilePixel.contribSum.toRGB();
            mergePixel.rgb[0] += rgb.r;
            mergePixel.rgb[1] += rgb.g;
            mergePixel.rgb[2] += rgb.b;
            mergePixel.filterWeightSum += tilePixel.filterWeightSum;
        }
    }

    void Film::SetImage(const Spectrum* img) const {
        int nPixels = croppedPixelBounds.Area();
        RGB rgb;
        for (int i = 0; i < nPixels; ++i) {
            Pixel& p = pixels[i];
            rgb = img[i].toRGB();
            p.rgb[0] = rgb.r;
            p.rgb[1] = rgb.g;
            p.rgb[2] = rgb.b;
            p.filterWeightSum = 1;
            p.splatRGB[0] = p.splatRGB[1] = p.splatRGB[2] = 0;
        }
    }

    void Film::AddSplat(const Point2f& p, const Spectrum& v) {
        if (!InsideExclusive((Point2i)p, croppedPixelBounds))
            return;
        RGB rgb = v.toRGB();
        float rgbArray[3] = { rgb.r,rgb.g,rgb.b };
        Pixel& pixel = GetPixel((Point2i)p);
        for (int i = 0; i < 3; ++i)
            pixel.splatRGB[i].Add(rgbArray[i]);
    }

    void Film::WriteImage(float splatScale) {
        std::unique_ptr<float[]> rgb(new float[3 * croppedPixelBounds.Area()]);
        int offset = 0;
        for (Point2i p : croppedPixelBounds) {
            Pixel& pixel = GetPixel(p);
            rgb[3 * offset] = pixel.rgb[0];
            rgb[3 * offset + 1] = pixel.rgb[1];
            rgb[3 * offset + 2] = pixel.rgb[2];
            float filterWeightSum = pixel.filterWeightSum;
            if (filterWeightSum != 0) {
                float invWt = (float)1 / filterWeightSum;
                rgb[3 * offset] = std::max((float)0, rgb[3 * offset] * invWt);
                rgb[3 * offset + 1] = std::max((float)0, rgb[3 * offset + 1] * invWt);
                rgb[3 * offset + 2] = std::max((float)0, rgb[3 * offset + 2] * invWt);
            }
            float splatRGB[3] = { pixel.splatRGB[0], pixel.splatRGB[1], pixel.splatRGB[2] };
            rgb[3 * offset] += splatScale * splatRGB[0];
            rgb[3 * offset + 1] += splatScale * splatRGB[1];
            rgb[3 * offset + 2] += splatScale * splatRGB[2];
            rgb[3 * offset] *= scale;
            rgb[3 * offset + 1] *= scale;
            rgb[3 * offset + 2] *= scale;
            ++offset;
        }
        WriteEXR(std::move(rgb), croppedPixelBounds.Diagonal()[0], croppedPixelBounds.Diagonal()[1], filename);
    }

} // namespace lightfold