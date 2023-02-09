#include <core/integrator.h>

#include <core/bsdf.h>
#include <core/stats.h>
#include <display/progress.h>

namespace lightfold {

    Spectrum UniformSampleAllLights(const Interaction& it, const Scene& scene,
        Sampler& sampler, const std::vector<int>& nLightSamples, bool handleMedia) {
        ProfilePhase p(Prof::DirectLighting);
        Spectrum L(0.f);
        for (size_t j = 0; j < scene.lights.size(); ++j) {
            // Accumulate contribution of _j_th light to _L_
            const std::shared_ptr<Light>& light = scene.lights[j];
            int nSamples = nLightSamples[j];
            const Point2f* uLightArray = sampler.Get2DArray(nSamples);
            const Point2f* uScatteringArray = sampler.Get2DArray(nSamples);
            if (!uLightArray || !uScatteringArray) {
                // Use a single sample for illumination from _light_
                Point2f uLight = sampler.Get2D();
                Point2f uScattering = sampler.Get2D();
                L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
                    handleMedia);
            }
            else {
                // Estimate direct lighting using sample arrays
                Spectrum Ld(0.f);
                for (int k = 0; k < nSamples; ++k)
                    Ld += EstimateDirect(it, uScatteringArray[k], *light, uLightArray[k],
                        scene, sampler, handleMedia);
                L += Ld / nSamples;
            }
        }
        return L;
    }

    Spectrum UniformSampleOneLight(const Interaction& it, const Scene& scene,
        Sampler& sampler, bool handleMedia, const Distribution1D* lightDistrib) {
        ProfilePhase p(Prof::DirectLighting);
        // Randomly choose a single light to sample, _light_
        int nLights = int(scene.lights.size());
        if (nLights == 0) return Spectrum(0.f);
        int lightNum;
        float lightPdf;
        if (lightDistrib) {
            lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
            if (lightPdf == 0) return Spectrum(0.f);
        }
        else {
            lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
            lightPdf = float(1) / nLights;
        }
        const std::shared_ptr<Light>& light = scene.lights[lightNum];
        Point2f uLight = sampler.Get2D();
        Point2f uScattering = sampler.Get2D();
        return EstimateDirect(it, uScattering, *light, uLight, scene, sampler, handleMedia)
            / lightPdf;
    }

    Spectrum EstimateDirect(const Interaction& it, const Point2f& uScattering,
        const Light& light, const Point2f& uLight, const Scene& scene, Sampler& sampler,
        bool handleMedia) {
        Spectrum Ld(0.f);
        // Sample light source with multiple importance sampling
        Tangent3f wi;
        float lightPdf = 0, scatteringPdf = 0;
        VisibilityTester visibility;
        Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
        if (lightPdf > 0 && !Li.IsBlack()) {
            // Compute BSDF or phase function's value for light sample
            Spectrum f;
            if (it.IsSurfaceInteraction()) {
                // Evaluate BSDF for light sampling strategy
                const SurfaceInteraction& isect = (const SurfaceInteraction&)it;
                f = isect.bsdf->distF(isect.wo, wi) *
                    AbsDot(wi, isect.shading.n);
                scatteringPdf = isect.bsdf->Pdf(isect.wo, wi);
            }
            else {
                // Evaluate phase function for light sampling strategy
                const MediumInteraction& mi = (const MediumInteraction&)it;
                float p = mi.phase->p(mi.wo, wi);
                f = Spectrum(p);
                scatteringPdf = p;
            }
            if (!f.IsBlack()) {
                // Compute effect of visibility for light source sample
                if (handleMedia) {
                    Li *= visibility.Tr(scene, sampler);
                }
                else {
                    if (!visibility.Unoccluded(scene)) {
                        Li = Spectrum(0.f);
                    }
                }

                // Add light's contribution to reflected radiance
                if (!Li.IsBlack()) {
                    if (IsDeltaLight(light.flags))
                        Ld += f * Li / lightPdf;
                    else {
                        float weight =
                            PowerHeuristic(1, lightPdf, 1, scatteringPdf);
                        Ld += f * Li * weight / lightPdf;
                    }
                }
            }
        }

        // Sample BSDF with multiple importance sampling
        if (!IsDeltaLight(light.flags)) {
            Spectrum f;
            if (it.IsSurfaceInteraction()) {
                // Sample scattered direction for surface interactions
                const SurfaceInteraction& isect = (const SurfaceInteraction&)it;
                f = isect.bsdf->Sample_F(isect.wo, &wi, uScattering, &scatteringPdf);
                f *= AbsDot(wi, isect.shading.n);
            }
            else {
                // Sample scattered direction for medium interactions
                const MediumInteraction& mi = (const MediumInteraction&)it;
                float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
                f = Spectrum(p);
                scatteringPdf = p;
            }
            if (!f.IsBlack() && scatteringPdf > 0) {
                // Account for light contributions along sampled direction _wi_
                float weight = 1;
                /*if (specular) {
                    lightPdf = light.Pdf_Li(it, wi);
                    if (lightPdf == 0) return Ld;
                    weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
                }*/

                // Find intersection and compute transmittance
                SurfaceInteraction lightIsect;
                Ray ray = it.SpawnRay(wi);
                Spectrum Tr(1.f);
                bool foundSurfaceInteraction =
                    handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
                    : scene.Intersect(ray, &lightIsect);

                // Add light contribution from material sampling
                Spectrum Li(0.f);
                if (foundSurfaceInteraction) {
                    if (lightIsect.primitive->GetAreaLight() == &light)
                        Li = lightIsect.Le(-wi);
                }
                else
                    Li = light.Le(ray);
                if (!Li.IsBlack()) Ld += f * Li * Tr * weight / scatteringPdf;
            }
        }
        return Ld;
    }

    std::unique_ptr<Distribution1D> ComputeLightPowerDistribution(
        const Scene& scene) {
        if (scene.lights.empty()) return nullptr;
        std::vector<float> lightPower;
        for (const auto& light : scene.lights)
            lightPower.push_back(light->Power().Value());
        return std::unique_ptr<Distribution1D>(
            new Distribution1D(&lightPower[0], lightPower.size()));
    }

    void SamplerIntegrator::Render(const Scene& scene) {
        Preprocess(scene, *sampler);
        // Render image tiles in parallel

        // Compute number of tiles, _nTiles_, to use for parallel rendering
        Bounds2i sampleBounds = camera->film->GetSampleBounds();
        Vector2i sampleExtent = sampleBounds.Diagonal();
        const int tileSize = 16;
        Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
            (sampleExtent.y + tileSize - 1) / tileSize);
        ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
        ParallelInit();
        ParallelFor2D([&](Point2i tile) {
            // Render section of image corresponding to _tile_
            // Get sampler instance for tile
            int seed = tile.y * nTiles.x + tile.x;
            std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

            // Compute sample bounds for tile
            int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

            // Get _FilmTile_ for tile
            std::unique_ptr<FilmTile> filmTile = camera->film->GetFilmTile(tileBounds);

            // Loop over pixels in tile to render them
            for (Point2i pixel : tileBounds) {
                ProfilePhase pp(Prof::StartPixel);
                tileSampler->StartPixel(pixel);
                if (!InsideExclusive(pixel, pixelBounds))
                    continue;
                do {
                    CameraSample cameraSample = tileSampler->GetCameraSample(pixel);
                    RayDifferential ray;
                    float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(
                        1 / std::sqrt((float)tileSampler->samplesPerPixel));
                    Spectrum L(0.f);
                    if (rayWeight > 0)
                        L = Li(ray, scene, *tileSampler);
                    if (L.HasNaNs()) {
                        std::cout << "Not-a-number radiance value returned." << std::endl;
                        L = Spectrum(0.f);
                    }
                    else if (L.Value() < -1e-5) {
                        std::cout << "Negative luminance value, returned." << std::endl;
                        L = Spectrum(0.f);
                    }
                    else if (std::isinf(L.Value())) {
                        std::cout << "Infinite luminance value returned." << std::endl;
                        L = Spectrum(0.f);
                    }
                    filmTile->AddSample(cameraSample.pFilm, L, rayWeight);
                } while (tileSampler->StartNextSample());
            }
            camera->film->MergeFilmTile(std::move(filmTile));
            reporter.Update();
        }, nTiles);
        reporter.Done();
        std::cout << "Rendering finished" << std::endl;

        // Save final image after rendering
        camera->film->WriteImage();
        ParallelCleanup();
    }

} // namespace lightfold