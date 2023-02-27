#include <material/matte.h>

#include <core/bsdf.h>
#include <lobe/diffuse.h>

namespace lightfold {

    void MatteMaterial::ComputeScatteringFunctions(SurfaceInteraction* si,
        TransportMode mode, bool allowMultipleLobes) const {
        // Evaluate textures for _MatteMaterial_ material and allocate BRDF
        Spectrum col = color->Evaluate(*si);
        float rough = roughness->Evaluate(*si);
        si->bsdf = new BSDF(*si);
        if (!col.IsBlack()) {
            si->bsdf->Add(std::make_shared<Diffuse>(col, rough));
        }
    }

} // namespace lightfold