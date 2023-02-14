#include <material/matte.h>
#include <bsdf/diffuse.h>

namespace lightfold {

    void MatteMaterial::ComputeScatteringFunctions(SurfaceInteraction* si,
        TransportMode mode, bool allowMultipleLobes) const {
        // Evaluate textures for _MatteMaterial_ material and allocate BRDF
        Spectrum col = color->Evaluate(*si);
        float rough = roughness->Evaluate(*si);
        if (!col.IsBlack()) {
            si->bsdf = new Diffuse(*si, col, rough);
        }
    }

} // namespace lightfold