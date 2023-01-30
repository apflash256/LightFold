#pragma once

#include <cmath>
#include <cstring>
#include <cstdint>

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

    inline float NextFloatUp(float v) {
        // Handle infinity and negative zero for _NextFloatUp()_
        if (std::isinf(v) && v > 0.) return v;
        if (v == -0.f) v = 0.f;

        // Advance _v_ to next higher float
        uint32_t ui = FloatToBits(v);
        if (v >= 0)
            ++ui;
        else
            --ui;
        return BitsToFloat(ui);
    }

    inline float NextFloatDown(float v) {
        // Handle infinity and positive zero for _NextFloatDown()_
        if (std::isinf(v) && v < 0.) return v;
        if (v == 0.f) v = -0.f;
        uint32_t ui = FloatToBits(v);
        if (v > 0)
            --ui;
        else
            ++ui;
        return BitsToFloat(ui);
    }

} // namespace lightfold