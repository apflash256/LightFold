#pragma once

#include <type_traits>

namespace lightfold {

    class RGB {
    public:
        float r;
        float g;
        float b;
    };

    class Spectra {
    public:
        virtual bool IsBlack() const = 0;
        virtual RGB toRGB() const = 0;
    };

    template <class C, std::enable_if_t<std::is_base_of<Spectra, C>::value>* = nullptr>
    inline const C operator*(float a, const C& s) {
        return s * a;
    }

    class BW : public Spectra {
    public:
        BW(float intensity = 0) : intensity(intensity) {}
        bool IsBlack() const { return intensity == 0; }
        RGB toRGB() const { return{ intensity, intensity, intensity }; }
        BW operator+(const BW& s) const { return { intensity + s.intensity }; }
        BW& operator+=(const BW& s) { intensity += s.intensity; }
        BW operator-(const BW& s) const { return { intensity - s.intensity }; }
        BW operator/(const BW& s) const { return { intensity / s.intensity }; }
        BW operator*(const BW& s) const { return { intensity * s.intensity }; }
        BW& operator*=(const BW& s) { intensity *= s.intensity; }
        BW operator*(float a) const { return { intensity * a }; }
        friend inline BW operator*(float a, const BW& s);
        BW& operator*=(float a) { intensity *= a; }
        BW operator/(float a) const { return { intensity / a }; }
        BW& operator/=(float a) { intensity /= a; }
        bool operator==(const BW& s) const { return { intensity == s.intensity }; }
        bool operator!=(const BW& s) const { return { intensity != s.intensity }; }
    private:
        float intensity;
    };

} // namespace lightfold