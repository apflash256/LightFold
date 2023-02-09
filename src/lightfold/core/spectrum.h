#pragma once

#include <cmath>

namespace lightfold {

    class RGB {
    public:
        float r;
        float g;
        float b;
    };

    class Spectra {
    public:
        // Spectra Interface
        virtual RGB ToRGB() const = 0;
        virtual void FromRGB(const RGB& rgb) = 0;
        virtual bool IsBlack() const = 0;
        virtual bool HasNaNs() const = 0;
        virtual float Value() const = 0;
    };

    class BW : public Spectra {
    public:
        // BW Public Methods
        BW(float intensity = 0) : intensity(intensity) {}

        bool IsBlack() const { return intensity == 0; }
        RGB ToRGB() const { return{ intensity, intensity, intensity }; }
        void FromRGB(const RGB& rgb) {
            this->intensity = (rgb.r + rgb.g + rgb.b) / 3;
        }
        bool HasNaNs() const {
            return std::isnan(intensity);
        }
        float Value() const {
            return intensity;
        }

        BW operator+(const BW& s) const { return { intensity + s.intensity }; }
        BW& operator+=(const BW& s) {
            intensity += s.intensity;
            return *this;
        }
        BW operator-(const BW& s) const { return { intensity - s.intensity }; }
        BW operator/(const BW& s) const { return { intensity / s.intensity }; }
        BW operator*(const BW& s) const { return { intensity * s.intensity }; }
        BW& operator*=(const BW& s) {
            intensity *= s.intensity;
            return *this;
        }
        BW operator*(float a) const { return { intensity * a }; }
        friend inline BW operator*(float a, const BW& s) {
            return s * a;
        }
        BW& operator*=(float a) {
            intensity *= a;
            return *this;
        }
        BW operator/(float a) const { return { intensity / a }; }
        BW& operator/=(float a) {
            intensity /= a;
            return *this;
        }
        bool operator==(const BW& s) const { return { intensity == s.intensity }; }
        bool operator!=(const BW& s) const { return { intensity != s.intensity }; }

    private:
        // BW Private Data
        float intensity;
    };

    using Spectrum = BW;

} // namespace lightfold