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

    class RGBSpectrum : public Spectra {
    public:
        // RGBSpectrum Public Methods
        RGBSpectrum(float intensity = 0) : r(intensity), g(intensity), b(intensity) {}
        RGBSpectrum(float r, float g, float b) : r(r), g(g), b(b) {}

        bool IsBlack() const { return (r == 0) && (g == 0) && (b == 0); }
        RGB ToRGB() const { return{ r, g, b }; }
        void FromRGB(const RGB& rgb) {
            r = rgb.r;
            g = rgb.g;
            b = rgb.b;
        }
        bool HasNaNs() const {
            return std::isnan(r) || std::isnan(g) || std::isnan(b);
        }
        float Value() const {
            return (r + g + b) / 3;
        }

        RGBSpectrum operator+(const RGBSpectrum& s) const { return { r + s.r, g + s.g, b + s.b }; }
        RGBSpectrum& operator+=(const RGBSpectrum& s) {
            r += s.r;
            g += s.g;
            b += s.b;
            return *this;
        }
        RGBSpectrum operator-(const RGBSpectrum& s) const { return { r - s.r, g - s.g, b - s.b }; }
        RGBSpectrum operator/(const RGBSpectrum& s) const { return { r / s.r, g / s.g, b / s.b }; }
        RGBSpectrum operator*(const RGBSpectrum& s) const { return { r * s.r, g * s.g, b * s.b }; }
        RGBSpectrum& operator*=(const RGBSpectrum& s) {
            r *= s.r;
            g *= s.g;
            b *= s.b;
            return *this;
        }
        RGBSpectrum operator*(float a) const { return { r * a , g * a, b * a }; }
        friend inline RGBSpectrum operator*(float a, const RGBSpectrum& s) {
            return s * a;
        }
        RGBSpectrum& operator*=(float a) {
            r *= a;
            g *= a;
            b *= a;
            return *this;
        }
        RGBSpectrum operator/(float a) const { return { r / a , g / a, b / a }; }
        RGBSpectrum& operator/=(float a) {
            r /= a;
            g /= a;
            b /= a;
            return *this;
        }
        bool operator==(const RGBSpectrum& s) const { return (r == s.r) && (g == s.g) && (b == s.b); }
        bool operator!=(const RGBSpectrum& s) const { return (r != s.r) || (g != s.g) || (b != s.b); }

        float r, g, b;
    };

    using Spectrum = RGBSpectrum;

} // namespace lightfold