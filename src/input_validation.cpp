#include "input_validation.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

void require_input(bool condition, const std::string& field, const std::string& reason) {
    if (!condition) throw InputError(field, reason);
}
[[noreturn]] void report_input_error(const std::exception& error) {
    std::fprintf(stderr, "Input error: %s\n", error.what());
    std::exit(EXIT_FAILURE);
}
int maximum_mlfma_levels() {
    // dec2bin2dec multiplies a signed long by two three times per level,
    // including its final unused multiplier. Keep every intermediate representable.
    return std::min((std::numeric_limits<long>::digits-1)/3,
                    std::numeric_limits<int>::digits);
}
long checked_index_sum(long a, long b, long limit, const std::string& field) {
    require_input(a >= 0 && b >= 0 && a <= limit && b <= limit-a,
                  field,"array size exceeds integer index capacity");
    return a+b;
}
long checked_index_product(long a, long b, long limit, const std::string& field) {
    require_input(a >= 0 && b >= 0 && (b == 0 || a <= limit/b),
                  field,"array size exceeds integer index capacity");
    return a*b;
}
void validate_interpolation_dimensions(int coarse, int fine, int points) {
    require_input(coarse > 0 && fine > 0 && points > 0 && points <= coarse,
                  "num_interpol", "stencil requires 1 <= num_interpol <= coarse sample count (" + std::to_string(coarse) + ")");
    const long long limit = std::numeric_limits<int>::max();
    require_input(3LL*coarse <= limit && 2LL*fine*points <= limit,
                  "num_interpol", "extended locations or interpolation workspace exceed integer index capacity");
}
void validate_farfield_mode(double mode) {
    require_input(std::isfinite(mode), "dir_ff", "expected a finite selector/step");
    if (mode <= 0.0) return; // Negative: direct cuts; zero: monostatic.
    // Preserve the IBC routine's minimum step; validate its actual allocation/loops.
    const double step = std::max(0.1, mode);
    const double rows = 180.0/step+1.0, columns = 360.0/step+1.0;
    require_input(rows*columns <= std::numeric_limits<int>::max() &&
                  std::ceil(rows)*std::ceil(columns) <= std::floor(rows*columns),
                  "dir_ff", "direct-grid step is unsupported: choose a step dividing 180 and 360 degrees exactly");
}

MlfmaDimensions validate_mlfma_dimensions(int requested, double extent,
    double wavelength, double k, int digits, int pc, double approximation,
    int interpolation, double farfield) {
    require_input(std::isfinite(extent) && extent > 0.0, "POINTS", "geometry extent must be finite and positive");
    require_input(std::isfinite(wavelength) && wavelength > 0.0 && std::isfinite(k) && k > 0.0,
                  "Frequency", "derived wavelength and wave number must be finite and positive");
    require_input(requested >= 1 && requested <= maximum_mlfma_levels(),
                  "num_mlfmalevels", "expected 1 (automatic) or explicit levels within integer index capacity");
    MlfmaDimensions result;
    result.levels = requested;
    if (requested == 1) {
        const double ratio = extent/wavelength;
        require_input(std::isfinite(ratio) && ratio > 0.0, "num_mlfmalevels", "extent/wavelength is not representable");
        const double resolved = std::max(2.0, std::floor(std::log2(ratio))+3.0);
        require_input(std::isfinite(resolved) && resolved <= maximum_mlfma_levels(),
                      "num_mlfmalevels", "automatic selection exceeds integer index capacity; maximum is " + std::to_string(maximum_mlfma_levels()));
        result.levels = static_cast<int>(resolved);
    }
    result.boxes.resize(result.levels); result.terms.resize(result.levels);
    if (pc == 2) result.approximate_terms.resize(result.levels);
    result.boxes[0] = extent;
    for (int level = 0; level < result.levels; ++level) {
        if (level) result.boxes[level] = result.boxes[level-1]/2;
        require_input(std::isfinite(result.boxes[level]) && result.boxes[level] > 0.0,
                      "num_mlfmalevels", "box size underflow at level " + std::to_string(level));
        // Preserve the existing truncation formula and rounding.
        const double terms = std::ceil(1.73*std::abs(k)*result.boxes[level] +
            2.16*std::pow(digits,0.666666666666666)*std::pow(std::abs(k)*result.boxes[level],0.333333333333333));
        require_input(std::isfinite(terms) && terms >= 0.0 &&
                      terms < std::sqrt(std::numeric_limits<int>::max()/2.0)-1.0,
                      "truncation", "harmonic/sample count is not representable at level " + std::to_string(level));
        result.terms[level] = static_cast<int>(terms);
        if (pc == 2) {
            const double reduced = approximation*result.terms[level];
            require_input(std::isfinite(reduced) && reduced >= 0.0,
                          "truncation", "val_approx produces invalid harmonics at level " + std::to_string(level));
            result.approximate_terms[level] = static_cast<int>(reduced);
        }
    }
    auto validate_grid = [&](const std::vector<int>& terms, bool accelerated_field) {
        for (int level = 0; level < result.levels; ++level) {
            const int theta = terms[level]+1, phi = theta*2;
            const bool used = accelerated_field || (level >= 2 && result.levels >= 4);
            if (!used) continue;
            validate_interpolation_dimensions(theta, theta, interpolation);
            validate_interpolation_dimensions(phi, phi, interpolation);
            // Neighboring levels can have different fine grids in both directions.
            if (level+1 < result.levels) {
                const int next_theta = terms[level+1]+1, next_phi = next_theta*2;
                validate_interpolation_dimensions(theta, next_theta, interpolation);
                validate_interpolation_dimensions(phi, next_phi, interpolation);
                validate_interpolation_dimensions(next_theta, theta, interpolation);
                validate_interpolation_dimensions(next_phi, phi, interpolation);
            }
        }
        if (accelerated_field) {
            validate_interpolation_dimensions(terms[0]+1,181,interpolation);
            validate_interpolation_dimensions(2*(terms[0]+1),360,interpolation);
        }
    };
    validate_grid(result.terms, false); // This variant has no accelerated far-field output.
    if (pc == 2 && approximation > 0.0) validate_grid(result.approximate_terms, false);
    return result;
}
