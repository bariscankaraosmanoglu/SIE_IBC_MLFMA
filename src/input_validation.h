#ifndef INPUT_VALIDATION_H
#define INPUT_VALIDATION_H
#include <stdexcept>
#include <string>
#include <vector>

class InputError : public std::runtime_error {
public:
    InputError(const std::string& field, const std::string& reason)
        : std::runtime_error(field + ": " + reason) {}
};
void require_input(bool condition, const std::string& field, const std::string& reason);
[[noreturn]] void report_input_error(const std::exception& error);
int maximum_mlfma_levels();
void validate_interpolation_dimensions(int coarse, int fine, int points);
void validate_farfield_mode(double mode);
long checked_index_sum(long a, long b, long limit, const std::string& field);
long checked_index_product(long a, long b, long limit, const std::string& field);

struct MlfmaDimensions {
    int levels;
    std::vector<double> boxes;
    std::vector<int> terms, approximate_terms;
};
MlfmaDimensions validate_mlfma_dimensions(int requested_levels, double extent,
    double wavelength, double wave_number, int digits, int pc, double approximation,
    int interpolation, double farfield);
#endif
