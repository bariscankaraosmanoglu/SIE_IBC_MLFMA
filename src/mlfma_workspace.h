#ifndef MLFMA_WORKSPACE_H
#define MLFMA_WORKSPACE_H

#include <complex>
#include <map>
#include <tuple>
#include <vector>

// Only sampling/interpolation metadata is retained, never per-RWG patterns.
// The solver's global aggregation arrays already require sequential operator
// calls. Workers share immutable metadata and own independent scratch buffers.
struct AngularGrid {
    std::vector<double> theta, theta_weights, phi, phi_weights;
    std::vector<double> sin_theta, cos_theta, sin_phi, cos_phi;
};
struct InterpolationTable {
    std::vector<double> theta_weights, phi_weights;
    std::vector<int> theta_indices, phi_indices;
};
struct ChildShiftTable {
    double parent_width, wavenumber;
    std::size_t samples;
    // Eight octants, theta-major/phi-minor, shared by both polarizations.
    // Upward shifts use these values; downward shifts use their conjugates.
    std::vector<std::complex<double>> upward;
    // Returns nullptr if center subtraction is not sufficiently close to the
    // canonical octree offset. Callers then preserve the direct shift formula.
    const std::complex<double>* for_shift(const double displacement[3]) const;
};
struct InterpolationScratch {
    std::vector<std::complex<double>> inflated_theta, half, inflated_phi;
    void prepare(int coarse_theta, int coarse_phi, int fine_theta);
};
struct MlfmaScratch {
    std::vector<std::complex<double>> fine_theta, fine_phi, coarse_theta, coarse_phi;
    InterpolationScratch interpolation;
    void prepare(int coarse_thetas, int coarse_phis, int fine_thetas, int fine_phis);
};
class MlfmaWorkspace {
    using GridKey = std::tuple<int,int,bool>;
    using TableKey = std::tuple<int,int,int,int,int,bool>;
    using ShiftKey = std::tuple<int,int,bool,double,double>;
    // Node-based maps keep pointers/references stable while another grid is added.
    std::map<GridKey, AngularGrid> grids_;
    std::map<TableKey, InterpolationTable> tables_;
    std::map<ShiftKey, ChildShiftTable> shifts_;
    std::vector<MlfmaScratch> workers_;
public:
    AngularGrid& grid(int theta, int phi, bool linear_theta = false);
    InterpolationTable& table(int coarse_theta, int coarse_phi, int fine_theta,
                              int fine_phi, int order, bool linear_fine = false);
    const ChildShiftTable& child_shifts(int theta, int phi, double parent_width,
                                      double wavenumber, bool linear_theta = false);
    // Every worker calls this before obtaining scratch; single's implicit
    // barrier publishes pool growth using the actual team, including serialized teams.
    void prepare_team();
    MlfmaScratch& worker();
    void clear();
    std::size_t retained_bytes() const;
    std::size_t grid_count() const { return grids_.size(); }
    std::size_t table_count() const { return tables_.size(); }
    std::size_t shift_count() const { return shifts_.size(); }
    std::size_t worker_count() const { return workers_.size(); }
};
MlfmaWorkspace& mlfma_workspace();

// Existing interpolation entry points remain available; production MLFMA paths
// use these overloads with caller-owned, reusable scratch.
void interpolation_c2f(int,int,int,int,int,int*,int*,double*,double*,
                       std::complex<double>*,std::complex<double>*,InterpolationScratch&);
void interpolation_f2c(int,int,int,int,int,int*,int*,double*,double*,
                       std::complex<double>*,std::complex<double>*,InterpolationScratch&);
#endif
