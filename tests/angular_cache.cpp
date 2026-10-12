// Independent phase formulas, cache identity, sparse octants, and invalidation.
#include "mlfma_workspace.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <omp.h>

double pi = std::acos(-1.0);
void theta_locations_weights(int,double*,double*);
void theta_locations_weights_lin(int,double*,double*);
void phi_locations_weights(int,double*,double*);
void interpolationweights(int,int,int,int,double*,double*,double*,int*);

static void require(bool ok, const char* message) {
    if (!ok) { std::fprintf(stderr,"FAIL %s\n",message); std::exit(1); }
}

int main() {
    MlfmaWorkspace workspace;
    double maximum = 0;
    for (bool linear : {false,true}) for (int nt : {1,3,7,13}) {
        if (linear && nt==1) continue;
        const int np=2*nt;
        auto& grid=workspace.grid(nt,np,linear);
        std::vector<double> theta(nt), tw(nt), phi(np), pw(np);
        if (linear) theta_locations_weights_lin(nt,tw.data(),theta.data());
        else theta_locations_weights(nt,tw.data(),theta.data());
        phi_locations_weights(np,pw.data(),phi.data());
        require(theta==grid.theta && tw==grid.theta_weights && phi==grid.phi && pw==grid.phi_weights,
                "sampling tables preserve locations and weights");
        for (double width : {0.125,0.123456789}) for (double k : {1.0,20.95845021951682}) {
            const auto& table=workspace.child_shifts(nt,np,width,k,linear);
            require(&table==&workspace.child_shifts(nt,np,width,k,linear),"phase cache identity");
            require(table.upward.size()==std::size_t(8*nt*np),"exact eight-octant payload");
            for (int octant=0;octant<8;++octant) {
                double shift[3];
                for (int axis=0;axis<3;++axis) shift[axis]=(octant&(1<<axis) ? 1 : -1)*width/4;
                const auto* phases=table.for_shift(shift);
                require(phases==table.upward.data()+octant*table.samples,"octant indexing");
                for (int t=0;t<nt;++t) for (int p=0;p<np;++p) {
                    const double dot=shift[0]*std::sin(theta[t])*std::cos(phi[p])
                                    +shift[1]*std::sin(theta[t])*std::sin(phi[p])+shift[2]*std::cos(theta[t]);
                    const auto expected=std::polar(1.0,-k*dot);
                    const auto downward=std::polar(1.0,k*dot);
                    const auto cached=phases[t*np+p];
                    maximum=std::max(maximum,std::abs(cached-expected));
                    require(std::abs(cached-expected)<=1e-14 &&
                            std::abs(std::conj(cached)-downward)<=1e-14,"both shift signs");
                    require(std::abs(cached*std::conj(cached)-1.0)<=1e-14,"unit phase magnitude");
                }
            }
        }
    }
    const auto& base=workspace.child_shifts(7,14,0.125,20.0);
    require(&base!=&workspace.child_shifts(7,14,0.25,20.0),"box width key");
    require(&base!=&workspace.child_shifts(7,14,0.125,21.0),"wavenumber key");
    require(&base!=&workspace.child_shifts(7,14,0.125,20.0,true),"sampling convention key");
    require(&base!=&workspace.child_shifts(8,16,0.125,20.0),"angular shape key");
    // Occupied-child order is unrelated to octant number, including missing octants.
    for (int octant : {7,0,5}) {
        double shift[3];
        for (int a=0;a<3;++a) {
            const double parent=17.25-a;
            const double child=parent+(octant&(1<<a) ? 1 : -1)*0.125/4;
            shift[a]=child-parent;
        }
        require(base.for_shift(shift)==base.upward.data()+octant*base.samples,"sparse translated children");
    }
    const auto& rounded=workspace.child_shifts(7,14,0.123456789,20.0);
    double shift[3];
    for (int a=0;a<3;++a) shift[a]=(1e9+0.123456789/4)-1e9;
    require(rounded.for_shift(shift)==nullptr,"cancellation uses direct fallback");
    shift[0]=0; shift[1]=shift[2]=0.125/4;
    require(base.for_shift(shift)==nullptr,"zero displacement cannot select an octant");
    shift[0]=std::numeric_limits<double>::quiet_NaN();
    require(base.for_shift(shift)==nullptr,"nonfinite displacement cannot select an octant");

    // Rebuild stencils independently, including distinct linear output sampling.
    for (bool linear : {false,true}) for (int order : {1,2,4}) {
        const int ct=5, cp=10, ft=9, fp=18;
        const auto& coarse=workspace.grid(ct,cp);
        const auto& fine=workspace.grid(ft,fp,linear);
        const auto& table=workspace.table(ct,cp,ft,fp,order,linear);
        require(&table==&workspace.table(ct,cp,ft,fp,order,linear),"interpolation cache identity");
        auto theta=coarse.theta, phi=coarse.phi, theta_fine=fine.theta, phi_fine=fine.phi;
        std::vector<double> wt(2*ft*order), wp(2*fp*order);
        std::vector<int> it(wt.size()), ip(wp.size());
        interpolationweights(ct,ft,order,0,theta.data(),theta_fine.data(),wt.data(),it.data());
        interpolationweights(cp,fp,order,1,phi.data(),phi_fine.data(),wp.data(),ip.data());
        require(wt==table.theta_weights && wp==table.phi_weights && it==table.theta_indices &&
                ip==table.phi_indices,"cached stencils equal direct construction");
    }
    const auto payload=workspace.retained_bytes();
    int bad=0;
    const double canonical[3]={0.125/4,-0.125/4,0.125/4};
    #pragma omp parallel for num_threads(32) reduction(+:bad)
    for (int j=0;j<10000;++j)
        bad += base.for_shift(canonical)!=base.upward.data()+5*base.samples;
    require(bad==0 && payload==workspace.retained_bytes(),"shared phase reads do not replicate storage");
    for (double invalid : {0.0,-1.0,std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        bool caught=false;
        try { workspace.child_shifts(7,14,invalid,20); } catch (const std::invalid_argument&) { caught=true; }
        require(caught,"invalid width rejected before map lookup");
        caught=false;
        try { workspace.child_shifts(7,14,0.125,invalid); } catch (const std::invalid_argument&) { caught=true; }
        require(caught,"invalid wavenumber rejected before map lookup");
    }
    bool caught=false;
    try { workspace.child_shifts(std::numeric_limits<int>::max(),2,0.125,20); }
    catch (const std::length_error&) { caught=true; }
    require(caught,"phase allocation bound checked before constructing grid");
    caught=false;
    try { workspace.child_shifts(3,6,std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max()); }
    catch (const std::overflow_error&) { caught=true; }
    require(caught,"nonfinite phase rejected before retaining table");
    const auto saved=base.upward;
    workspace.clear();
    require(workspace.retained_bytes()==0 && workspace.shift_count()==0,"phase cache cleanup");
    require(saved==workspace.child_shifts(7,14,0.125,20).upward,"phase rebuild identical");
    workspace.clear(); workspace.clear();
    std::printf("CHECK angular_cache max_phase_error %.17g cleanup_and_rebuild passed\n",maximum);
}
