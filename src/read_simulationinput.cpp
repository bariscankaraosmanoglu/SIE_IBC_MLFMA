#include "solver_globals.h"
#include "input_validation.h"
#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>

using namespace std;

extern char modelname[];
extern int num_mlfmalevels;
extern int num_digits;
extern int num_interpol;
extern int source_type;
extern int val_itermax;
extern double val_soltol;
extern int num_media;
extern complex<double> val_impedance[];
extern double val_freq;
extern int num_freqs;
extern double freq_start, freq_end;
extern complex<double> *source_CSP;
extern complex<double> *source_current;
extern complex<double> *source_DG;
extern double *source_dipole;
extern int *num_source;
extern int val_pol;
extern int val_pc;
extern double val_approx;
extern int val_itermax_pc;
extern int val_typeff;
extern int *num_ffdir_tp;
extern double dir_ff;
extern double dir_cfc[];
extern int parallel_omp_on;
extern int num_threads;

extern int num_nodes;
extern int num_tris;
extern double *node_pos;
extern int *tri_nodes;
extern int *tri_media;

extern double pi;
extern double val_e;
extern double val_epsilonO;
extern double val_muO;
extern double val_omega;
extern double val_c;
extern double val_lambdaO;
extern double val_kO;
extern double val_etaO;


namespace {
class Reader {
    ifstream file{"input.pnu"};
    size_t line = 0;
    std::streamoff bytes = 0;
public:
    vector<string> record;
    Reader() {
        require_input(bool(file), "input.pnu", "cannot open input file");
        file.seekg(0,ios::end); bytes=file.tellg(); file.seekg(0);
        require_input(bool(file) && bytes >= 0,"input.pnu","cannot determine input file size");
    }
    string where(const string& field) const {
        return "input.pnu:" + to_string(line) + ": " + field;
    }
    bool next() {
        string text;
        while (getline(file,text)) {
            ++line; record.clear();
            istringstream words(text); string word;
            while (words >> word) record.push_back(word);
            if (!record.empty()) return true;
        }
        require_input(!file.bad(),where("input"),"read failure");
        return false;
    }
    void fields(const string& field, size_t count) {
        if (record.size() != count) throw InputError(where(field),"expected " + to_string(count) + " fields; found " + to_string(record.size()));
    }
    void expected(const string& field, size_t count) {
        if (!next()) throw InputError(where(field),"unexpected end of file");
        if (record[0] != field) throw InputError(where(field),"expected parameter '" + field + "'; found '" + record[0] + "'");
        fields(field,count);
    }
    void seek(const string& field, size_t count) {
        while (next()) {
            if (record[0] == field) { fields(field,count); return; }
            require_input(record[0] != "POINTS" && record[0] != "CELLS",where(field),"unexpected section '" + record[0] + "'");
        }
        throw InputError(where(field),"required section not found");
    }
    void records_available(int count, int minimum_bytes, const string& field) {
        const std::streamoff position = file.tellg();
        check(position >= 0 && count <= (bytes-position)/minimum_bytes,field,
              "declared record count exceeds the remaining input data");
    }
    int integer(size_t index, const string& field) const {
        const string& token = record[index];
        char* end = nullptr; errno = 0;
        const long long value = strtoll(token.c_str(),&end,10);
        if (errno == ERANGE || end == token.c_str() || *end != '\0' || value < INT_MIN || value > INT_MAX)
            throw InputError(where(field),"invalid integer '" + token + "'");
        return static_cast<int>(value);
    }
    double real(size_t index, const string& field) const {
        const string& token = record[index];
        char* end = nullptr; errno = 0;
        const double value = strtod(token.c_str(),&end);
        if (end == token.c_str() || *end != '\0' || !isfinite(value))
            throw InputError(where(field),"expected finite number; found '" + token + "'");
        return value;
    }
    void check(bool condition, const string& field, const char* reason) const {
        if (!condition) throw InputError(where(field),reason);
    }
    void check(bool condition, const string& field, const string& reason) const {
        if (!condition) throw InputError(where(field),reason);
    }
};
struct Input {
    string name;
    int formulation, source, frequencies, sources[2], polarization, pc, levels, digits;
    int interpolation, inner_iterations, iterations, parallel, threads, materials, nodes, triangles;
    double start, end, pw[4], field_mode, approximation, tolerance, contribution[2];
    complex<double> impedance[50]{};
    unique_ptr<double[]> points;
    unique_ptr<int[]> cells, media;
};
Input parse_input() {
    Reader r; Input in;
    r.seek("SimulationInputs",1);
    r.expected("model_name",2); in.name = r.record[1];
    r.check(in.name.size() < 200,"model_name","maximum model name length is 199 bytes");
    r.expected("Formulations",2); in.formulation = r.integer(1,"Formulations");
    r.check(in.formulation == 1,"Formulations",
            "this IBC executable implements mode 1; mode 0 CFIE assembly is unavailable");
    r.expected("Frequency",4);
    in.start=r.real(1,"Frequency"); in.end=r.real(2,"Frequency"); in.frequencies=r.integer(3,"Frequency");
    r.check(in.start > 0.0 && in.end > 0.0,"Frequency","frequency endpoints must be positive");
    r.check(in.frequencies == 1,"Frequency","frequency sweeps are unsupported; this executable requires count 1 and uses the start frequency");
    r.expected("source_type",2); in.source=r.integer(1,"source_type");
    const char* source_names[] = {"plane wave","delta gap","dipole","current source","complex source point"};
    r.check(in.source == 0,"source_type","unsupported illumination " + to_string(in.source) +
            (in.source >= 0 && in.source < 5 ? " (" + string(source_names[in.source]) + ")" : "") + "; this executable supports mode 0 (plane wave)");
    r.expected("source_PW",8);
    in.pw[0]=r.real(1,"source_PW"); in.pw[1]=r.real(2,"source_PW"); in.sources[0]=r.integer(3,"source_PW");
    in.pw[2]=r.real(4,"source_PW"); in.pw[3]=r.real(5,"source_PW"); in.sources[1]=r.integer(6,"source_PW");
    in.polarization=r.integer(7,"source_PW");
    r.check(in.sources[0] > 0 && in.sources[1] > 0 && 1LL*in.sources[0]*in.sources[1] <= INT_MAX,
            "source_PW","illumination counts must be positive and their product must fit an integer");
    r.check(in.polarization == 0 || in.polarization == 1,"source_PW","polarization must be 0 or 1");
    r.check(isfinite(in.pw[1]-in.pw[0]) && isfinite(in.pw[3]-in.pw[2]),"source_PW","angle range is not representable");
    r.expected("Preconditioners",2); in.pc=r.integer(1,"Preconditioners");
    r.check(in.pc >= 0 && in.pc <= 2,"Preconditioners","supported modes are 0, 1, and 2");
    r.expected("num_mlfmalevels",2); in.levels=r.integer(1,"num_mlfmalevels");
    r.check(in.levels >= 1 && in.levels <= maximum_mlfma_levels(),"num_mlfmalevels",
            "expected 1 (automatic) or explicit levels up to " + to_string(maximum_mlfma_levels()) + " for the current integer indices");
    r.expected("dir_ff",2); in.field_mode=r.real(1,"dir_ff"); validate_farfield_mode(in.field_mode);
    r.expected("num_digits",2); in.digits=r.integer(1,"num_digits");
    r.check(in.digits > 0,"num_digits","expected a positive accuracy parameter");
    r.expected("val_approx",2); in.approximation=r.real(1,"val_approx");
    r.check(in.pc != 2 || (in.approximation >= 0.0 && in.approximation <= 1.0),"val_approx","preconditioner 2 requires 0 <= val_approx <= 1");
    r.expected("num_interpol",2); in.interpolation=r.integer(1,"num_interpol");
    r.check(in.interpolation > 0 && in.interpolation <= INT_MAX/2,"num_interpol","expected positive stencil width within integer index capacity");
    r.expected("val_itermax_pc",2); in.inner_iterations=r.integer(1,"val_itermax_pc");
    r.check(in.pc != 2 || in.inner_iterations > 0,"val_itermax_pc","preconditioner 2 requires a positive inner iteration limit");
    r.expected("val_soltol",2); in.tolerance=r.real(1,"val_soltol");
    r.check(in.tolerance > 0.0 && in.tolerance < 1.0,"val_soltol","expected 0 < relative tolerance < 1");
    r.expected("dir_cfc",3); in.contribution[0]=r.real(1,"dir_cfc"); in.contribution[1]=r.real(2,"dir_cfc");
    r.expected("val_itermax",2); in.iterations=r.integer(1,"val_itermax");
    r.check(in.iterations > 0,"val_itermax","expected a positive outer iteration limit");
    r.expected("parallel_omp",3); in.parallel=r.integer(1,"parallel_omp"); in.threads=r.integer(2,"parallel_omp");
    r.check(in.parallel == 0 || in.parallel == 1,"parallel_omp","enable flag must be 0 or 1");
    r.check(in.threads > 0 && in.threads <= INT_MAX/2,"parallel_omp","thread count must be positive and its range indices representable");
    r.expected("##",2); r.check(r.record[1] == "End-Param","End-Param","expected '## End-Param'");
    r.expected("ModelParts",2); in.materials=r.integer(1,"ModelParts");
    r.check(in.materials > 0 && in.materials <= 50,"ModelParts","expected between 1 and 50 material parts");
    for (int part = 0; part < in.materials; ++part) {
        r.check(r.next(),"ModelParts","unexpected end of material list"); r.fields("ModelParts",3);
        const string field = "ModelParts " + r.record[0];
        r.check(r.record[0].size() < 200,field,"maximum part name length is 199 bytes");
        in.impedance[part]={r.real(1,field),r.real(2,field)};
        const complex<double> value=in.impedance[part];
        r.check(value != complex<double>{},field,"IBC coefficient must be nonzero (reciprocal required)");
        const complex<double> inverse=1.0/value;
        const complex<double> direct_root=sqrt(value), inverse_root=sqrt(inverse);
        r.check(isfinite(inverse.real()) && isfinite(inverse.imag()) &&
                isfinite(direct_root.real()) && isfinite(direct_root.imag()) &&
                isfinite(inverse_root.real()) && isfinite(inverse_root.imag()),field,
                "IBC reciprocal and square roots must be representable");
    }
    r.seek("POINTS",2); in.nodes=r.integer(1,"POINTS");
    r.check(in.nodes >= 3 && in.nodes <= INT_MAX/3,"POINTS","point count must be at least 3 and fit coordinate index capacity");
    r.records_available(in.nodes,5,"POINTS");
    in.points.reset(new double[size_t(in.nodes)*3]);
    for (int point = 0; point < in.nodes; ++point) {
        r.check(r.next(),"POINTS","unexpected end of point list at point " + to_string(point)); r.fields("POINTS",3);
        const string field="POINTS point " + to_string(point);
        for (int axis = 0; axis < 3; ++axis) in.points[size_t(point)*3+axis]=r.real(axis,field);
    }
    r.seek("CELLS",2); in.triangles=r.integer(1,"CELLS");
    r.check(in.triangles > 0 && in.triangles <= INT_MAX/6,"CELLS","triangle count must be positive and fit triangle/RWG index capacity");
    r.records_available(in.triangles,7,"CELLS");
    in.cells.reset(new int[size_t(in.triangles)*3]); in.media.reset(new int[in.triangles]);
    set<array<int,3>> unique_triangles;
    for (int triangle = 0; triangle < in.triangles; ++triangle) {
        const string field = "CELLS triangle " + to_string(triangle);
        r.check(r.next(),field,"unexpected end of triangle list"); r.fields(field,4);
        in.media[triangle]=r.integer(0,field);
        r.check(in.media[triangle] >= 0 && in.media[triangle] < in.materials,field,"material index is outside ModelParts");
        array<int,3> nodes;
        for (int corner = 0; corner < 3; ++corner) {
            nodes[corner]=r.integer(corner+1,field);
            r.check(nodes[corner] >= 0 && nodes[corner] < in.nodes,field,"node index is outside POINTS");
            in.cells[size_t(triangle)*3+corner]=nodes[corner];
        }
        r.check(nodes[0] != nodes[1] && nodes[1] != nodes[2] && nodes[0] != nodes[2],field,"repeated node");
        double edges[3][3];
        for (int edge = 0; edge < 3; ++edge) {
            const int first=edge<2 ? nodes[0] : nodes[1];
            const int second=edge<2 ? nodes[edge+1] : nodes[2];
            for (int axis = 0; axis < 3; ++axis) edges[edge][axis]=in.points[size_t(second)*3+axis]-in.points[size_t(first)*3+axis];
            const double length=sqrt(edges[edge][0]*edges[edge][0]+edges[edge][1]*edges[edge][1]+edges[edge][2]*edges[edge][2]);
            r.check(isfinite(length) && length > 0.0,field,"edge length is zero or not representable");
        }
        const double cross[3]={edges[0][1]*edges[1][2]-edges[0][2]*edges[1][1],
                               edges[0][2]*edges[1][0]-edges[0][0]*edges[1][2],
                               edges[0][0]*edges[1][1]-edges[0][1]*edges[1][0]};
        const double area2=sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
        r.check(isfinite(area2) && area2 > 0.0,field,"triangle area is zero or not representable");
        sort(nodes.begin(),nodes.end());
        r.check(unique_triangles.insert(nodes).second,field,"duplicate triangle");
    }
    // The shared file format may carry later sections. Reject surplus numeric
    // cell records, while retaining compatibility with named trailing metadata.
    if (r.next()) {
        const char first=r.record[0][0];
        r.check(!(first >= '0' && first <= '9') && first != '-' && first != '+' &&
                r.record[0] != "CELLS" && r.record[0] != "POINTS","CELLS","unexpected extra mesh record after declared triangles");
    }
    return in;
}
}

void read_simulationinput() {
    try {
        Input in=parse_input();
        // Validate constants before exposing any parsed arrays to geometry/assembly.
        const double omega=2*pi*in.start, c=1/sqrt(val_epsilonO*val_muO), wavelength=c/in.start;
        const double k=2*pi/wavelength, eta=sqrt(val_muO/val_epsilonO);
        require_input(isfinite(omega) && omega > 0.0 && isfinite(wavelength) && wavelength > 0.0 &&
                      isfinite(k) && k > 0.0,"Frequency","derived electromagnetic constants are not representable");
        strcpy(modelname,in.name.c_str());
        num_mlfmalevels=in.levels; num_digits=in.digits; num_interpol=in.interpolation;
        source_type=in.source; num_freqs=in.frequencies; freq_start=in.start; freq_end=in.end;
        copy(in.pw,in.pw+4,source_PW); copy(in.sources,in.sources+2,num_source); val_pol=in.polarization;
        val_pc=in.pc; dir_ff=in.field_mode; val_approx=in.approximation; val_itermax_pc=in.inner_iterations;
        val_soltol=in.tolerance; copy(in.contribution,in.contribution+2,dir_cfc); val_itermax=in.iterations;
        parallel_omp_on=in.parallel; num_threads=in.threads; num_media=in.materials;
        copy(in.impedance,in.impedance+in.materials,val_impedance);
        num_nodes=in.nodes; num_tris=in.triangles;
        node_pos=in.points.release(); tri_nodes=in.cells.release(); tri_media=in.media.release();
        val_freq=in.start; val_omega=omega; val_c=c; val_lambdaO=wavelength; val_kO=k; val_etaO=eta;
        printf("model_name %s\nFormulations %d\nFrequency %.17g %.17g %d\nsource_type %d\n",
               modelname,in.formulation,freq_start,freq_end,num_freqs,source_type);
        printf("source_PW %.17g %.17g %d %.17g %.17g %d %d\n",source_PW[0],source_PW[1],num_source[0],source_PW[2],source_PW[3],num_source[1],val_pol);
        printf("Preconditioners %d\nnum_mlfmalevels %d\ndir_ff %.17g\nnum_digits %d\nval_approx %.17g\nnum_interpol %d\nval_itermax_pc %d\nval_soltol %.17g\n",
               val_pc,num_mlfmalevels,dir_ff,num_digits,val_approx,num_interpol,val_itermax_pc,val_soltol);
        printf("dir_cfc %.17g %.17g\nval_itermax %d\nparallel_omp %d %d\nModelParts %d\n",
               dir_cfc[0],dir_cfc[1],val_itermax,parallel_omp_on,num_threads,num_media);
        for(int part=0;part<num_media;++part) printf("val_impedance[%d]: %.17g %.17g\n",part,val_impedance[part].real(),val_impedance[part].imag());
        printf("Number of points: %d\nNumber of cells: %d\nFrequency: %f GHz\nIllumination Type: Plane Wave\n",
               num_nodes,num_tris,val_freq/1e9);
        if(freq_end != freq_start) printf("Single-frequency count 1 uses Frequency start; end is retained in the shared input format.\n");
        printf("PI: %1.20f\nOmega: %f rad Hz\nEpsilon: %1.4e F/m\nMu: %1.4e H/m\nSpeed of Light: %1.4e m/s\nWavelength %f m\nWave number: %f rad/m\nWave Impedance: %f ohm\n\n",
               pi,val_omega,val_epsilonO,val_muO,val_c,val_lambdaO,val_kO,val_etaO);
        cout.flush();
    } catch(const InputError& error) { report_input_error(error); }
      catch(const bad_alloc&) { report_input_error(InputError("input.pnu","unable to allocate mesh storage")); }
}
