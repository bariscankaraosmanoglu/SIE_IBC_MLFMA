#!/usr/bin/env python3
"""Bounded IBC port checks; requires AGENTS.md approval before runtime execution."""
import argparse
from array import array
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
UNITS={
 'krylov_solvers':['gmres.cpp','fgmres.cpp','fgmres_in.cpp','krylov.cpp'],
 'krylov_preconditioner':['gmres.cpp','fgmres.cpp','fgmres_in.cpp','pc_solution.cpp','pc_solution_in.cpp','krylov.cpp'],
 'krylov_stage':['gmres.cpp','fgmres.cpp','stage_solution.cpp','krylov.cpp'],
 'block_factorization':['matrixBlock.cpp'],
 'angular_cache':['mlfma_workspace.cpp','sample_locations_weights.cpp','interpolationweights.cpp','input_validation.cpp'],
 'surface_components':['surface_components.cpp'],
 'thread_partitions':['calculate_loadbalance.cpp']}

def fixture():
    points=[];faces=[]
    vertices=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    triangles=[(0,2,4),(2,1,4),(1,3,4),(3,0,4),(2,0,5),(1,2,5),(3,1,5),(0,3,5)]
    for material,center in enumerate((-0.6,0.6)):
        first=len(points)
        points += [(center+0.003*x,0.007+0.003*y,0.011+0.003*z) for x,y,z in vertices]
        faces += [(material,*[first+i for i in f]) for f in triangles]
    lines=['SimulationInputs','model_name ibc_port_reference','Formulations 1',
           'Frequency 1000000000 1000000000 1','source_type 0','source_PW 37 37 1 19 21 3 0',
           'Preconditioners 2','num_mlfmalevels 1','dir_ff 0','num_digits 6','val_approx 0.8',
           'num_interpol 4','val_itermax_pc 10','val_soltol 0.001','dir_cfc 90 0',
           'val_itermax 100','parallel_omp 1 1','## End-Param','ModelParts 2',
           'part0 0.001 0','part1 0.2 0.1',f'POINTS {len(points)}']
    lines += [' '.join(map(str,p)) for p in points]
    lines += [f'CELLS {len(faces)}']+[' '.join(map(str,f)) for f in faces]
    return '\n'.join(lines)+'\n'

def binary(path):
    values=array('d');values.frombytes(path.read_bytes())
    assert len(values)%2==0 and all(math.isfinite(x) for x in values)
    return [complex(values[i],values[i+1]) for i in range(0,len(values),2)]
def relative(a,b):
    assert len(a)==len(b)
    return math.sqrt(math.fsum(abs(x-y)**2 for x,y in zip(a,b))/(math.fsum(abs(x)**2 for x in a) or 1))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir',type=Path,default=ROOT/'src')
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--baseline-source',type=Path)
    parser.add_argument('--sanitize',action='store_true')
    parser.add_argument('--build-only',action='store_true')
    parser.add_argument('--skip-build',action='store_true')
    args=parser.parse_args()
    source=args.source_dir.resolve()
    output=(args.output_dir or Path(tempfile.mkdtemp(prefix='mlfma-ibc-check-',dir='/tmp'))).resolve()
    output.mkdir(parents=True,exist_ok=True)
    flags=['-std=gnu++17','-Wall','-fexceptions','-fopenmp','-I'+str(source)]
    flags += ['-O1','-g','-fsanitize=address,undefined','-fno-omit-frame-pointer'] if args.sanitize else ['-O3']
    env=dict(os.environ,OMP_DYNAMIC='FALSE',OMP_NUM_THREADS='8',
             ASAN_OPTIONS='detect_leaks=0',UBSAN_OPTIONS='halt_on_error=1',TMPDIR=str(output))
    env.pop('OMP_THREAD_LIMIT',None)
    commands=json.loads((output/'commands.json').read_text()) if args.skip_build and (output/'commands.json').exists() else []
    results={}
    def record(command,**kw):
        commands.append({'command':list(map(str,command)),**kw})
        (output/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    def build(name,files,build_flags=flags):
        command=['g++',*build_flags,*map(str,files),'-o',str(output/name)]
        if not args.skip_build:
            with (output/(name+'-build.log')).open('w') as log:
                result=subprocess.run(command,env=env,stdout=log,stderr=log,timeout=300)
            record(command,exit_status=result.returncode,kind='build')
            assert result.returncode==0,(name,output/(name+'-build.log'))
    for name,files in UNITS.items():
        build(name,[ROOT/'tests'/f'{name}.cpp',*[source/f for f in files]])
    all_sources=sorted(p for p in source.glob('*.cpp') if p.name!='main.cpp')
    build('ibc_operators',[ROOT/'tests/ibc_operators.cpp',*all_sources])
    build('run_mlfma',[source/'main.cpp',*all_sources])
    if args.baseline_source:
        assert not args.sanitize,'Guarded original operator reference is normal-build only.'
        baseline=args.baseline_source.resolve()
        text=(baseline/'calculate_Znearmatrix.cpp').read_text()
        before='                                val_imp = val_impedance[tri_media[basistri]];\n                                if(basistri!=-1){'
        after='                                if(basistri!=-1){\n                                    val_imp = val_impedance[tri_media[basistri]];'
        assert before in text
        guarded=output/'baseline-guarded-assembly.cpp';guarded.write_text(text.replace(before,after))
        files=[p for p in baseline.glob('*.cpp') if p.name not in ('main.cpp','calculate_loadbalance.cpp','calculate_Znearmatrix.cpp')]
        build('baseline_operators',[ROOT/'tests/ibc_operators.cpp',guarded,*sorted(files)],
              ['-O3','-std=gnu++17','-fopenmp','-DIBC_LEGACY_REFERENCE','-I'+str(baseline)])
    if args.build_only:
        print('Build-only checks complete; no tests or solver executed:',output);return
    def run(name,binary_name,arguments=(),text=None,case_env=env,expected=0):
        directory=output/name;directory.mkdir(exist_ok=True)
        if text is not None:(directory/'input.pnu').write_text(text)
        command=[str(output/binary_name),*map(str,arguments)]
        with (directory/'run.log').open('w') as log:
            result=subprocess.run(command,cwd=directory,env=case_env,stdout=log,stderr=log,timeout=90)
        record(command,cwd=str(directory),exit_status=result.returncode,kind='runtime')
        assert result.returncode==expected,(name,result.returncode,directory/'run.log')
        results[name]={'exit_status':result.returncode}
        print(name,'PASS',flush=True);return directory
    for name in UNITS:
        if name=='block_factorization':
            for mode in ('empty','first-pivot','later-pivot','repeated','singular','near-singular','nonfinite','complex','scaled-small','scaled-large'):
                run('unit-block-'+mode,name,(mode,))
        else:run('unit-'+name,name)
    text=fixture();(output/'fixture.pnu').write_text(text)
    cases=[('serial',1,1,env),('disabled8',0,8,env)]
    cases += [(f'threads{n}',1,n,env) for n in (2,8,16,32)]
    cases += [('limited32',1,32,dict(env,OMP_THREAD_LIMIT='2')),('dynamic32',1,32,dict(env,OMP_DYNAMIC='TRUE'))]
    names=('near11.bin','near12.bin','near21.bin','near22.bin','near.bin','far_j.bin','far_m.bin','approximate.bin')
    reference=None
    for name,enabled,threads,case_env in cases:
        directory=run(name,'ibc_operators',(enabled,threads),text,case_env)
        if reference is None:reference=directory
        errors={file:relative(binary(reference/file),binary(directory/file)) for file in names}
        assert max(errors.values())<=1e-12,(name,errors)
        results[name]['relative_errors']=errors
    if args.baseline_source:
        directory=run('guarded-original-serial','baseline_operators',(1,1),text)
        errors={file:relative(binary(directory/file),binary(reference/file)) for file in names}
        assert max(errors.values())<=1e-12,errors
        results['guarded-original-serial']['relative_errors']=errors
    # Full solver: all 48 unknowns, all 3 illuminations, existing convergence tolerance.
    matrix=binary(reference/'accurate_matrix.bin');rhs=binary(reference/'rhs.bin');n=len(rhs)
    assert n==48 and len(matrix)==n*n
    for pc in (0,1,2):
        case=text.replace('Preconditioners 2',f'Preconditioners {pc}')
        directory=run('solve-mode'+str(pc),'run_mlfma',text=case,expected=1 if pc==1 else 0)
        rows=[line.split() for line in (directory/'result_solver_status.txt').read_text().splitlines() if line and not line.startswith('#')]
        assert len(rows)==3 and all(r[2]=='converged' and float(r[5])<=float(r[7]) for r in rows)
        if pc==1:
            # Left-preconditioned convergence is not physical convergence. This
            # fixture deliberately exposes the difference; do not tighten input
            # tolerances or hide the solver's correct failure status.
            assert all(r[1]=='physical_residual_not_converged' and float(r[6])>float(r[7]) for r in rows)
        else:
            assert all(r[1]=='converged' and float(r[6])<=float(r[7]) for r in rows)
        coefficients=[complex(*map(float,line.split())) for line in (directory/'result_rwgcoefficient.txt').read_text().splitlines()]
        assert len(coefficients)==n
        residual=[rhs[r]-sum(matrix[r*n+c]*coefficients[c] for c in range(n)) for r in range(n)]
        physical=math.sqrt(sum(abs(v)**2 for v in residual)/sum(abs(v)**2 for v in rhs))
        assert abs(physical-float(rows[-1][6]))<=1e-12,(pc,physical,rows[-1])
        assert (physical>0.001) if pc==1 else (physical<=0.001)
        for file in ('result_currentontris_J.txt','result_currentontris_M.txt','result_cfcontris.txt','result_farefield.txt'):
            values=[float(x) for line in (directory/file).read_text().splitlines() if not line.startswith('%') for x in line.split()]
            assert values and all(math.isfinite(x) for x in values),file
        results['solve-mode'+str(pc)]['independent_physical_residual']=physical
    invalid={
      'formulation0':('Formulations 1','Formulations 0'),
      'source1':('source_type 0','source_type 1'),
      'zero-impedance':('part0 0.001 0','part0 0 0'),
      'nan-impedance':('part0 0.001 0','part0 nan 0'),
      'infinite-impedance':('part0 0.001 0','part0 1 inf'),
      'reciprocal-overflow':('part0 0.001 0','part0 1e-320 0'),
      'part-fields':('part0 0.001 0','part0 0.001'),
      'frequency-sweep':('Frequency 1000000000 1000000000 1','Frequency 1000000000 1000000000 2'),
      'threads0':('parallel_omp 1 1','parallel_omp 1 0'),
      'preconditioner3':('Preconditioners 2','Preconditioners 3'),
      'tolerance0':('val_soltol 0.001','val_soltol 0'),
      'grid7':('dir_ff 0','dir_ff 7'),
      'invalid-part':('0 0 2 4','2 0 2 4'),
      'invalid-node':('0 0 2 4','0 0 2 12'),
      'truncated-mesh':('CELLS 16','CELLS 99999'),
      'stencil-too-wide':('num_interpol 4','num_interpol 999')}
    for name,(a,b) in invalid.items():
        assert a in text,a
        directory=run('invalid-'+name,'run_mlfma',text=text.replace(a,b),expected=1)
        assert 'Input error:' in (directory/'run.log').read_text()
    # A zero approximate far operator is a supported near-only inner solve.
    directory=run('solve-near-only','run_mlfma',text=text.replace('val_approx 0.8','val_approx 0'))
    rows=[line.split() for line in (directory/'result_solver_status.txt').read_text().splitlines() if line and not line.startswith('#')]
    assert len(rows)==3 and all(r[1]=='converged' and float(r[6])<=0.001 for r in rows)
    if args.sanitize:
        run('leak-cleanup','ibc_operators',(1,2),text,dict(env,ASAN_OPTIONS='detect_leaks=1'))
    (output/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    print('All bounded IBC checks passed:',output)
if __name__=='__main__':main()
