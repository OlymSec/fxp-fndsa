#!/usr/bin/env python3
"""Build unmodified fndsa in temporary storage; check the protocol in both directions."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',type=Path,required=True);p.add_argument('--cc',default='cc');p.add_argument('--reference',type=Path,default=ROOT/'tests/reference');p.add_argument('--output',type=Path);a=p.parse_args()
    ref=a.reference.resolve();build=a.build.resolve();a.output=a.output or build/'interop';a.output.mkdir(parents=True,exist_ok=True)
    files=[x for x in (ref/'vendor/fndsa_fx32_64').iterdir() if x.is_file()]
    before={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in files}
    rb=a.output/'reference-build';log=a.output/'build.log'
    def run(cmd):
        with log.open('a') as f:subprocess.run(list(map(str,cmd)),cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
    run(['make','-C',ref,'-j4','CC='+a.cc,'BUILD='+str(rb),'all'])
    objects=sorted(x for x in rb.glob('*.o') if x.name!='fndsa.o')
    target=subprocess.check_output([a.cc,'-dumpmachine'],text=True)
    clang='clang' in subprocess.check_output([a.cc,'--version'],text=True).lower()
    flags=['-O2','-std=gnu11','-I'+str(ref/'vendor/fndsa_fx32_64'),'-I'+str(ROOT/'include')]
    if ('aarch64' in target or 'arm64' in target) and not clang:flags+=['-DFNDSA_NEON=0']
    for corrected in [False,True]:
        objs=objects
        if corrected:
            src=(ref/'vendor/fndsa_fx32_64/util_mrm.c').read_text()
            for tag in ['H1','H2']:
                old=f'shake_inject(&sc, "Falcon MRM: {tag}", 2);';assert src.count(old)==1
                src=src.replace(old,f'shake_inject(&sc, "Falcon MRM: {tag}", 14);')
            fixed=a.output/'util_mrm_spec.c';fixed.write_text(src)
            run([a.cc,*flags,'-c',fixed,'-o',a.output/'util_mrm_spec.o'])
            objs=[a.output/'util_mrm_spec.o' if x.name=='util_mrm.o' else x for x in objects]
        binary=a.output/('specified' if corrected else 'unmodified')
        run([a.cc,*flags,'-DCORRECTED_MRM='+str(int(corrected)),ROOT/'tests/test_fndsa_interop.c',build/'libansa_fxp_fndsa.a',*objs,'-pthread','-o',binary])
        result=subprocess.check_output([str(binary)],text=True);(a.output/(binary.name+'.log')).write_text(result);print(result,end='')
    assert all(hashlib.sha256(Path(n).read_bytes()).hexdigest()==h for n,h in before.items())
    (a.output/'result.json').write_text(json.dumps({'status':'PASS','reference_sources_unchanged':True,'reference_sha256':before,'MRM_reference_adjustment':'Only the two domain-tag lengths changed from 2 to 14 in a temporary object; library sources unchanged'},indent=2)+'\n')
if __name__=='__main__':main()
