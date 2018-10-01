from glob import glob
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

def makekpp(mech, path):
    import os.path
    parampath = '../{}/src/dsmacc_Parameters.f90'.format(path)
    outpath = '{}.pyf'.format(mech)
    if os.path.exists(outpath):
        param_mtime = os.path.getmtime(parampath)
        out_mtime = os.path.getctime(outpath)
        if out_mtime > param_mtime: return
    with open(parampath) as paramf: 
        lines = [l for l in paramf.read().split('\n') if ' NVAR =' in l or ' NSPEC =' in l or ' NREACT =' in l]
    params = {'mech': mech}
    data = '\n'.join(lines).replace('  INTEGER, PARAMETER :: ', '')
    exec(data, None, params)

    with open('kpp.pyf.tmp') as intmp, open(outpath, 'w') as outf: 
      outf.write(intmp.read() % params)

_sources = ['../{}/src/dsmacc_Precision.f90', '../{}/src/dsmacc_Parameters.f90', '../{}/src/dsmacc_Function.f90', '../{}/src/dsmacc_Global.f90', '../{}/src/dsmacc_JacobianSP.f90', '../{}/src/dsmacc_Jacobian.f90', '../{}/src/dsmacc_LinearAlgebra.f90', '../{}/src/dsmacc_Rates.f90', '../{}/src/dsmacc_Integrator.f90', '../{}/src/dsmacc_Monitor.f90', '../{}/src/dsmacc_Util.f90', '../{}/src/dsmacc_Initialize.f90', '../{}/src/dsmacc_Main.f90', '../{}/src/dsmacc_StoichiomSP.f90', '../{}/src/dsmacc_Stoichiom.f90', '../{}/src/dsmacc_Model.f90', 'pyint.f90', 'pyglob.f90', 'pyrate.f90', 'pymon.f90']

def addext(mech, path):
    makekpp(mech, path)
    sources = ['{}.pyf'.format(mech)] + [s.format(path) for s in _sources]
    mod = Extension('dsmacc.{}'.format(mech),
                sources = sources,
                library_dirs = ['../tuv_new/', '../UCI_fastJX72e/'],
                libraries = ['tuv', 'fastjx'],
                include_dirs = ['.', '../src'],
                extra_f90_compile_args = ['-cpp', '-O0'])
    return mod
    
crimod = addext('cri', 'cri')
geoschemmod = addext('geoschem', 'geoschem')

iso = Extension('dsmacc._isoropia',
            sources = ['../ISOROPIA/isoropia_dsmacc.F90'],
            library_dirs = ['../ISOROPIA/'],
            libraries = ['Isoropia'],
            include_dirs = ['.', '../ISOROPIA/'],
            extra_f90_compile_args = ['-cpp', '-O0'])

            
if __name__ == '__main__':
    setup (name = 'dsmacc',
           version = '1.0',
           description = 'Python interface to KPP chemistry solvers',
           author = 'Barron H. Henderson',
           author_email = 'barronh@gmail.com',
           url = 'https://docs.python.org/extending/building',
           package_dir = {'': 'lib'},
           packages = ['dsmacc'],
           ext_modules = [crimod, geoschemmod, iso])
