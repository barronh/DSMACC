from glob import glob
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

def makekpp():
    with open('../src/dsmacc_Parameters.f90') as paramf: 
        lines = [l for l in paramf.read().split('\n') if ' NVAR =' in l or ' NSPEC =' in l or ' NREACT =' in l]
    params = {}
    data = '\n'.join(lines).replace('  INTEGER, PARAMETER :: ', '')
    exec(data, None, params)

    with open('kpp.pyf.tmp') as intmp, open('kpp.pyf', 'w') as outf: 
      outf.write(intmp.read() % params)

makekpp()

kppmod = Extension('dsmacc.kpp',
                sources = ['kpp.pyf', '../src/dsmacc_Precision.f90', '../src/dsmacc_Parameters.f90', '../src/dsmacc_Function.f90', '../src/dsmacc_Global.f90', '../src/dsmacc_JacobianSP.f90', '../src/dsmacc_Jacobian.f90', '../src/dsmacc_LinearAlgebra.f90', '../src/dsmacc_Rates.f90', '../src/dsmacc_Integrator.f90', '../src/dsmacc_Monitor.f90', '../src/dsmacc_Util.f90', '../src/dsmacc_Initialize.f90', '../src/dsmacc_Main.f90', '../src/dsmacc_StoichiomSP.f90', '../src/dsmacc_Stoichiom.f90', '../src/dsmacc_Model.f90', 'pyint.f90', 'pyglob.f90', 'pyrate.f90', 'pymon.f90'],
                library_dirs = ['../tuv_new/', '../UCI_fastJX72e/'],
                libraries = ['tuv', 'fastjx'],
                include_dirs = ['.', '../src'],
                extra_f90_compile_args = ['-cpp', '-O2'])
if __name__ == '__main__':
    setup (name = 'dsmacc',
           version = '1.0',
           description = 'Python interface to KPP chemistry solvers',
           author = 'Barron H. Henderson',
           author_email = 'barronh@gmail.com',
           url = 'https://docs.python.org/extending/building',
           package_dir = {'': 'lib'},
           packages = ['dsmacc'],
           ext_modules = [kppmod])
