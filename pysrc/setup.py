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
                sources = ['kpp.pyf', 'pyint.f90', 'pyglob.f90', 'pyrate.f90', 'pymon.f90'] + glob('../src/dsmacc_*.f90'),
                library_dirs = ['../tuv_new/', '../UCI_fastJX72e/'],
                libraries = ['tuv', 'fastjx'],
                include_dirs = ['.', '../src'],
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
           ext_modules = [kppmod])

# def configuration(parent_package='', top_path='.'):
#     from numpy.distutils.misc_util import Configuration
#     config = Configuration('dsmacc',parent_package, top_path, )
#     config.add_extension('kpp',
#                 sources = ['kpp.pyf', 'pyint.f90', 'pyglob.f90', 'pyrate.f90', 'pymon.f90'] + glob('../src/dsmacc_*.f90'),
#                 library_dirs = ['../tuv_new/', '../UCI_fastJX72e/'],
#                 libraries = ['tuv', 'fastjx'],
#                 include_dirs = ['.', '../src'],
#                 extra_f90_compile_args = ['-cpp', '-O0'])
#     return config
# 
# if __name__ == '__main__':
#     from numpy.distutils.core import setup
#     setup(**configuration(top_path='.').todict())