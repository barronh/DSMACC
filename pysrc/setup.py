from Cython.Distutils import build_ext
from distutils.extension import Extension
from distutils.core import setup
import numpy

def configuration(parent_package = '', top_path = None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('.', parent_package, top_path)
    config.add_extension('kpp', sources = ['kpp.pyf', 'pyint.f90', 'pyglob.f90', 'pyrate.f90', 'pymon.f90', '../src/dsmacc_Global.f90', '../src/dsmacc_Parameters.f90'], extra_objects = '../src/dsmacc_Function.o ../src/dsmacc_Global.o ../src/dsmacc_Initialize.o ../src/dsmacc_Integrator.o ../src/dsmacc_Jacobian.o ../src/dsmacc_JacobianSP.o ../src/dsmacc_LinearAlgebra.o ../src/dsmacc_Main.o ../src/dsmacc_Model.o ../src/dsmacc_Monitor.o ../src/dsmacc_Parameters.o ../src/dsmacc_Precision.o ../src/dsmacc_Rates.o ../src/dsmacc_Stoichiom.o ../src/dsmacc_StoichiomSP.o ../src/dsmacc_Util.o  ../tuv_new/TUV.o ../tuv_new/functs.o ../tuv_new/grids.o ../tuv_new/la_srb.o ../tuv_new/numer.o ../tuv_new/odo3.o ../tuv_new/odrl.o ../tuv_new/orbit.o ../tuv_new/qys.o ../tuv_new/rdetfl.o ../tuv_new/rdinp.o ../tuv_new/rdxs.o ../tuv_new/rtrans.o ../tuv_new/rxn.o ../tuv_new/savout.o ../tuv_new/setaer.o ../tuv_new/setalb.o ../tuv_new/setcld.o ../tuv_new/setno2.o ../tuv_new/seto2.o ../tuv_new/setsnw.o ../tuv_new/setso2.o ../tuv_new/sphers.o ../tuv_new/swbiol.o ../tuv_new/swchem.o ../tuv_new/swphys.o ../tuv_new/vpair.o ../tuv_new/vpo3.o ../tuv_new/vptmp.o ../tuv_new/wshift.o ../src/constants.o'.split(), include_dirs = '../src', module_dirs = '../src')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    cfg = configuration(top_path = '')
    setup(**cfg.todict())
