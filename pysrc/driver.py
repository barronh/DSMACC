"""
This file is an example of how to use dsmacc in python
"""
from dsmacc import model, dynenv

# mod = model();
# alternatively, use a file ('test.dat') to make a dynamic
# physical environment. Add variables as you see fit.
# 
import pandas as pd
import io
envdata = pd.read_csv(io.StringIO("""JDAY_GMT,TEMP,PBL
2000172.2500,298.,100.
2000172.2916,290.,100.
2000172.3333,275.,200.
2000172.3750,260.,500.
2000172.4166,245.,1000.
2000172.4583,230.,1500.
2000172.5000,215.,1500.
2000172.5416,200.,1500.
2000172.5833,190.,500.
""")).to_dict(orient = 'list')

emissdata = pd.read_csv(io.StringIO("""JDAY_GMT,NO
2000172.2500,4e-14
2000172.2916,4e-14
2000172.3333,9e-14
2000172.3750,13e-14
2000172.4166,20e-14
2000172.4583,13e-14
2000172.5000,9e-14
2000172.5416,4e-14
2000172.5833,4e-14
""")).to_dict(orient = 'list')

mod = dynenv(envdata = envdata, emissdata = emissdata, bkgdata = {})
mod.run(2000172.25, 8, 180, conc_ppb = dict(O3 = 30., C5H8 = 850., NO = 128.), globvar = dict(LAT = 30., LON = 0, TEMP = 298., PRESS = 101325.))
