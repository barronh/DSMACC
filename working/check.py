"""import xml.etree.ElementTree
e = xml.etree.ElementTree.parse('mcm_subset.xml').getroot()
species = [s for s in e[0] if 'species_number' in s.keys() and len(s) > 0]

smiles = dict()
for s in species:
    smiles[s.get('species_name')] = [i for i in s[0].itertext()][0]
    
    
Nspecies = {k: v for k, v in smiles.items() if v.count('N') > 0}

PANspecies = {k: v for k, v in Nspecies.items() if 'O=N(=O)OO' in v or 'OON(=O)=O' in v or 'OO[N+](=O)[O-]' in v}
NTRspecies = {k: v for k, v in Nspecies.items() if ('O[N+](=O)[O-]' in v or 'O=N(=O)O' in v or 'ON(=O)=O' in v or 'O[N+](=O)[O-])' in v) and k not in PANspecies}

print({k: v for k, v in Nspecies.items() if not k in PANspecies and not k in NTRspecies})
print('N', len(Nspecies))
print('PAN', len(PANspecies))
print('PAN = ' + ' + '.join([k for k in PANspecies.keys()]))
print('NTR', len(NTRspecies))
print('NTR = ' + ' + '.join([k for k in NTRspecies.keys()]))
"""

from matplotlib import use; use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
rate = pd.read_csv('Rate_timeseries_1.dat', delimiter = '!')
rate.rename(columns = lambda x: x.strip(), inplace = True)
conc = pd.read_csv('Spec_timeseries_1.dat', delimiter = '!')
conc.rename(columns = lambda x: x.strip(), inplace = True)
rate.eval('JDAY_LST = JDAY_GMT + (LON_degE / 15.)/24 - 2006000', inplace = True)
for key in conc.columns:
    if key not in 'TIME JDAY_GMT LAT_degN LON_degE PRESS_hPa TEMP_K THETA H2O CFACTOR M N2 O2 JNO2FACT JO1DFACT'.split() and not 'Unnamed' in key:
        conc.eval('%s = %s / CFACTOR * 1e9' % (key, key), inplace = True)
conc.eval('JDAY_LST = JDAY_GMT + (LON_degE / 15.)/24 - 2006000', inplace = True)
conc.eval('NOx = NO2 + NO', inplace = True)
conc.eval('RN = NOx / HNO3', inplace = True)

plt.close()
ax = rate.plot(x = 'JDAY_LST', y = ['NO2 --> O + NO', 'NO3 --> NO', 'NO3 --> O + NO2', 'H2O2 --> 2 OH'])
rate.plot(x = 'JDAY_LST', y = ['THETA'], secondary_y = True, ax = ax)
plt.savefig('rates.png')

plt.close()
conc.plot(x = 'JDAY_LST', y = 'NO NO2 NO3 HNO3'.split(), stacked = True)
plt.ylabel('cumulative ppbv')
plt.savefig('concs.png')
