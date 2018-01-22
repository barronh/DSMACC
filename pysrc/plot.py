import pandas as pd;
data = pd.read_csv('output.dat', delimiter = ',')
data.eval('JDAY=time/3600./24.', inplace = True)
data.eval('NOx=NO+NO2', inplace = True)
ntrs = [k for k in data.columns if 'NO3' in k and k not in ('NO3', 'HNO3')]
data.eval('NTR1 = %s' % ' + '.join(ntrs[:25]), inplace = True)
data.eval('NTR2 = %s' % ' + '.join(ntrs[25:]), inplace = True)
data.eval('NTR = NTR1 + NTR2', inplace = True)
data.eval('N2O5N = 2*N2O5', inplace = True)
data.eval('NXOYN = N2O5N + NO3', inplace = True)
pans = [k for k in data.columns if 'PAN' in k] + 'PPN PHAN PBENZ IPBENZ'.split(' ')
data.eval('PANS = %s' % ' + '.join(pans), inplace = True)
if 'PBLH' in data.columns:
    ax = data.plot(x = 'JDAY', y = ['O3', 'C5H8', 'PBLH'], linestyle = '--')
else:
    ax = data.plot(x = 'JDAY', y = ['O3', 'C5H8'], linestyle = '--')
ax.set_ylabel('O3, C5H8 ppb')
ax = data.plot(x = 'JDAY', y = 'NOx NXOYN HONO HNO3 PANS HO2NO2 NTR NOA NA'.split(), secondary_y = True, ax = ax, stacked = True)
ax.set_ylabel('$\sum_i{NOy_i}$')
ax.figure.savefig('cri_isop.png')
