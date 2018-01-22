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
data.eval('C5H8 = C5H8 * 100.', inplace = True)

haspbl = 'PBL' in data.columns
hastemp = 'TEMP' in data.columns
if haspbl and hastemp:
    pos = [.1, .1, .55, .7]
    anchor = (.75, 1.05)
elif haspbl or hastemp:
    pos = [.1, .1, .7, .7]
    anchor = (.65, 1.05)
pax = data.plot(x = 'JDAY', y = ['O3', 'C5H8'], linestyle = '-', legend=False, color = ['k', 'g'])
pax.set_position(pos)
sax = pax.twinx()
sax.set_position(pos)
patches, labels = pax.get_legend_handles_labels()
data.plot(x = 'JDAY', y = 'NOx NXOYN HONO HNO3 PANS HO2NO2 NTR NOA NA'.split(), secondary_y = False, ax = sax, stacked = True, legend = False, linestyle = '--')
patchesnax, labelsnax = sax.get_legend_handles_labels()
patches = patches + patchesnax
labels = labels + labelsnax 
offset = 60
if haspbl:
    pblax = pax.twinx()
    s = pblax.spines['right']
    s.set_position(('outward', offset))
    offset += 80
    s.set_color('b')
    pblax.set_position(pos)
    data.plot(x = 'JDAY', y = ['PBL'], linestyle = '-', legend=False, ax = pblax, color = 'b')
    patchesnax, labelsnax = pblax.get_legend_handles_labels()
    patches = patches + patchesnax
    labels = labels + labelsnax 
    pblax.set_ylabel('PBL', color = 'b')
if hastemp:
    tempax = pax.twinx()
    s = tempax.spines['right']
    s.set_position(('outward', offset))
    s.set_color('r')
    tempax.set_position(pos)
    data.plot(x = 'JDAY', y = ['TEMP'], linestyle = '-', legend=False, ax = tempax, color = 'r')
    patchesnax, labelsnax = tempax.get_legend_handles_labels()
    patches = patches + patchesnax
    labels = labels + labelsnax 
    tempax.set_ylabel('TEMP', color = 'r')

pax.set_ylabel('O3 ppb, 50x C5H8 ppb')
sax.set_ylabel('$\sum_i{NOy_i}$')
pax.legend(patches, labels, loc = 'lower center', bbox_to_anchor = anchor, ncol = 5, handlelength = None)
pax.figure.savefig('cri_isop.png')
