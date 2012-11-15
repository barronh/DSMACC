def readfile(path):
    header, data = file(path, 'r').read().strip().split('\n')
    
    header = map(str.strip, header.split('!')[:-1])
    data = map(float, data.split('!')[:-1])
    return dict(zip(header, data))

spec_ref = readfile('Spec_1.dat.check')
spec_test = readfile('../Spec_1.dat')

abs_tol = 1e-25
rel_tol = 1e-5
print 'Species'.ljust(16) + ' Pass'
checks = []
for key, refval in spec_ref.iteritems():
    checkval = spec_test[key]
    diff = checkval - refval
    pctdiff = diff / refval
    
    allowable = abs(refval) * rel_tol + abs_tol
    check = abs(diff) <= allowable
    checks.append(check)
    print key.ljust(16), check

sum_checks = sum(checks)
len_checks = len(checks)
print '-'*25
print 'Passed'.ljust(16), '%d/%d' % (sum_checks, len_checks)

if len_checks > sum_checks:
    print '!!!!Failed:', len_checks - sum_checks

assert(all(checks))
