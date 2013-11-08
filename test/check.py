endc = '\033[0m'
def readfile(path):
    lines = file(path, 'r').read().strip().split('\n')
    header = lines.pop(0)
    header = map(str.strip, header.split('!')[:-1])
    data = dict([(h, []) for h in header])
    for line in lines:
        for k, v in zip(header, map(float, line.split('!')[:-1])):
            data[k].append(v)
    return data

def check(check_path, test_path):
    from math import copysign
    spec_ref = readfile(check_path)
    spec_test = readfile(test_path)
    
    abs_tol = 1e-5
    rel_tol = 1e-4
    print 'Pass  ' + 'Species'.ljust(16)
    checks = []
    for key, refvals in spec_ref.iteritems():
        if key in spec_test:
            checkvals = spec_test[key]
            diffs = [checkval - refval for checkval, refval in zip(checkvals, refvals)]
            allowables = [abs(refval) * rel_tol + abs_tol for refval in refvals]
            tests = [diff <= allowable for diff, allowable in zip(map(abs, diffs), allowables)]
            check = all(tests)
        else:
            check = 0
        checks.append(check)
        beginc = '\033[92m'
        if not check:
            beginc = '\033[91m'

        print beginc + str('PASS' if check else 'FAIL').ljust(5), key.ljust(16) + endc
        if check == False:
            ti = 0
            print beginc + 'Total , Failed , Pct', len(refvals), ',', len(refvals) - sum(tests), ',', round((1-sum(tests)/float(len(refvals)))*100, 1), '%' + endc
            for r, c, d, a, t in zip(refvals, checkvals, diffs, allowables, tests):
                if not t:
                    ti += 1
                    print beginc + '  Fail: %s' % ti + endc
                    print beginc + '  Ref    : %s' % r + endc
                    print beginc + '  Check  : %s' % c + endc
                    print beginc + '  Diff   : %s' % d + endc
                    print beginc + '  Allowed: %s' % copysign(a, d) + endc
                    print
    
    sum_checks = sum(checks)
    len_checks = len(checks)
    counts = ' %5d/%d' % (sum_checks, len_checks)
    failed = len_checks > sum_checks
    
    
    
    if failed:
        beginc = '\033[91m'
    print beginc + '-'*30 + endc
    print '\033[92m' + 'Passed'.ljust(16) + counts + endc
    print  beginc + 'Failed'.ljust(16) + ' %5d/%d' % ((len_checks - sum_checks), len_checks) + endc
    print ''    
    
    return checks


if __name__ == '__main__':
    import sys
    check_path, test_path = sys.argv[1:]
    check(check_path, test_path)
     
