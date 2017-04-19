from __future__ import print_function
endc = '\033[0m'
def readfile(path):
    lines = open(path, 'r').read().strip().split('\n')
    header = lines.pop(0)
    header = [hv.strip() for hv in header.split('!')[:-1]]
    data = dict([(h, []) for h in header])
    for line in lines:
        for k, v in zip(header, [float(lv) for lv in line.split('!')[:-1]]):
            data[k].append(v)
    return data

def check(check_path, test_path):
    from math import copysign
    spec_ref = readfile(check_path)
    spec_test = readfile(test_path)
    
    abs_tol = 1e-5
    rel_tol = 1e-4
    print('Pass  ' + 'Species'.ljust(16))
    checked = []
    checks = []
    rdiffvals = []
    for key, refvals in spec_ref.items():
        if key in spec_test:
            checkvals = spec_test[key]
            diffs = [checkval - refval for checkval, refval in zip(checkvals, refvals)]
            rdiffs = [0 if c == 0 else d/c for d, c in zip(diffs, checkvals)]
            rdiff = sum(rdiffs)/len(rdiffs)
            allowables = [abs(refval) * rel_tol + abs_tol for refval in refvals]
            tests = [diff <= allowable for diff, allowable in zip([abs(d) for d in diffs], allowables)]
            #if key == 'ANO3':
            #    import pdb; pdb.set_trace()
            check = all(tests)
        else:
            check = 0
            rdiff = 0
        checked.append(key)
        checks.append(check)
        rdiffvals.append(rdiff)
        beginc = '\033[92m'
        if not check:
            beginc = '\033[91m'

        print(beginc + str('PASS' if check else 'FAIL').ljust(5), key.ljust(16) + endc)
        if check == False:
            print(beginc + 'Total , Failed , Pct')
            print(len(refvals), ',', len(refvals) - sum(tests), ',', round((1-sum(tests)/float(len(refvals)))*100, 1), '%' + endc)
            for ti, (r, c, d, a, t) in enumerate(zip(refvals, checkvals, diffs, allowables, tests)):
                if not t:
                    beginc = '\033[91m'
                    print(beginc + '  Fail: %s' % ti + endc)
                else:
                    beginc = '\033[92m'
                    print(beginc + '  Pass: %s' % ti + endc)
                print(beginc + '  Ref    : %s' % r + endc)
                print(beginc + '  Check  : %s' % c + endc)
                print(beginc + '  Diff   : %s' % d + endc)
                print(beginc + '  Allowed: %s' % copysign(a, d) + endc)
                print()
    
    sum_checks = sum(checks)
    len_checks = len(checks)
    counts = ' %5d/%d' % (sum_checks, len_checks)
    failed = len_checks > sum_checks
    
    print('Species MNB   %')
    for key, check, rdiff in zip(checked, checks, rdiffvals):
         if check:
              print('\033[92m' + key + (' %.3e' % rdiff) + (' %.0f%%' % (rdiff*100)) + endc)
         else:
              print('\033[91m' + key + (' %.3e' % rdiff) + (' %.0f%%' % (rdiff*100)) + endc)
    print()
    if failed:
        beginc = '\033[91m'
    print(beginc + '-'*30 + endc)
    print('\033[92m' + 'Passed'.ljust(16) + counts + endc)
    print( beginc + 'Failed'.ljust(16) + ' %5d/%d' % ((len_checks - sum_checks), len_checks) + endc)
    print(''    )
    
    return checks


if __name__ == '__main__':
    import sys
    check_path, test_path = sys.argv[1:]
    check(check_path, test_path)
     
