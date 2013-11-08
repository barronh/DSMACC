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
        print str(check).ljust(5), key.ljust(16)
        if check == False:
            ti = 0
            print 'Faild , Total , Pct', len(refvals) - sum(tests), ',', len(refvals), ',', round((1-sum(tests)/float(len(refvals)))*100, 1), '%'
            for r, c, d, a, t in zip(refvals, checkvals, diffs, allowables, tests):
                if not t:
                    ti += 1
                    print '  Fail:', ti
                    print '  Ref    :', r
                    print '  Check  :', c
                    print '  Diff   :', d
                    print '  Allowed:', copysign(a, d)
                    print
    
    sum_checks = sum(checks)
    len_checks = len(checks)
    print '-'*25
    print 'Passed'.ljust(16), '%d/%d' % (sum_checks, len_checks)
    
    if len_checks > sum_checks:
        print '!!!!Failed:', len_checks - sum_checks
    else:
        print ''    
    print ''    
    
    return checks


if __name__ == '__main__':
    import sys
    check_path, test_path = sys.argv[1:]
    check(check_path, test_path)
     
