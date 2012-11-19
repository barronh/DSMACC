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
    spec_ref = readfile(check_path)
    spec_test = readfile(test_path)
    
    abs_tol = 1e-25
    rel_tol = 1e-5
    print 'Species'.ljust(16) + ' Pass'
    checks = []
    for key, refvals in spec_ref.iteritems():
        checkvals = spec_test.get(key, [-9999. for i in refvals])
        diffs = [checkval - refval for checkval, refval in zip(checkvals, refvals)]
        allowables = [abs(refval) * rel_tol + abs_tol for refval in refvals]
        check = all([diff <= allowable for diff, allowable in zip(map(abs, diffs), allowables)])
        checks.append(check)
        print key.ljust(16), check
    
    sum_checks = sum(checks)
    len_checks = len(checks)
    print '-'*25
    print 'Passed'.ljust(16), '%d/%d' % (sum_checks, len_checks)
    
    if len_checks > sum_checks:
        print '!!!!Failed:', len_checks - sum_checks
    
    assert(all(checks))
    return checks


if __name__ == '__main__':
    import sys
    check_path, test_path = sys.argv[1:]
    check(check_path, test_path)
     
