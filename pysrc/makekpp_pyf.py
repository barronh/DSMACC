#!/usr/bin/env python
with open('../src/dsmacc_Parameters.f90') as paramf:
    lines = [
        _l for _l in paramf.read().split('\n')
        if ' NVAR =' in _l or ' NSPEC =' in _l or ' NREACT =' in _l
    ]
    data = '\n'.join(lines).replace('  INTEGER, PARAMETER :: ', '')
    exec(data)

with open('kpp.pyf.tmp') as intmp, open('kpp.pyf', 'w') as outf:
    outf.write(intmp.read() % globals())
