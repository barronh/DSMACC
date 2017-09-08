#!/usr/bin/env python
with open('../src/dsmacc_Parameters.f90') as paramf: 
  lines = [l for l in paramf.read().split('\n') if ' NVAR =' in l or ' NSPEC =' in l or ' NREACT =' in l]
  data = '\n'.join(lines).replace('  INTEGER, PARAMETER :: ', '')
  exec(data)

with open('kpp.pyf.tmp') as intmp, open('kpp.pyf', 'w') as outf: 
  outf.write(intmp.read() % globals())