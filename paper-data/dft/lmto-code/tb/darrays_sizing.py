


def f(realsz, bsizes, psizes):
   blocks = (realsz - 1) / bsizes + 1
   bl_prc = (blocks - 1) / psizes + 1
   gsizes = bl_prc * psizes * bsizes
   return gsizes




for i in range(130):
   print '%5i %5i' %  (i, f(i, 32, 4))
