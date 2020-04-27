
import re
from small_utils import loadfile

lapack_path = '/home/dmt/dev/lapack-3.4.1/SRC'
scalapack_path = '/home/dmt/dev/scalapack/scalapack-2.0.2/SRC'

for m in ' ', 'p':
   for t in 'dsy', 'zhe':
      for p in 'ev', 'gv':
         for a in 'x', 'r', 'd':
            r=m+t+p+a
            r = r.strip()
            path = lapack_path if m == ' ' else scalapack_path

            try:
               s = loadfile(path+'/'+r+'.f')
               s = re.sub(r'^[*Cc!].*(?m)','',s)
               pr = re.findall(r'subroutine \s+ (%s \s* \( .*? \)) (?six)'%r,s)[0]
            except:
               pr = ''


            pr = re.sub(r'\n \s{5} [^\s](?ix)',' ',pr)
            pr = re.sub(r'\s+(?ix)',' ',pr)
            print '%s: %s' % (r, pr)


