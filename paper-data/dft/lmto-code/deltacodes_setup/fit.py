#!/usr/bin/env python3
import numpy
import subprocess
import re
import sys
import os.path
from scipy.optimize import curve_fit

import warnings
warnings.simplefilter('ignore') # scipy's fitting emits redundant warnings occasionally

# which system to look at?
if len(sys.argv) > 1 :
    extn=sys.argv[1].lower()
else:
    exit()

def bm(x,e0,v0,b0,b1):
    vzox = (v0/x)**(2.0/3.0)
    return e0 + (9.0/16.0)*v0*b0*( ((vzox-1.0)**3.0)*b1 + ((vzox-1.0)**2.0)*(6.0-4.0*vzox) )

# for output in A^3/atom, eV for energies
rydEv = 13.6056919
bohrAngVol = 0.1481845

lsout = subprocess.check_output(['ls'])
names = 'scflog_[0-9]\.[0-9][0-9]' + '\.' +extn +'\\\\n'
scffiles = re.findall(names,str(lsout))  # returns a list of files with scflog[letters,digits,underscores,points].extn

if scffiles != []:

    v = []
    e = []

    for fn in scffiles:
        fn = re.sub('\\\\n','',fn)  # strip off trailing newline from the ls command above
        f = open(fn, 'r')
        
        eryd = 0 
        nbas = 0
        vau = 0
        for line in f:
            # line by line allows subsitution in arb. large files
            res = re.sub('^c.*ehf=(-\d*\.\d*).*','\\1',line)
            if res != line:
                eryd = res

            res = re.sub('.*nbas =(\s*\d*).*','\\1',line)
            if res != line:
                nbas = res

            res = re.sub('.*vol =(\s*\d*\.\d*).*','\\1',line)
            if res != line:
                vau = res

        if eryd != 0 and nbas != 0 and vau !=0:
            v.append(float(vau)*bohrAngVol/float(nbas))
            e.append(float(eryd)*rydEv/float(nbas))
        else:
            sys.stderr.write(fn+' did not converge\n')

    if len(v) > 5:
        vguess = sum(v)/len(v)
        eguess = min(e)
        try:
            (a,b) = curve_fit(bm, v, e, p0=(eguess,vguess,10.0,4.5),method='trf',jac='3-point')  # fitting from scipi
        except RuntimeError:
            sys.stderr.write('BM fitting failed for '+extn+'\n')
        except OptimizeWarning:
            sys.stderr.write('BM fitting unreliable for '+extn+'\n')
        else:
            e0 = a[0]
            v0 = a[1]
            b0 = a[2]
            b1 = a[3]

            erre0 = numpy.sqrt(b[0][0])/e0
            errv0 = numpy.sqrt(b[1][1])/v0
            errb0 = numpy.sqrt(b[2][2])/b0
            errb1 = numpy.sqrt(b[3][3])/b1

            #if errv0 < 1e-2:
            if errv0 < 5e-2:
                # 160.217 eV/GPa
                print('{0:4s}{1:12.6f}{2:12.6f}{3:12.6f}'.format(extn.capitalize(),v0,b0*160.217733,b1))
            else:
                print('{0:4s}{1:30s}{2:12.6f}'.format(extn.capitalize(),' poor convergence: error V0: ',errv0),file=sys.stderr)

    else:
        sys.stderr.write('too few datapoints for '+extn+'\n')

else:
    sys.stderr.write('data missing for '+extn+'\n')
