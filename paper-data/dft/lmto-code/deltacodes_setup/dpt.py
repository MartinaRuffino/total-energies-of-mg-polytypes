#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, re
#from eosfit import BM as bmfit
from numpy import empty, array, linspace, zeros, inf
from scipy.optimize import curve_fit
from colorsys import hsv_to_rgb


#import matplotlib
#matplotlib.use('Agg')
#from pylab import plot, savefig, clf


rydEv       = 13.6056919
bohrAngVol  =  0.1481845
echarge     =  1.60217733e-19

pt = '''
 H                                                  He
 Li Be                               B  C  N  O  F  Ne
 Na Mg                               Al Si P  S  Cl Ar
 K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
 Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
 Cs Ba    Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po    Rn

                                                 Lu
'''

ptl = pt.split()

def loadfile(f, mode='r'):
    fd = open(f, mode)
    s = fd.read()
    fd.close()
    return s

def writefile(f, s, mode='w'):
    if (f != sys.stdout):
        fd = open(f, mode)
    else:
        fd = sys.stdout
    fd.write(s)
    if (fd != sys.stdout):
        fd.close()
    return 0


e_ptr = re.compile(r'(?xms) ^c .*? ehk= \s* ([-0-9.]+)')
v_ptr = re.compile(r'(?xs) Cell \s+ vol \s* = \s* ([0-9.]+)')
n_ptr = re.compile(r'(?xs) nbas \s* = \s* (\d+)')
tbl3_ptr = re.compile(r'(?xm)^(\w{1,2})\s+(-?\d*\.?\d*)\s+(-?\d*\.?\d*)\s+(-?\d*\.?\d*)')


def ptcrd(e):
    e = e.ljust(2)
    i = pt.find(e)
    p = pt[pt.find('H'):i].count('\n')
    g = (len(pt[pt[:i].rfind('\n'):i])-1)//3

    return g,p



# these svg funs shall be bundled into a class at some point
def svg(w,h):
    return '''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="{w}cm" height="{h}cm" version="1.1" xmlns="http://www.w3.org/2000/svg">
<desc>Errs</desc>'''.format(h=h,w=w)


def rect(x,y,w,h,c):
    return '<rect x="{x}cm" y="{y}cm" width="{w}cm" height="{h}cm" fill="{c}" stroke="black"/>'.format(x=x,y=y,w=w,h=h,c=c)

def circle(x,y,r,fill='none',stroke='none'):
    return '<circle cx="{x}cm" cy="{y}cm" r="{r}cm" fill="{fill}" stroke="{stroke}"/>'.format(x=x,y=y,r=r,fill=fill,stroke=stroke)

def text(t,x,y,ff,fs):
    return '<text x="{x}cm" y="{y}cm" font-size="{fs}" font-family="{ff}">{t}</text>'.format(t=t,x=x,y=y,ff=ff,fs=fs) #font-family="{ff}"

def hsv_to_html(h,s,v):
    return '#'+''.join(['%02x'%(int(i*255)) for i in hsv_to_rgb(h,s,v)])


def ptdraw(vals):
    sz = 2,2

    gsz = 18*sz[0],8*sz[1]

    s = svg(*gsz)

    #s += rect(0,0,gsz[0],gsz[1],'none')


    for e in pt.split():
        c = list(ptcrd(e))
        v = vals[e]

        if (v != None):
            h = min(max(0,(3.5-v)/9.),1./3.)
            fill = hsv_to_html(h,0.5,1)
            stroke = 'none'
            v_str = '%.1f'%v
        else:
            stroke = hsv_to_html(0,0.5,1)
            fill = 'none'
            v_str = 'n/a'

        s += '\n'
        s += '\n' + rect(c[0]*sz[0], c[1]*sz[1], sz[0], sz[1], 'none')
        s += '\n' + text(e, (c[0]+0.1)*sz[0], (c[1]+0.4)*sz[1], 'Futura, sans-serif', '20pt')
        s += '\n' + text(v_str, (c[0]+0.1)*sz[0], (c[1]+0.9)*sz[1], 'Futura, sans-serif', '10pt')
        s += '\n' + circle((c[0]+0.75)*sz[0], (c[1]+0.75)*sz[1], 0.4, fill, stroke)

    s += '\n</svg>'

    return s

def parse_dlt(s):
    de = {}
    for e in pt.split(): de[e] = None

    vs = tbl3_ptr.findall(s)
    for v in vs:
        e = v[0]
        de[e] = float(v[1])

    return de

def get_vne(s):
    v = float(v_ptr.search(s).group(1))
    n = int(n_ptr.search(s).group(1))
    #c sc=1 ehf=-4.6273341 ehk=-4.6273368
    e = float(e_ptr.search(s).group(1))

    return v,n,e

# from Jerome
def bm(x,e0,v0,b0,b1):
    vzox = (v0/x)**(2/3.)
    return e0 + (9./16.)*v0*b0*(((vzox-1)**3)*b1 + ((vzox-1)**2)*(6-4*vzox))

# from Jerome
def bmfit(ve):

    assert(ve.shape[0] > 4)

    v, e = ve[:,0], ve[:,1]

    vguess = sum(v)/len(v)
    eguess = min(e)

    try:
        (a,b) = curve_fit(bm, v, e, p0=(eguess,vguess,10.0,4.5), method='trf', jac='3-point')  # fitting from scipi
        e0, v0, b0, b1 = tuple(a[:4])
        errv0 = b[1][1]**0.5/v0
    except RuntimeError:
        sys.stderr.write('\nBM fitting failed\n')
        v0, b0, b1, errv0 = 0, 0, 0, inf
    except OptimizeWarning:
        sys.stderr.write('\nBM fitting unreliable\n')
        v0, b0, b1, errv0 = 0, 0, 0, inf

    return v0, b0, b1, errv0


#def get_vbp(flns, el=''):
def get_vbp(flns, el):
    n = len(flns)

# this is a bit too much, el is now mandatory argument
    #if el == '':
        ## try to guess the element name from the common patterns in the filenames
        #fln0 = array(list(flns[0]))
        #m = fln0 == fln0
        #for i in range(1,n):
            #m = (fln0 == array(list(flns[i]))) & m
        #comm = flns[0][:m.tolist().index(False)]
        #comm = comm.split('/')[-1].strip('-_ ')
        #els = re.findall(r'\w+',comm)
        #for e in els:
            #if e.title() in ptl:
                #el = e
                #break

    ve = empty((n,2))

    for i,fln in enumerate(flns):
        s = loadfile(fln)
        vol,nat,en = get_vne(s)

        ve[i,:] = array([vol*bohrAngVol, en*rydEv])/nat

    if el.lower() == 'he':
        print(ve*nat/array([bohrAngVol, rydEv]))

    #x = ve[:,1]
    #if not (x[0] > x[1] and x[1] > x[2] and x[2] > x[3] and x[3] < x[4] and x[4] < x[5] and x[5] < x[6]):
        #sys.stderr.write(' likely not converged\n')
    #clf()
    #plot(ve[:,0],ve[:,1])
    #savefig('plot-'+el.lower()+'.pdf')

    volume, bulk_modulus, bulk_deriv, rs = bmfit(ve)

    rsp = ''
    if rs != inf:
        rsp = '%s %.5f %.5f %.3f' % (el.title().ljust(2), volume, (bulk_modulus * echarge * 1.0e21), bulk_deriv)
    if rsp == '':
        sys.stderr.write('\nfailed el ' + el + '\n')
        sys.stderr.flush()

    return rsp

def elscl(el, ptr):
    sys.stdout.write(' '+el)

    el = el.lower()
    flns = list([ptr.format(e=el,s='%4.2f'%scl) for scl in linspace(0.94,1.06,7)])

    for fln in flns:
        if not os.path.exists(fln):
            sys.stderr.write('\nfile %s not found\n' % fln)
            return ''
        if e_ptr.search(loadfile(fln)) == None:
            sys.stderr.write('\nfile %s not converged\n' % fln)
            return ''

    s = get_vbp(flns, el)

    sys.stdout.flush()
    sys.stderr.flush()

    return s

class elscl_t(object):
    def __init__(self, ptr):
        self.ptr = ptr
    def __call__(self, el):
        return elscl(el, self.ptr)

def read_summary(f):
    d = tbl3_ptr.findall(loadfile(f))
    return dict([(l[0],array(l[1:],dtype=float)) for l in d])

def dcdiff(fls):
# Birch–Murnaghan equation of state
# E(v) = 9/16 * v b ((v^(2/3) - 1)^3 b' + (v^(2/3) - 1)^2 (6 - 4 v^(2/3)))
# expanded in powers of v

# diff abs = 1000 * (1/((1.06 - 0.94) * (v2 + v1)/2) *
#    sum_ij=0,3 ((e2_i - e1_i) * (e2_j - e1_j)) / (1 - (i+j)*2/3) * ((v2 + v1)/2)**(1 - (i+j)*2/3) * (1.06**(1 - (i+j)*2/3) - 0.94**(1 - (i+j)*2/3))
# )**(1/2)

# diff rel = 100 * (
#    (sum_ij=0,3 ((e2_i - e1_i) * (e2_j - e1_j)) / (1 - (i+j)*2/3) * ((v2 + v1)/2)**(1 - (i+j)*2/3) * (1.06**(1 - (i+j)*2/3) - 0.94**(1 - (i+j)*2/3)))/
#    (sum_ij=0,3 ((e2_i+e1_i)/2 * (e2_j+e1_j)/2) / (1 - (i+j)*2/3) * ((v2 + v1)/2)**(1 - (i+j)*2/3) * (1.06**(1 - (i+j)*2/3) - 0.94**(1 - (i+j)*2/3)))
# )**(1/2)

    assert(len(fls) == 2)

    #r = list([loadtxt(fl, dtype={'names': ('element', 'V0', 'B0', 'BP'), 'formats': ('S2', float, float, float)}) for fl in fls])
    r = list([read_summary(fl) for fl in fls])

    #ce = sorted(r[0].keys() & r[1].keys()) # common elements
    ce = sorted(set(r[0].keys()) & set(r[1].keys())) # add set() for py2 compatibility

    ne = len(ce)
    fac =  1.0e9 / echarge / 1.0e30

    v = [0,0]
    b = [0,0]
    ev = [0,0]

    for i in 0,1:
        vl = array([r[i][el] for el in ce])
        v[i] = vl[:,0]
        b[i] = vl[:,1]*fac
        d = vl[:,2]

        cv = 9./16.* b[i] * v[i]
        ev[i] = [cv*(6-d), cv*v[i]**(2/3.)*(3*d-16), cv*v[i]**(4/3.)*(14-3*d), cv*v[i]**2*(d-4)]

    scm = 0.94
    scx = 1.06

    vav = (v[1] + v[0])*0.5

    np = len(ev[0])
    em = list([(ev[1][i] - ev[0][i]) for i in range(np)])
    ep = list([(ev[1][i] + ev[0][i])*0.5 for i in range(np)])

    x = zeros(ne)
    y = zeros(ne)
    for i in range(4):
        for j in range(4):
            p = 1 - (i+j)*2/3.
            vap = vav**p * (scx**p - scm**p)
            x += em[i] * em[j] * vap / p
            y += ep[i] * ep[j] * vap / p

    Da = 1000 * (x/((scx - scm) * vav))**0.5
    Dr = 100 * (x/y)**(1/2.)

    D1 = Da / (vav/30 * (b[1] + b[0])/(fac*200))

    #return Da, Dr, D1

    s = '''--------------------
# delta between %s and %s by dpt.py
# dlt [meV/atom]  rltv dlt [%%]  dlt1 [meV/atom]
--------------------
''' % tuple(fls)

    # crude sorting for print
    cdx = list([ce.index(el) for el in ptl if el in ce])

    s += '\n'.join(['%s %10.5f %10.5f %8.3f'%(ce[i].title().ljust(2), Da[i], Dr[i], D1[i])  for i in cdx])
    s += '\n--------------------'

    daix = Da.argmax(); drix = Dr.argmax(); d1ix = D1.argmax()
    daim = Da.argmin(); drim = Dr.argmin(); d1im = D1.argmin()

    s += '\nnp.mean  %8.3f %5.1f %8.3f' % (Da.mean(), Dr.mean(), D1.mean())
    s += '\nnp.std   %8.3f %5.1f %8.3f' % (Da.std(), Dr.std(), D1.std())
    s += '\nnp.max   %8.3f %5.1f %8.3f (%s, %s, %s)' % (Da[daix], Dr[drix], D1[d1ix], ce[daix], ce[drix], ce[d1ix])
    s += '\nnp.min   %8.3f %5.1f %8.3f (%s, %s, %s)' % (Da[daim], Dr[drim], D1[d1im], ce[daim], ce[drim], ce[d1im])

    s += '\n--------------------\n'

    return s

if __name__ == '__main__':
    import argparse
    from multiprocessing import Pool
    from subprocess import Popen, PIPE

    aprs = argparse.ArgumentParser(description='deltacores periodic table (dpt) plotter')
    aprs.add_argument('-j', type=int, help='limit number of processes spawned (default is to use all visible cpus)')
    aprs.add_argument('-p', '--pattern', default='dc/{e}-{s}.d/scflog.{s}.{e}.log', help='pattern to use to find the log files, for example: "dc/{e}-{s}.d/scflog.{s}.{e}.log"')
    aprs.add_argument('-e', '--eosfile', default='eos-lmf.txt', help='where to summarise eq of state data')
    aprs.add_argument('-r', '--reference', default='eos-ref.txt', help='reference to compare with')
    aprs.add_argument('-d', '--delta', default='eos-delta-lmf-ref.txt', help='where to write the delta between the eosfile and the ref file')


    args = aprs.parse_args()

    nproc = args.j
    elptr = args.pattern
    eosfile = args.eosfile
    reffile = args.reference
    dltfile = args.delta
    dptfile = re.sub(r'\.txt$','', dltfile) + '.svg'

    assert(eosfile != reffile)
    assert(eosfile != dltfile)
    assert(reffile != dltfile)

    elscl_loc = elscl_t(elptr)

    p = Pool(nproc)
    sys.stdout.write('processing... ')
    sys.stdout.flush()
    eoss = p.map(elscl_loc, ptl)
    p.close()
    del p
    sys.stdout.write('\n')
    sys.stdout.flush()
    sys.stderr.flush()

# the same as above but strctly serial
    #eoss = list(map(elscl_loc, ptl))

    s = '\n'.join([eos for eos in eoss if eos != ''])
    writefile(eosfile, s)
    print(eosfile + ' written')

    if (os.path.exists(reffile)):
        o = dcdiff((eosfile, reffile))
        writefile(dltfile, o)
        print(dltfile + ' written')

        de = parse_dlt(o)
        s = ptdraw(de)
        writefile(dptfile, s)
        print(dptfile + ' written')
    else:
        print(reffile + ' not found, skipping delta and table')

