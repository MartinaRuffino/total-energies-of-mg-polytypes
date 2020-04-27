#!/usr/bin/env python

# shall run with both py2 and py3

import os, sys, re
import subprocess

def xc(cmd):
    return subprocess.check_output(cmd.split())

def readfile(fln):
    f = open(fln, 'r')
    s = f.read()
    f.close()
    return s

class CtlEntry:
    def __init__(self,name,value):
        self.name = name
        self.value = value

    def prtkv(self,f):
        if self.value == "":
            buf = " " + self.name 
        else:
            buf = " " + self.name + "="

        if type(self.value) is tuple:
            for i in self.value:
                buf = buf + " " + str(i)  # this adds an extra preceeding space for tuple values :-/
        else:
            buf = buf + str(self.value)
        f.write(buf)

class ctrl:
    def __init__(self):
        self.sec = []
        self.sub = []

    def add(self,secname,name,value):
        secpresent = 0
        for i in self.sec:
            if i[0] == secname:
                secpresent = 1
                i[1].append(CtlEntry(name,value))
        if secpresent == 0:
            self.sec.append((secname,[CtlEntry(name,value)]))

    def addsub(self,secname,subsec,name,value):
        secpresent = 0
        subpresent = 0

        for i in self.sub:
            if i[0] == secname:
                secpresent = 1

                for j in i[1]:
                    if j[0] == subsec:
                        subpresent = 1
                        j[1].append(CtlEntry(name,value))

                if subpresent == 0:
                    i[1].append((subsec,[CtlEntry(name,value)]))

        if secpresent == 0:
            self.sub.append((secname,[(subsec,[CtlEntry(name,value)])]))

    def prt(self,f):
        for i in self.sec:
            if i[0] == "% const":
                # special case for const declarations
                f.write("\n" + i[0] ) ## includes extra space between blocks
            else:
                f.write("\n" + i[0] + "\n") ## includes extra space between blocks
            ctr = 1
            for j in i[1]:
                j.prtkv(f)
            for j in self.sub:
                if j[0] == i[0]:  # any subsections in this part
                    for k in j[1]:
                        f.write("\n " + k[0] + "[")
                        for l in k[1]:
                            l.prtkv(f)
                        f.write("]")
                    
        for j in self.sub:  # if a subsection exists on its own then it would be missed
            printed = 0
            for i in self.sec:
                if j[0] == i[0]:  # any subsections in this part
                    printed = 1
            if printed == 0:
                f.write("\n" + j[0]) ## includes extra space between blocks
                for k in j[1]:
                    f.write("\n " + k[0] + "[")
                    for l in k[1]:
                        l.prtkv(f)
                    f.write("]")
        f.write("\n") # terminating newline
 
def genctl(opts,kpts):

    a = ctrl() # class instantiation

    pt = '''
     H                                                                                            He
     Li Be                                                                         B  C  N  O  F  Ne
     Na Mg                                                                         Al Si P  S  Cl Ar
     K  Ca                                           Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
     Rb Sr                                           Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
     Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
    '''

    keyval = []
    zctr = 1
    for i in pt.split():
        keyval.append((i,zctr))
        zctr = zctr + 1
    atnum = dict(keyval)

    capelem = opts['elem'].capitalize()

    zval = atnum[capelem]

    # defaults for lmax/lmaxa
#    if ( zval > 54 ):
#        lmx = 4
#    elif ( zval > 18 ):
    if ( zval > 18 ):
        lmx = 3
    elif ( zval > 2 ):
        lmx = 2
    else:
        lmx = 1
    
    # default for lmxa also
    lmxa = lmx + 1

    rmax = 10.0
    if float(opts['rmt']) > rmax:
        rmt = rmax
    else:
        rmt = opts['rmt']

    a.add("% const","sc",'1.0')
    a.add("% const","minsc",'0.94') 

    a.add("VERS","LM",7)
    a.add("VERS","FP",7)
    a.add("IO","VERBOS",30)

    a.addsub("HAM","AUTOBAS","LMTO",5)
    a.addsub("HAM","AUTOBAS","MTO",14)
    a.addsub("HAM","AUTOBAS","LOC",1)
    a.addsub("HAM","AUTOBAS","PFLOAT",11)
    a.addsub("HAM","AUTOBAS","PNU",opts['pnu'])
    a.addsub("HAM","AUTOBAS","QLOC",0.002)
    a.addsub("HAM","AUTOBAS","ELOC",-2.2)
    a.addsub("HAM","AUTOBAS","EH",-0.5)

    if 'pack' in opts and opts['pack'][0] < 0.3:
        a.add("HAM","PWMODE",11)
        pwemax = 2.0
        if opts['pack'][0] < 0.12 or opts['pack'][0] < opts['pack'][1]:
            pwemax = 4.0
        a.add("HAM","PWEMAX", pwemax)
       

    #a.addsub("HAM","AUTOBAS","EIN",2.0)
    #a.addsub("HAM","AUTOBAS","EOUT",2.0)

    # for antiferromagnetic cases, print out the second mom also
    af = 1
    nsp = 1
    if 'mmom' in opts:
        comp = len(opts['mmom'])
        #if comp == 0:
        #    nsp = 1
        #    af = 1
        if comp == 1:
            nsp = 2
        if comp == 2:
            nsp = 2
            af = 2
        if comp > 2:
            print("logic error")
       
    a.add("HAM","GMAXDA",opts['gmax']) # change to GMAXDA to better handle cell scaling
    a.add("HAM","XCFUN",(0,101,130))
    a.add("HAM","NSPIN",nsp)
    a.add("HAM","ELIND",-1.0)
    a.add("HAM","TOL",1e-16)

    if opts['elem'] not in kpts:
        print('kpts missing for '+ capelem)
        exit()

    a.add("BZ","NKABC",kpts[opts['elem']])
    a.add("BZ","METAL",5)

    a.add("EWALD","TOL",1e-16)

    a.add("ITER","NIT",500)
    a.add("ITER","CONV",1e-5)
    a.add("ITER","CONVC",1e-4)
    a.add("ITER","MIX","B6,b=0.9,n=6;B6,b=0.2") # beware, MIX doesn't work like a proper subcategory (and spaces are fatal)

    a.add("SYMGRP","find","")

    a.add("SITE","FILE","site")
    a.add("STRUC","FILE","site")

    if 'alat' in opts:
        a.add("STRUC","DALAT",'{'+str(alat)+'*(sc^(1/3)-1.0)}')

    for isp in range(0,af):
        if isp == 0:
            a.add("SPEC","ATOM",capelem)
        else:
            a.add("SPEC","\n","") # newline for clarity
            a.add("SPEC","ATOM",capelem+"X")
        a.add("SPEC","R",'{'+str(rmt)+'*minsc^(1/3)}') # touching spheres at smallest volume
        a.add("SPEC","Z",zval)
        a.add("SPEC","LMX",lmx)
        a.add("SPEC","LMXA",lmxa)
        a.add("SPEC","KMXA",5)
        a.add("SPEC","KMXV",5)
        a.add("SPEC","A",0.01)

        if 'pz' in opts:
            a.add("SPEC","\n PZ",(opts['pz'][0],opts['pz'][1],opts['pz'][2],opts['pz'][3]))
    
        if nsp == 2:
            a.add("SPEC","MMOM",(opts['mmom'][isp][0],opts['mmom'][isp][1],opts['mmom'][isp][2],opts['mmom'][isp][3],))

    f = open('ctrl.' + opts['elem'], 'w')  # setup basic control file
    a.prt(f)
    f.close()

if __name__ == '__main__':
    import argparse

    aprs = argparse.ArgumentParser(description='init2ctrl generates ctrl.ext from init.ext for lmf')
    aprs.add_argument('-l', '--launcher', default='', help='special process launcher, for example "mpirun -n 1" if your mpi setup is broken and does not work in singleton mode')
    aprs.add_argument('-k', '--kpts', default='kpts.lmto.summary', help='kpoints list file')
    aprs.add_argument('element', help='will use init.element file to bootstrap lmf calculation')

    args = aprs.parse_args()

    extn = args.element
    launcher = args.launcher
    kpts_file = args.kpts

    init_file = 'init.' + extn

    if not os.path.isfile(init_file):
        print("%s not found" % init_file)
        exit()

    kpts = dict([l.split(None,1) for l in readfile(kpts_file).splitlines()])

    initfile = readfile(init_file)
    mmom_m = re.findall(r'(?xs) MMOM \s* = \s* ( \-?\d*\.?\d* (?: [\s,]+ \-?\d*\.?\d*)*)', initfile)
    mmom = []
    for l in mmom_m:
        tmom = list(map(lambda s: int(float(s)), l.replace(',',' ').split()))
        tmom.extend([0]*max(0,(4-len(tmom))))
        mmom.append(tmom)


    readin = xc(launcher + " blm init." + extn).decode() # (this could also be used for rmt)

# obtain new muffin tin radius using blm
# same expression works for lmchk output, too
# readin = subprocess.check_output(["lmchk","--getwsr",extn])
# spec  name        old rmax    new rmax     ratio
#   1   K           4.324832    4.918254    1.137213

    rmt = re.search(r'(?xs) old \s+ rmax \s+  new \s+ rmax \s+ ratio \n \s*? 1 \s* \w+ \s* (\d*\.?\d*)', readin).group(1)

# print out the packing fraction
#Initial sphere packing = 74%
    pack = list(map(lambda s: float(s)/100., re.findall(r'(?xs) Initial \s+ sphere \s+ packing \s* = \s* (\d*\.?\d*) \s* \% \s* scaled \s+ to \s* (\d*\.?\d*) \s* \%', readin)[0]))

# generate a trial ctrl then run lmfa
    opts = dict(gmax=1, elem=extn, rmt=rmt, pnu=1, mmom=mmom)
    genctl(opts,kpts)

# get new gmax estimate
    readin = xc(launcher + " lmfa --usebasp ctrl." + extn).decode()
    #GMAX=8.5 (valence)  8.5 (local orbitals)
    pvals_m = re.search(r'(?xs) GMAX \s* = \s* (\d*\.?\d*) (?: \s+ \(valence\) \s+ (\d*\.?\d*) \s+ \(local\s+orbitals\))?', readin)
    pval = float(pvals_m.group(1))
    ploc = pvals_m.group(2)
    pval = max(pval, 6.0)
    if ploc != None: pval = max(pval, float(ploc))

# get PNU for setting high-lying PZ, requires that PNU=1 in first setup
    pz = [0,0,0,0]

# in this case, the Pnu string from lmfa always contains spdf,g values
    parray = list(map(float, re.search(r'(?xs) Autogenerated \s+ Pnu\s*: \s* ((\d*\.?\d*\s+)+)', readin).group(1).split()))
    parray.extend([0]*max(0, 4-len(parray)))

# check if an l=2 p value exists; p>0.5 (ish) for occupied d cases
    ipzd = int(parray[2])
    fpzd = parray[2] - float(ipzd)
    if fpzd > 0.51: pz[2] = ipzd+1.3

# for Lu--W, add a 5f HLLO, too
    ipzf = int(parray[3])
    fpzf = parray[3] - float(ipzf)
    if fpzf > 0.51: pz[3] = ipzf+1.3


# for DALAT machinery...
#% site-data vn=3.0 fast io=15 nbas=3 alat=14.43758894 plat= 1.0 0.0 0.0 0.9189977 0.3942629 0.0 0.9189977 0.1888104 0.3461125
    sitedata = readfile('site.' + extn)
    alat = re.search(r'(?xs) alat \s* = \s* (\d*\.?\d*)', sitedata).group(1)

    opts = dict(gmax=pval, elem=extn, rmt=rmt, alat=alat, pz=pz, pnu=0, mmom=mmom, pack=pack)
    genctl(opts,kpts)

# clean up junk from blm
    subprocess.check_output('rm -f actrl.{x} log.{x} rst.{x} moms.{x} mixm.{x}'.format(x=extn).split())
