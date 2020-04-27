#!/usr/bin/env python

# dmt 2014-2017

import os,sys,re,time
from subprocess import Popen, PIPE
from functools import reduce


#if sys.version_info.major > 2: raw_input = input
if sys.version_info[0] > 2:
    raw_input = input
else:
    from codecs import open

def loadfile(fn):
    f = open(fn,mode='r',encoding='utf8')
    s = f.read()
    f.close()
    return s

def writefile(fn,s,mode='w'):
    f = open(fn,mode=mode,encoding='utf8')
    f.write(s)
    f.close()

def getmakeflags(ff, prepend = ''):
# useful makefile escapes:
#   qhdr = module use modpath\\nmodule load special/12 \$$var
    makeprepend = 'vars0 := $(.VARIABLES)\n' + prepend + '\n'
    makeappend = '''
nvars := $(filter-out $(vars0) vars0, $(.VARIABLES))
sflgs := $(filter lessflags%, $(nvars))
sfils := $(filter lessfiles%, $(nvars))

echoflags:
\t@echo cc      : $(cc)
\t@echo cxx     : $(cxx)
\t@echo cflags  : $(cflags)
\t@echo omp     : $(omp)
\t@echo nvcc    : $(nvcc)
\t@echo ccomp   : $(ccomp)
\t@echo fc      : $(fc)
\t@echo ffxc    : $(ffxc)
\t@echo fflags  : $(fflags)
\t@echo ppflags : $(ppflags)
\t@echo soft    : $(soft)
\t@echo ar      : $(ar)
\t@echo ldflags : $(ldflags)
\t@echo modflag : $(modflag)
\t@echo knrcc   : $(knrcc)
\t@echo gwflags : $(gwflags)
\t@echo prefix  : $(prefix)
\t@echo bldpath : $(bldpath)
\t@echo mpirun  : $(mpirun)
\t@echo qcmd    : $(qcmd)
\t@echo qhdr    : $(qhdr)

echospecials:
\t$(foreach v, $(sflgs), $(info $(v) : $($(v))))
\t$(foreach v, $(sfils), $(info $(v) : $($(v))))

'''

    tmp = '.flags.mk'
    sff = loadfile(ff)
    writefile(tmp,makeprepend + sff + makeappend)

    p = Popen(['make','--no-print-directory','-r','-R','-f',tmp,'echoflags'],stdout=PIPE,shell=False,env=os.environ,cwd=None)
    o,e = p.communicate()
    del p

    o = o.strip().splitlines()
    flags = {}
    for l in o:
        k,v = tuple(l.decode('ascii').split(':',1))
        flags[k.strip()] = v.strip()

    p = Popen(['make','--no-print-directory','-r','-R','-f',tmp,'echospecials'],stdout=PIPE,shell=False,env=os.environ,cwd=None)
    o,e = p.communicate()
    del p

    o = o.strip().splitlines()
    sflags = {}
    for l in o:
        if l.startswith(b'lessflags'):
            k,v = tuple(l.decode('ascii').split(':',1))
            sflags[k.strip()[len('lessflags'):]] = v.strip()

    special = {}
    for l in o:
        if l.startswith(b'lessfiles'):
            k,v = tuple(l.decode('ascii').split(':',1))
            for f in v.split():
                special[f] = sflags[k.strip()[len('lessfiles'):]]

    flags['special'] = special

    #os.remove(tmp)

    return flags


def getflags(srcpath):

    flags0 = dict(
        cc      = '',
        cxx     = '',
        cflags  = '',
        omp     = '',
        nvcc    = '',
        ccomp   = '',
        fc      = '',
        ffxc    = '',
        fflags  = '',
        ppflags = '',
        soft    = '',
        ar      = 'ar',
        ldflags = '',
        modflag = '',
        knrcc   = '',
        gwflags = '',
        prefix  = '/opt/lm',
        bldpath = 'bld',
        mpirun  = '',
        qcmd    = 'tm -n NPROC -N NAME -o OUTFILE --',
        qhdr    = '',
        special = {}
    )
    #flf = 'flags.py'
    #if os.path.isfile(flf):
        #from flags import flags
        ##print flags
        #for k in flags0:
            #if not flags.has_key(k):
                #flags[k] = flags0[k]
        #return flags

    del flags0['special']
    sflags0 = reduce(lambda s, p: s + p[0] + ' = ' + str(p[1]) + '\n', flags0.items(), '')

    flf = 'flags.mk'
    if not os.path.isfile(flf):
        sys.stdout.flush(); sys.stderr.flush()
        sys.stderr.write('\nflags.mk not found. Enter args for genflags.py below:\n')
        Popen([srcpath + '/genflags.py','shortuse'],shell=False,stdout=sys.stderr).wait()
        sys.stderr.write(srcpath.rstrip('/') + '/genflags.py ')
        flargs = raw_input()
        fflf = open(flf,'w')
        Popen([srcpath + '/genflags.py']+flargs.split(),shell=False,stdout=fflf).wait()
        fflf.close()
        sys.stderr.write('Flags saved in flags.mk for reuse.\n')
    return getmakeflags(flf, prepend = sflags0)


def fndfl(fl,flds):
    #return the index of the path in which fl is found
    for i in range(len(flds)):
        if os.path.isfile(flds[i] + '/' + fl):
            return i
    return -1

#def fname(p):
    '''os.path.basename'''
    #return p[p.rfind('/')+1:]

#def addprefix(p,o):
    #return ' '.join([(p+i) for i in o.split()])
#def addsuffix(s,o):
    #return ' '.join([(i+s) for i in o.split()])

class animation_t(object):
    def __init__(self, skip=1):
        self.c = 0
        self.s = skip
        self.f = sys.stdout
        self.a = r'\|/-'
        self.n = len(self.a)
        self.d = '.'

        self.f.write(self.d)

    def step(self):
        if self.c % self.s == 0 and self.c > 0:
            d = self.d if (self.c//self.s)%self.n == 1 else ''
            self.f.write('\b'+d+r'\|/-'[(self.c//self.s)%self.n])
            self.f.flush()
        self.c += 1

    def __del__(self):
        self.f.write('\b.\n')
        self.f.flush()


def bld(s, r, o, si = [], oi = [], so = [], vs = {}):
    '''
        s: sources
        r: rule
        o: explicit outputs
        si: implicit sources
        oi: implicit outputs
        so: order only sources
        vs: environment flags to override
    '''

    #s = s.split()
    #o = o.split()
    #si = si.split()
    #oi = oi.split()
    #so = so.split()
    #try:

    outs = ' '.join(o)

    if oi != []:
        outs += ' | ' + ' '.join(oi)

    t = 'build {outs}: {rule} {ins}'.format(outs = outs, rule = r, ins = ' '.join(s))
    #except:
        #print >>sys.stderr, o, oi,s,si
        #sys.exit(-1)
    if si != []: # implicit sources (inputs/sources not showing up in $in)
        t += ' | ' + ' '.join(si)

    if so != []: # order only sources
        t += ' || ' + ' '.join(so)

    for k,v in vs.items():
        t += '\n    {k} = {v}'.format(k = k, v = v)

    return [t]


def findall_to_lower(p, d):
    return list(i.lower() for i in p.findall(d))

#from subprocess import Popen,PIPE

#def cgrep(p,f):
    #t = Popen(['grep','-HoPe',p,f],shell=False,stdout=PIPE)
    #return t.communicate()[0]

def laddprefix(p, o):
    if hasattr(o, 'split'): o = o.split()
    return list(map(lambda i: p+i, o))

def laddsuffix(s, o):
    if hasattr(o, 'split'): o = o.split()
    return list(map(lambda i: i+s, o))

def lsandwich(p, o, s):
    if hasattr(o, 'split'): o = o.split()
    return list(map(lambda i: p+i+s, o))

def sandwich(p, o, s):
    return ' '.join(lsandwich(p, o, s))

def frontreplace(p, q, o):
    if hasattr(o, 'split'): o = o.split()
    n = len(p)
    return list(map(lambda i: (q+i[n:]) if i.startswith(p) else i, o))

#def backreplace(p, q, o):
    #if hasattr(o, 'split'): o = o.split()
    #n = len(p)
    #return list(map(lambda i: (i[:-n]+q) if i.endswith(p) else i, o))

def uniq(l,d=1):
    #shoddy nonscrambling uniq:
    #   d= 1 preserve first occurence
    #   d=-1 preserve last occurence
    ul = []
    for i in l[::d]:
        if i not in ul:
            ul.append(i)
    ul = ul[::d]
    return ul


class nj_t(object):

    msi_ptr = re.compile(r'(?mix)^\s*use\s+(\w+)')
    isi_ptr = re.compile(r'''(?mix)^\s*include\s+['"]?([\w\-._]+)['"]?''')
    oi_ptr = re.compile(r'(?mix)^\s*module\s+(\w+)\s*(?:$|!)')

    def __init__(self, modpath = 'mods', bldpath = 'bld', libpath = 'lib', libsfx='.a', rules = None, fppstyle='ccomp'):
        srcpath = os.path.dirname(sys.argv[0])
        flags = getflags(srcpath)
        #sys.stderr.write(str(flags)+'\n')
        if flags['ffxc'] == '': flags['ffxc'] = flags['fc']

        if flags['bldpath'] != '': bldpath = flags['bldpath']
        #if os.path.realpath(srcpath) != os.path.realpath('.'): bldpath = '.'
        if not os.path.abspath(srcpath).startswith(os.path.abspath('.')): bldpath = '.'

        self.srcpath = srcpath
        self.bldpath = bldpath
        self.modpath = modpath
        self.libpath = libpath
        self.libsfx = libsfx
        self.flags = flags
        self.includes = {}
        self.emods = set(['iso_c_binding', 'iso_fortran_env', 'omp_lib'])
        self.eincs = set([]) # external includes
        self.gincs = set([]) # generated includes
        self.fppstyle = fppstyle

        self.mkclk = 0
        self.dclk = 0
        self.lclk = 0
        self.oiclk = 0
        self.isiclk = 0
        self.msiclk = 0

        self.direct = ''

        self.gen_rules(rules = rules)

        self.animation = animation_t(skip=32)

    def gen_rules(self, rules = None):
        self.rules = '''

# Export this env to change progress reporting
# NINJA_STATUS="[%c %r %s/%t %p] "

builddir = $bldpath

ccomp = $bldpath/ccomp

implicit_cflags    = -I$bldpath
implicit_cxxflags  = -I$bldpath
implicit_nvccflags = -I$bldpath
implicit_fflags    = -I$bldpath
implicit_ldflags   = -L$bldpath/$libpath -Wl,-rpath,'$$ORIGIN/$libpath:$$ORIGIN/../$libpath'
mods = $modflag $bldpath/$modpath -I$bldpath/$modpath

rule cp
    command = cp $in $out
rule install
    command = install -D $in $out

rule cxxc
    command = $cxx $implicit_cxxflags $cxxflags -c $in -o $out
rule cxxl
    command = $cxx $implicit_cxxflags $cxxflags $in $implicit_ldflags $libs $ldflags -o $out

rule nvcc
    command = $nvcc $implicit_nvccflags $nvccflags -c $in -o $out
rule nvcl
    command = $nvcc $implicit_nvccflags $nvccflags $in $implicit_ldflags $libs $ldflags -o $out

rule cc
    command = $cc $implicit_cflags $cflags -c $in -o $out
rule cl
    command = $cc $implicit_cflags $cflags $in $implicit_ldflags $libs $ldflags -o $out

rule fc
    command = $fc $implicit_fflags $mods $fflags -c $in -o $out
rule fxc
    command = $ffxc $implicit_fflags $mods $fflags -c $in -o $out
rule ffp
    command = $ccomp -c! $ppflags $in $out
rule fxp
    command = $ccomp $ppflags $in $out
rule fl
    command = $fc $implicit_fflags $fflags $in $implicit_ldflags $libs $ldflags -o $out

rule ar
    command = rm -f $out && $ar rc $out $in

rule cmd
    command = echo command not given

rule mkdir
    command = mkdir -p $out

rule nj
    command = python $in $out
    generator = 1
    pool = console

'''
        if rules != None: self.rules = rules
        for k,v in self.flags.items():
            if k != 'special' and k != 'bldpath':
                self.rules = k+' = '+v + '\n' + self.rules

        self.rules = 'libpath = ' + self.libpath + '\n' + self.rules
        self.rules = 'modpath = ' + self.modpath + '\n' + self.rules
        self.rules = 'bldpath = ' + self.bldpath + '\n' + self.rules
        self.rules = 'srcpath = ' + self.srcpath + '\n' + self.rules
        self.rules = 'ninja_required_version = 1.7\n' + self.rules

        self.js = []

    def grepoi(self,d):
        clk1 = time.clock()
        r = findall_to_lower(nj_t.oi_ptr, d)
        #r = re.findall(r'(?m)^\s*module\s+(\w+)\s*$',d)
        clk2 = time.clock()
        self.oiclk += clk2-clk1
        return r

    def grepmsi(self,d):
        clk1 = time.clock()
        r = findall_to_lower(nj_t.msi_ptr, d)
        clk2 = time.clock()
        self.msiclk += clk2-clk1
        return r

    def grepisi(self,d):
        clk1 = time.clock()
        r = nj_t.isi_ptr.findall(d)
        clk2 = time.clock()
        self.isiclk += clk2-clk1
        return r

    def bld(self, *args, **kwds):
        self.js.extend(bld(*args, **kwds))
        self.animation.step()

    def resolve_fincls(self, s,r,vs={}):
        if s in self.includes: return self.includes[s]

        clk1 = time.clock()
        d = loadfile(self.srcpath+'/'+s)
        clk2 = time.clock()
        self.lclk += clk2-clk1

        ois = set(self.grepoi(d)) - set(['procedure'])
        msi = set(self.grepmsi(d)) - self.emods
        oi = lsandwich('$bldpath/$modpath/', ois, '.mod')
        si = lsandwich('$bldpath/$modpath/', msi - ois, '.mod')
        #si += resolve_incls(s,r,vs=vs)

        #print s, set(msi_ptr.findall(d)) - ois - set(['iso_c_binding','procedure'])

        sis = list(set(self.grepisi(d)) - self.eincs)

        fsis = []
        if sis != []:
            #print 'sis',sis
            f = r+'flags'
            #print flags[f]
            incs = re.findall(r'\-I\s*([$\w._][$\w/._\-]*)', vs[f] if f in vs else self.flags[f]) + [os.path.dirname(s),'$bldpath']
            #if '$bldpath/$modpath' in incs: incs.remove('$bldpath/$modpath')
            incs = frontreplace('$bldpath','',incs)
            incs = frontreplace('/','',incs)
            incs = laddprefix(self.srcpath+'/',incs)
            #incs = map(lambda i: os.path.normpath(i),incs)
            #print incs

            #incs = frontreplace(bldpath,'',incs)
            #print 'sled', incs
            #print incs, sis
            for sisi in sis:
                fnd = fndfl(sisi,incs) # da se napravi da tyrsi iz receptite syshto.. za generirani fajlove
                #print incs[fnd]
                if fnd != -1:
                    ifl = os.path.normpath(os.path.join(frontreplace(self.srcpath+'/','',[incs[fnd]])[0],sisi))
                    if ifl not in self.includes: self.includes[ifl] = self.resolve_fincls(ifl,r,vs=vs)
                    fsis.append('$bldpath/'+ifl)
                else:
                    #if not found then we assume it will be generated in the build and blindly put a dependency on it. If it is not generated, the build will fail with clear enough message from the compiler.
                    fsis.append('$bldpath/'+sisi)
        si.extend(fsis)

        return si,oi

                    #if fnd.startswith(srcpath):
                        #fnd = '$srcpath'+fnd[len(srcpath):]

    def mk(self,s,r=None,o=None,si=None,oi=None,so=None,vs={}):
        '''
            s: sources
            r: rule
            o: explicit outputs
            si: implicit sources
            oi: implicit outputs
            so: order only sources
            vs: environment flags to override

        '''

        clk1 = time.clock()

        if hasattr(s,'split'): s = s.split()
        if o == None: assert(len(s) == 1)

        #modes:
        #   s:list, o:single(mandatory to pass)
        #   s:single, o:single(optional)

        #print s,o,si,oi,so

        if si == None: si = []
        if oi == None: oi = []
        if so == None: so = []
        lvs = vs.copy()
        if o == None: o = s[0]+'.o' # o not being passed (optional) implies s is single file

        oe = o.rsplit('.',1)[1] if o.rfind('/') < o.rfind('.') else ''

        #assert((oe != '' and oe != 'so') or r != None) # if oe=='' probably linking execuable, oe=='so' dyn lib, in bth cases explicit rule is needed or default may be used. .

        t = []
        b,se = '',''
        if oe not in ('a','so',''):
            assert(len(s) == 1)
            s = s[0]
            sb,se = s.rsplit('.',1)
            b = sb
            if (o != s+'.o'):
                b = o.rsplit('.',1)[0]
                if b.endswith('.'+se):
                    b = b.rsplit('.',1)[0]


            if r == None:
                if se[0].lower() == 'f':
                    r = 'f'
                elif se == 'c':
                    r = 'c'
                elif se in ('cc','cpp','cxx'):
                    r = 'cxx'
                elif se == 'cu':
                    r = 'nvc'
            if si == [] or oi == []:
                clk10 = time.clock()
                fsi,foi = self.resolve_fincls(s,r,vs=vs)
                clk20 = time.clock()
                self.dclk += clk20-clk10

                si.extend(fsi)
                oi.extend(foi)

            if list(filter(lambda i: i.endswith('.mod'), si+oi)) == []: lvs['mods'] = ''

            ff = se[1:] in ('90','95','03','08') # is free form

            if r == 'f' and not ff: r += 'x'
            r += 'c'

            ppflags = lvs['ppflags'] if 'ppflags' in lvs else ''
            if 'ppflags' in lvs: del lvs['ppflags']

            if self.fppstyle == 'ccomp' and (se == 'f' or se == 'for' or se[0] == 'F') and (ppflags.strip() != ''):
                pt = 'ffp' if ff else 'fxp'
                m = '$bldpath/'+b+'_PP.'+se.lower()
                self.bld(['$srcpath/'+s], pt, [m], so = ['$ccomp'], vs=dict(ppflags=ppflags))
                s = m
            else:
                s = '$srcpath/'+s
            s = [s]
        else:
            if r == None:
                if oe == 'a':
                    r = 'ar'
                elif oe in ('','so'):
                    r = 'fl'
                    ls = [i[2:] for i in s if i.startswith('-l')]
                    si.extend(lsandwich('$bldpath/'+self.libpath+'/lib',ls,self.libsfx))
                    lvs['libs'] = ' '.join(laddprefix('-l',ls))
                    if lvs['libs'] == '':
                        lvs['implicit_ldflags'] = ''
                        lvs.pop('libs')
                    s = [i for i in s if not i.startswith('-l')]
                    if oe == 'so':
                        lvs['ldflags'] = '-shared'

            s = laddprefix('$bldpath/',s)

        if oi != []: lvs['restat'] = 1
        if o in self.flags['special']: lvs['fflags'] = self.flags['special'][o]
        self.bld(s, r, ['$bldpath/'+o], si=si, oi=oi, so=so, vs=lvs)

        clk2 = time.clock()
        self.mkclk += clk2 - clk1

        #if (s[0].endswith('mpi_mod.F90')): print (si)

    def mks(self,s,vs={}):
        '''
            s: sources
            vs: environment flags to override

            rules, explicit, implicit and order only sources and outputs will hopefully be figured out automatically here.
        '''
        clk1 = time.clock()
        if hasattr(s, 'split'): s = s.split()
        clk2 = time.clock()
        self.mkclk += clk2 - clk1

        for i in s:
            #t.extend(self.mk(i,vs=vs))
            self.mk(i,vs=vs)

    def mkl(self,s,l,r=None,si=None,oi=None,so=None,vs={}):
        self.mk(s,r=r,o=self.libpath+'/lib'+l+self.libsfx,si=si,oi=oi,so=so,vs=vs)

    def mkincludes(self):
        '''issue copy targets for includes and other implicit sources dependencies'''
        for k,v in self.includes.items():
            if (k not in self.gincs):
                self.bld([os.path.join('$srcpath',k)], 'cp', [os.path.join('$bldpath',k)], si = v[0], oi = v[1])

    def printtimes(self):
        sys.stderr.write('mks: '   +str(self.mkclk)+'s\n')
        sys.stderr.write('fincls: '+str(self.dclk)+'s\n')
        sys.stderr.write('loads: ' +str(self.lclk)+'s\n')
        sys.stderr.write('goi: '   +str(self.oiclk)+'s\n')
        sys.stderr.write('gmsi: '  +str(self.msiclk)+'s\n')
        sys.stderr.write('gisi: '  +str(self.isiclk)+'s\n')
        sys.stderr.write('goi + gmsi + gisi: '+str(self.oiclk + self.msiclk + self.isiclk)+'s\n')


    def write(self, fln='build.ninja', lw=0):
        clk1 = time.clock()
        o = '\n'.join(uniq(self.js)) + '\n' + self.direct + '\n'
        clk2 = time.clock()
        if '-time' in sys.argv: sys.stderr.write('uniq: '+str(clk2-clk1)+'s\n')


        if lw > 0:
            # reset lw or if a word is longer than lw will bring a disaster
            # in principle it is possible to break midword but this is less hasle for now. if done care must be taken to to avoid fake $ escapes ($$).
            lw = max(lw,1+max(len(i) for i in o.split()))

            ow = []
            for l in o.split('\n'):
                l = l.rstrip()+' '
                i = 0
                n = len(l)

                while True:
                    j = 1+i+l[i:i+lw].rfind(' ')
                    e = '$' if j < n else ''
                    ow.append((l[i:j]+e).rstrip())
                    if j >= n: break
                    i = j

                ## simple but ugly and incorrect when it hits on $ (makes fake $$ escapes)
                #n = len(l)//lw
                #for i in range(n):
                    #ow.append(l[i*lw:(i+1)*lw]+'$')
                #ow.append(l[n*lw:])
            o = '\n'.join(ow)

        writefile(fln,o,mode='a')

        #print(includes)

