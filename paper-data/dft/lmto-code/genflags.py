#!/usr/bin/env python

from __future__ import print_function
import os, sys, re
from glob import glob as glob
from subprocess import Popen,PIPE


def getenv(n,d=''):
    return os.environ[n] if n in os.environ.keys() else d

def xc(c):
    if hasattr(c,'split'): c = c.split()
    p = Popen(c, stdout=PIPE, stderr=PIPE, shell=False)
    r = p.wait()
    s = p.communicate()
    del p
    return s + (r,)

def uniqf(l, f=lambda i: False):
    ul = []
    for i in l:
        if not ((i in ul) and f(i)):
            ul.append(i)
    return ul

plts = dict(
    gcc = dict(
        cc = 'gcc',
        cxx = 'g++',
        fc = 'gfortran',
        omp = '-fopenmp',
        cflags = dict(
            dbg = '-g -O0 -fstack-protector-all -Wstack-protector -Wall -pipe',
            opt = '-O3 -march=native -pipe -g1',
            opg = '-O3 -march=native -pipe -g -fstack-protector-all -Wstack-protector -Wall',
        ),
        fflags = dict(
             #-Wno-argument-mismatch can be used to quieten things a bit on v7+
            dbg = '-g -O0 -pipe -fbacktrace -ffpe-trap=invalid,zero,overflow -std=legacy -Wno-argument-mismatch -finit-integer=2147483647 -finit-real=snan -finit-character=35 -finit-logical=true -fstack-protector-all',
            # -fmax-stack-var-size=512 -fcheck=all
            opt = '-O3 -march=native -pipe -g1 -std=legacy -Wno-argument-mismatch',
            opg = '-g -O3 -march=native -pipe -fbacktrace -ffpe-trap=invalid,zero,overflow -std=legacy -Wno-argument-mismatch -finit-integer=2147483647 -finit-real=snan -finit-character=35 -finit-logical=true -fstack-protector-all',
        ),
        modflag = '-J',
        ldflags = '-L${MKLROOT}/lib/intel64 -lfftw3xf_gnu -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread ',
        scaflag = '-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_gf_lp64 -lmkl_lapack95_lp64 ',
        knrcc = 'gcc',
        xline = '-ffixed-line-length-132 -ffree-line-length-0',
    ),

    intel = dict(
        cc = 'icc',
        cxx = 'icpc',
        fc = 'ifort',
        omp = '-fopenmp',
        cflags = dict(
            dbg = '-O0 -g',
            opt = '-O3 -debug minimal -march=native', # -xAVX
            opg = '-O3 -debug minimal -march=native -traceback -fp-stack-check -fp-speculation safe',
        ),
        fflags = dict(
            dbg = '-O0 -g -traceback -fp-stack-check -init=snan -init=arrays', # -check all  -heap-arrays 512
            opt = '-O3 -debug minimal -march=native',
            opg = '-O3 -g -march=native -traceback -fp-stack-check -init=snan -init=arrays -fp-speculation safe',
        ),
        modflag = '-module',
        ldflags = '-L${MKLROOT}/lib/intel64 -lfftw3xf_intel -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread ',
        scaflag = '-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 ',
        knrcc = 'icc',
        xline = '-extend-source 132',
    ),

    cray = dict(
        fc = 'ftn',
        cc = 'cc',
        cxx= 'CC',
        knrcc = 'gcc',
        cflags = dict(
            dbg = '-O0 -g',
            opt = '-O3 -debug minimal',
        ),
        fflags = dict(
            dbg = '-O0 -g',
            opt = '-O3 -debug minimal',
        ),
        omp = '-fopenmp',
        modflag = '-module',
        ldflags = '',
        xline = '-extend-source 132',
    ),


    env = dict(
        fc = getenv('FC'),
        cc = getenv('CC'),
        cxx = getenv('CXX'),
        knrcc = getenv('CC'),
        cflags = dict(
            dbg = '-g -O0',
            opt = '-O3',
        ),
        fflags = dict(
            dbg = '-g -O0',
            opt = '-O3',
        ),
        omp = '-fopenmp',
        modflag = '-J' if getenv('FC').endswith('gfortran') else '-module',
        ldflags = getenv('LDFLAGS'),
        xline = '-ffixed-line-length-132' if getenv('FC').endswith('gfortran') else '-extend-source 132'
    ),
)

xtra = dict(
    nvcc = dict(
        dbg = 'nvcc -O0 -arch=sm_35 --ftz false --prec-div true --prec-sqrt true --generate-line-info -g -G',
        opt = 'nvcc -O3 -arch=sm_35 --ftz true --prec-div false --prec-sqrt false', ## --relocatable-device-code true
        opg = 'nvcc -O3 -arch=sm_35 --ftz true --prec-div false --prec-sqrt false --generate-line-info -g -G',
        libs = '-L${CUDADIR}/lib64 -lcublas -lcudart -lstdc++ -L${MAGMADIR}/lib64 -lmagma -lmagmablas -ldl ',
    ),

    #ppflags = '-dunix -dAUTO-ARRAY -dF90 -dNOQUAD -dBLAS -dBLAS3 -dLAPACK -dFFTW -dFFTW3 -dINTEL_IFORT -dFORTRAN2003 -dMPIK -dSCALAPACK -dELPA -dCUDA -dLS++'.split(),
    ppflags = '-dMPIK -dSCALAPACK -dELPA -dCUDA -dLS++'.split(),
    gwflags = '-DEXPAND_ISWAP -DEXPAND_VDV -DCOMMONLL -DUSE_X0KBLAS -DX0KBLAS_DIV -DMbytes_X0KBLAS_DIV=2 -DNWORD_RECORDSIZE=1 -DEXPAND_SORTEA'
)


#        -Wl,-rpath,${MKLROOT}/lib/intel64

#           -Wl,--start-group \
#           ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a \
#           ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a \
#           ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \
#           ${MKLROOT}/lib/intel64/libmkl_sequential.a \
#           ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
#           ${MKLROOT}/lib/intel64/libmkl_core.a \
#           -Wl,--end-group -lpthread



def get_conf(args):
    d = dict(
        cud = 'cuda' in args,
        mpi = (['']+[a for a in args if a.endswith('mpi')])[-1],
        #mod = 'dbg' if 'dbg' in args else 'opt',
        mod = 'opt',
        plt = 'gcc',
        elp = 'elpa' in args,
        lsx = 'ls++' in args,
        static = 'static' in args,
        native = 'native' in args,
        mkl = 'mkl' in args,
        knl = 'knl' in args,
        openblas = 'openblas' in args,
        sca = 'scalapack' in args,
        cov = 'cov' in args,
        trace = 'trace' in args,
        prefix = (['/opt/lm'] + [a[len('prefix='):] for a in args if a.startswith('prefix=')])[-1],
        dissect = 'dissect' in args,
        as_needed = 'as-needed' in args,
        idbc = 'idbc' in args,
        nohdf5 = 'nohdf5' in args,
        pic = 'pic' in args,
        iomp = 'iomp' in args,
    )

    if 'dbg' in args: d['mod'] = 'dbg'
    if 'opg' in args: d['mod'] = 'opg'

    plt = set(plts.keys()) & set(args)

    assert(len(plt) < 2) # more than 1 platforms passed

    if len(plt) == 0:
        if getenv('FC') == '' or getenv('CC') == '':
            print('Warning: no platform specified.. ${FC} and ${CC} are empty.. defaulting to gcc',file=sys.stderr)
            d['plt'] = 'gcc'
        else:
            print('Warning: no platform specified will try to use env variables FC, CC, LDFLAGS, etc..',file=sys.stderr)
            d['plt'] = 'env'
    else:
        d['plt'] = plt.pop()

    d['mkl'] = d['mkl'] or (d['knl'] and d['plt'] == 'cray')
    d['static'] = d['static'] or d['plt'] == 'cray'

    if d['mkl'] and d['openblas']:
        print('Warning: conflicting blas implementations specified mkl and openblas',file=sys.stderr)

    return d

def uniq(l):
# use 'set' when order does not matter
    ul = []
    for i in l:
        if i not in ul:
            ul.append(i)
    return ul

def libcontains_ff(fl,fn):
    if not os.path.exists(fl): return False
    cmd = 'nm --defined-only --demangle '
    if fl.endswith('.so'): cmd += '--dynamic '
    s = xc(cmd + fl)[0]
    return fn.encode() in s # will give lots of false positives

def libcontains_lf(L,l,fn):
    f = L+'/lib'+l[2:]
    r = False
    r |= libcontains_ff(f+'.so',fn)
    r |= libcontains_ff(f+'.a',fn)
    return r

def libexists(L,l):
    f = L+'/lib'+l[2:]+'.'
    return os.path.exists(f+'so') or os.path.exists(f+'a')

def ppsl(l):
    return ', '.join(l)[::-1].replace(' ,', ' and '[::-1])[::-1]

def findlibspath(name='', libs='', incs='', env=os.environ, v=True, report_failure=True, pkgconfig_names=''):

    incpaths = []
    libpaths = []
    if 'CPATH' in env: incpaths.extend(env['CPATH'].split(':'))
    if 'INCLUDE' in env: incpaths.extend(env['INCLUDE'].split(':'))
    if 'C_INCLUDE_PATH' in env: incpaths.extend(env['C_INCLUDE_PATH'].split(':'))
    incpaths.extend('/usr/include /usr/local/include /opt/local/include'.split())

    if 'LIBRARY_PATH' in env: libpaths.extend(env['LIBRARY_PATH'].split(':'))
    if 'LD_RUN_PATH' in env: libpaths.extend(env['LD_RUN_PATH'].split(':'))
    if 'LD_LIBRARY_PATH' in env: libpaths.extend(env['LD_LIBRARY_PATH'].split(':'))
    if 'LIBPATH' in env: libpaths.extend(env['LIBPATH'].split(':'))
    libpaths.extend('/usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib /opt/local/lib /usr/lib/x86_64-linux-gnu'.split())

    incpaths = [p.rstrip('/') for p in incpaths]
    libpaths = [p.rstrip('/') for p in libpaths]

    libs = libs.split()
    incs = incs.split()

    fincpaths = [p for p in incpaths if all(os.path.exists(p+'/'+i) for i in incs)][:1] if incs != [] else []
    flibpaths = [p for p in libpaths if all(libexists(p,l)          for l in libs)][:1] if libs != [] else []

    #fincpaths = uniq(fincpaths)
    #flibpaths = uniq(flibpaths)

    if len(incs) > 0:
        if len(fincpaths) == 0 and len(flibpaths) == 0:
            for pcn in pkgconfig_names.split():
                s = xc('pkg-config --exists ' + pcn)
                if s[2] == 0:
                    fincpaths = [i[2:] for i in xc('pkg-config --cflags-only-I ' + pcn)[0].split() if os.path.exists(i[2:])]
                    flibpaths = [l[2:] for l in xc('pkg-config --libs-only-L ' + pcn)[0].split() if os.path.exists(l[2:])]
                    break


    # TODO: this (to the end of the function) needs further cleanup
    if v:
        if len(incs) > 0:
            if len(fincpaths) == 0:
                if report_failure: print('could not find '+ppsl(incs), file=sys.stderr)
            else:
                print('found '+ppsl(incs) + ' in ' + ppsl(fincpaths), file=sys.stderr)

        if len(libs) > 0:
            if len(flibpaths) == 0:
                if report_failure: print('could not find ' + ppsl(libs), file=sys.stderr)
            else:
                print('found ' + ppsl(libs) + ' in ' + ppsl(flibpaths), file=sys.stderr)


    p = uniq([i.rsplit('/',1)[0] for i in fincpaths] + [l.rsplit('/',1)[0] for l in flibpaths])

    prefix = ''
    rincs = ''
    rlibs = ''

# this is nice but if plenty of system libs are used it ends up with a long list of the same prefixes for diff libs eventually.
    #if len(p) == 1:
        #prefix = name + ' = ' + p[0]
        #rincs = ' '.join([('-I${'+name+'}/'+i.rsplit('/',1)[1]) for i in fincpaths])
        #rlibs = ' '.join(['-L${'+name+'}/'+l.rsplit('/',1)[1] for l in flibpaths])
    #elif len(p) > 1:
        #print('could not find includes '+', '.join(incs) + ' and libs ' + ', '.join(libs) + ' under a single roof', file=sys.stderr)
        #rincs = ' '.join(['-I'+i for i in fincpaths])
        #rlibs = ' '.join(['-L'+l for l in flibpaths])

    if len(p) > 0:
        rincs = ' '.join(['-I'+i for i in fincpaths])
        rlibs = ' '.join(['-L'+l for l in flibpaths])
        if len(p) == 1:
            prefix = name + ' = ' + p[0]
        if (len(p) > 1):
            skip_warn = False
            # "modern" prefixes found in suse and ubuntu
            mpref_i = '-I/usr/include'
            mpref_l = '-L/usr/lib/x86_64-linux-gnu'
            if rincs.startswith(mpref_i) and rlibs.startswith(mpref_l): skip_warn = rincs[len(mpref_i):] == rlibs[len(mpref_l):]
# fedora's weirdness
            skip_warn |= rincs == '-I/usr/lib64/gfortran/modules' and rlibs == '-L/usr/lib64'
            skip_warn |= rincs == '-I/usr/include/openmpi-x86_64' and rlibs == '-L/usr/lib64/openmpi/lib'
            if not skip_warn: print('could not find includes '+', '.join(incs) + ' and libs ' + ', '.join(libs) + ' under a single roof', file=sys.stderr)

    return prefix,rincs,rlibs


def add_inclib_from_paths(paths,env):
    for p in paths:
        if p in env:
            env['CPATH'] = env[p]+'/include' + ((':' + env['CPATH']) if 'CPATH' in env else '')
            env['LIBRARY_PATH'] = ':'.join([env[p]+i for i in '/lib /lib64'.split()]) + ((':' + env['LIBRARY_PATH']) if 'LIBRARY_PATH' in env else '')




def gen_fldict(plt='gcc', mod='opt', mpi='', prefix='/opt/lm', cud=False, elp=False, lsx=False, static=False, native=False, mkl=False, sca=False, openblas=False, knl=False, cov=False, trace=False, dissect=False, as_needed=False, idbc=False, nohdf5=False, pic=False, iomp=False):
    import os, re

    ppflags = xtra['ppflags']
    if not cud:
        for k in ppflags:
            if k.startswith('-dCUDA'): ppflags.remove(k)
    if mpi == '' and plt!='cray':
        for k in ppflags:
            if k.startswith('-dMPI'): ppflags.remove(k)
        ppflags.remove('-dSCALAPACK')
        ppflags.remove('-dELPA')
        ppflags.remove('-dLS++')
        sca = False


    if not elp:
        if '-dELPA' in ppflags: ppflags.remove('-dELPA')
    if not lsx:
        if '-dLS++' in ppflags: ppflags.remove('-dLS++')

    d = dict(
        cc      = plts[plt]['cc'],
        cflags  = plts[plt]['cflags'][mod],
        fc      = plts[plt]['fc'],
        fflags  = plts[plt]['fflags'][mod] + ' ' + plts[plt]['xline'],
        modflag = plts[plt]['modflag'],
        ldflags = plts[plt]['ldflags'],
        nvcc    = xtra['nvcc'][mod],
        cxx     = plts[plt]['cxx'],
        ppflags = ' '.join(xtra['ppflags']),
        knrcc   = plts[plt]['knrcc'],
        gwflags = plts[plt]['xline'] + ' ' + xtra['gwflags'],
        omp     = plts[plt]['omp'],
        libxcpath = '',
        special = {},
        header = '',
        prefix = prefix,
    )

    if mkl and plt=='cray': d['ldflags'] = plts['intel']['ldflags']

    if not native:
        d['fflags'] = re.sub(r'(?ix)\-march=native\s?','',re.sub(r'(?ix)\-x\s*host\s?','',d['fflags']))
        d['cflags'] = re.sub(r'(?ix)\-march=native\s?','',re.sub(r'(?ix)\-x\s*host\s?','',d['cflags']))
        if knl:
            d['cflags'] += ' -march=knl'
            d['fflags'] += ' -march=knl'
            # to be removed when intel 18 is old
            if plt in 'intel cray'.split():
                d['cflags'] = re.sub(r'(?ix)\-march=knl\s?','-xMIC-AVX512',d['cflags'])
                d['fflags'] = re.sub(r'(?ix)\-march=knl\s?','-xMIC-AVX512',d['fflags'])

    if sys.platform.startswith('darwin'):
        d['cflags'] = re.sub(r'(?ix)\s*\-fstack\-protector\-all','',d['cflags'])
        d['fflags'] = re.sub(r'(?ix)\s*\-fstack\-protector\-all','',d['fflags'])

    if plt == 'env': return d

    incs = ''
    header = []


    env = os.environ.copy()
    add_inclib_from_paths('SCALAPACKPATH SCALAPACKROOT SCALAPACKDIR SCALAPACK_PATH SCALAPACK_ROOT SCALAPACK_DIR SCALAPACKBASE SCALAPACK_BASE SCALAPACKHOME SCALAPACK_HOME'.split(), env)
    if mpi == 'openmpi': env['LIBRARY_PATH'] = '/usr/lib64/openmpi/lib' + ((':' + env['LIBRARY_PATH']) if 'LIBRARY_PATH' in env else '')

    scalib = '-lscalapack'
    scapref,scainc,scaLib = findlibspath(name='SCALAPACK_ROOT', libs=scalib, env=env, v=True, report_failure=False, pkgconfig_names='scalapack scalapack-openmpi')
    if scaLib == '':
        scalib = '-lscalapack-openmpi'
        scapref,scainc,scaLib = findlibspath(name='SCALAPACK_ROOT', libs=scalib, env=env, v=True, report_failure=False, pkgconfig_names='scalapack scalapack-openmpi')
    d['scapath'] = ''

    mklroot = (['']+[os.environ[v] for v in 'MKL_HOME MKLROOT'.split() if v in os.environ]).pop()
    mklavail = (mklroot != '' or mkl) and (not openblas)
    #mklpref,mklinc,mkllib = findlibspath(name='MKLROOT', libs='-lmkl_core', incs='mkl_services.mod')

    if mklavail: sca = True

#Only leave -lfftw3xf_* in ldflags if MKL is from before 2014. Otherwise remove it since it is included in the base lib.
    if mklavail:
        mkl_y = min(map(int,re.findall(r'20\d\d',mklroot+'2099')))
        if not (mkl_y < 2014):
            d['ldflags'] = re.sub(r'-lfftw3xf_\w+\s','',d['ldflags'])
        if plt == 'gcc' and mkl_y < 2016:
            d['fflags'] += ' -fcray-pointer'
        if iomp:
            d['ldflags'] = re.sub(r'-lgomp\b','-liomp5 -L${MKLROOT}/../compiler/lib/intel64',d['ldflags'])
    else:
        if plt != 'cray' and plt != 'env':
            d['ldflags'] = ''

            env = os.environ.copy()

            env['CPATH'] = ':'.join([env[i] for i in 'FFTW3_INCLUDE FFTW3_INC FFTW_INCLUDE FFTW_INC'.split() if i in env]) + ((':' + env['CPATH']) if 'CPATH' in env else '')
            env['LIBRARY_PATH'] = ':'.join([env[i] for i in 'FFTW3_LIB FFTW3_LIBPATH FFTW_LIB FFTW_LIBPATH'.split() if i in env]) + ((':' + env['LIBRARY_PATH']) if 'LIBRARY_PATH' in env else '')

            fftw3pref,fftw3inc,fftw3lib = findlibspath(name='FFTW3PATH', libs='-lfftw3', incs='fftw3.f03', env = env, pkgconfig_names='fftw3')
            if fftw3lib != '':
                #if fftw3pref != '/': header.append(fftw3pref)
                incs += ' '+fftw3inc
                d['ldflags'] += ' '+fftw3lib
                d['ldflags'] += ' -lfftw3'

            env = os.environ.copy()
            if sys.platform.startswith('darwin'):
                env['LIBRARY_PATH'] = ((env['LIBRARY_PATH'] + ':') if 'LIBRARY_PATH' in env else '') + '/usr/local/opt/openblas/lib'

            obll = ''
            openblaspref=''
            openblasinc = ''
            openblaslib = ''
            for i in '-lopenblas_openmp -lopenblaso -lopenblas'.split():
                openblaspref,openblasinc,openblaslib = findlibspath(libs=i, v=openblas, env=env, report_failure=False, pkgconfig_names='openblas openblas-openmp')
                obll = ' ' + i
                if openblaslib != '': break

            if openblaslib != '':
                if not sys.platform.startswith('darwin'):
                    if not libcontains_lf(openblaslib,obll,'zgetri'): d['ldflags'] += ' -llapack'
                d['ldflags'] += ' ' + openblaslib + obll
            else:
                lpkpref,lpkinc,lpklib = findlibspath(name='LAPACKPATH', libs='-llapack -lblas')
                #header.append(lpkpref)
                if lpklib != '': d['ldflags'] += ' ' + lpklib
                d['ldflags'] += ' -llapack -lblas'

            if sca: plts[plt]['scaflag'] = '-lscalapack'
            if scaLib != '':
                plts[plt]['scaflag'] = scaLib + ' ' + scalib
                #header.append(scapref)


        else:
            d['ldflags'] = re.sub(r'-lfftw3xf_\w+\s','',d['ldflags'])


    if re.search(r'-l?mkl', d['ldflags'] + d['fflags']) != None :
        d['gwflags'] += ' -DMKL'
        if re.search(r'-lmkl', d['ldflags']) != None and re.search(r'-mkl', d['fflags']) == None:
            ##d['fflags'] += ' -I${MKLROOT}/include/intel64/lp64'
            #incs += ' -I${MKLROOT}/include/intel64/lp64'
            incs += ' -I${MKLROOT}/include'
            header.append('MKLROOT = '+mklroot)


    #print (d['ldflags'])

    env = os.environ.copy()

    env['CPATH'] = ':'.join([env[i] for i in 'LIBXC_INCLUDE LIBXC_INC'.split() if i in env]) + ((':' + env['CPATH']) if 'CPATH' in env else '')
    if plt == 'gcc': env['CPATH'] += ':/usr/lib64/gfortran/modules'

    env['LIBRARY_PATH'] = ':'.join([env[i] for i in 'LIBXC_LIB LIBXC_LIBPATH'.split() if i in env]) + ((':' + env['LIBRARY_PATH']) if 'LIBRARY_PATH' in env else '')

    add_inclib_from_paths('LIBXCPATH LIBXCROOT LIBXCDIR LIBXCBASE LIBXC_PATH LIBXC_ROOT LIBXC_DIR LIBXC_BASE'.split(), env)

    xclib = '-lxcf90 -lxc'
    xcpref,xcinc,xcLib = findlibspath(libs=xclib, incs='xc_f90_lib_m.mod', env = env, pkgconfig_names='libxc')
    if xcLib == '':
        xclib = '-lxc'
        xcpref,xcinc,xcLib = findlibspath(libs=xclib, incs='xc_f90_lib_m.mod', env = env, pkgconfig_names='libxc')

    if xcLib != '':
        #header.append(xcpref)
        incs += ' ' + xcinc
        d['ldflags'] += ' '+xcLib
    d['ldflags'] += ' ' + xclib

    d['hdf5path'] = ''

    if not (plt == 'cray' or nohdf5):
        hdf5env = os.environ.copy()
        if mpi == 'openmpi':
            hdf5env['CPATH'] = '/usr/include/openmpi-x86_64:/usr/include/hdf5/openmpi' + ((':' + hdf5env['CPATH']) if 'CPATH' in hdf5env else '')
            hdf5env['LIBRARY_PATH'] = '/usr/lib/x86_64-linux-gnu/hdf5/openmpi' + ((':' + env['LIBRARY_PATH']) if 'LIBRARY_PATH' in env else '')

        hdf5env['CPATH'] = ':'.join([hdf5env[i] for i in 'HDF5_INCLUDE HDF5_INC'.split() if i in hdf5env]) + ((':' + hdf5env['CPATH']) if 'CPATH' in hdf5env else '')
        hdf5env['LIBRARY_PATH'] = ':'.join([hdf5env[i] for i in 'HDF5_LIB HDF5_LIBPATH'.split() if i in hdf5env]) + ((':' + hdf5env['LIBRARY_PATH']) if 'LIBRARY_PATH' in hdf5env else '')

        add_inclib_from_paths('HDF5PATH HDF5ROOT HDF5_PATH HDF5_ROOT H5PATH H5ROOT H5_PATH H5_ROOT HDF5DIR HDF5_DIR H5DIR H5_DIR HDF5BASE HDF5_BASE'.split(), hdf5env)

        hdf5pref,hdf5inc,hdf5lib = findlibspath(name='HDF5_ROOT', libs='-lhdf5', incs='hdf5.h', env=hdf5env, pkgconfig_names='hdf5 hdf5-openmpi')

        #if hdf5inc != '': d['fflags'] += ' '+hdf5inc
        if hdf5lib != '':
            #header.append(hdf5pref)
            incs += ' ' + hdf5inc
            d['ldflags'] += ' '+hdf5lib
        d['ldflags'] += ' -lhdf5'
        d['gwflags'] += ' -DUSE_HDF5'

    mpigw = False
    #mpigw = True

    if mpigw:
        if mpi != '' or plt == 'cray': d['gwflags'] += ' -DUSE_MPI'
        if sca or scaLib != '': d['gwflags'] += ' -DUSE_SCALAPACK'

    if static: d['ldflags'] = '-static -Wl,--start-group '  + d['ldflags'] + ' -ldl -Wl,--end-group'

    if mpi != '' or plt == 'cray':
        if plt != 'cray':
            if mpi == 'openmpi':
                d['cc'] = 'mpicc'
                d['fc'] = 'mpif90'
                d['cxx'] = 'mpicxx'
                if sca or scaLib != '':
                    d['ldflags'] = plts[plt]['scaflag'] + ' ' + d['ldflags']
                else:
                    d['ppflags'] = d['ppflags'].replace('-dSCALAPACK','')

            elif mpi == 'intelmpi':
                d['cc'] = 'mpiicc'
                d['fc'] = 'mpiifort'
                d['cxx'] = 'mpiicpc'
                d['ldflags'] = plts[plt]['scaflag'].replace('openmpi','intelmpi') + ' ' + d['ldflags']
            d['mpirun'] = 'mpirun'

            if dissect:
                ptr = re.compile(r'\s+-(?:Wl|L|l)[^\s]*')
                c = d['fc'].split(None,1)[0]; d['fc'] = d['fc'].split(None,1)[1:]
                if (d['fc'] == []): d['fc'] = ''
                s = xc(c + ' -show')[0].decode()
                d['ldflags'] = ' '.join(ptr.findall(s) + d['ldflags'].split())
                d['fc'] = ptr.sub('',s).strip() + ' ' + d['fc']

        else:
            d['mpirun'] = 'aprun'
    else:
        d['mpirun'] = 'env'

    if cud:
        d['ldflags'] += xtra['nvcc']['libs']

    if d['ldflags'] != '':
        d['ldflags'] = ' '.join(uniqf(d['ldflags'].split(),f=lambda i: i.startswith('-L')))

    if as_needed:
        #d['ldflags'] = '-Wl,-copy-dt-needed-entries -Wl,-as-needed ' + d['ldflags']
        d['ldflags'] = '-Wl,-as-needed ' + d['ldflags']

    if idbc:
        d['fflags'] += ' -gdwarf-3'
        d['cflags'] += ' -gdwarf-3'

    if pic:
        d['fflags'] += ' -fpic' # -fPIC for unlimited table size on some arches
        d['cflags'] += ' -fpic'

    if incs != '':
        incs = ' ' + ' '.join(uniqf(incs.split(),f=lambda i: i.startswith('-I')))
        d['fflags'] += incs
        incs = re.sub(r'(?ix)\-I\s*/include(?:\s|$)','',incs)
        incs = re.sub(r'(?ix)\-I\s*/usr/include(?:\s|$)','',incs)
        d['cflags'] += incs

    if cov:
        d['fflags'] += ' --coverage -fprofile-update=atomic'
        d['cflags'] += ' --coverage -fprofile-update=atomic'

    d['header'] = '\n'.join(header)

    d['special'] = {}

    if (plt == 'intel' or plt == 'cray') and mod[:2] == 'op':
        d['special'] = {
            'v7input/m_rdctrl.f.o'    : '-O1',
            'v7input/m_rdctrlchk.f.o' : '-O1',
            'subs/strxops.f.o'        : '-O2',
            'subs/spcgrp.f.o'         : '-O2',
            #'subs/dosmsh.f.o'        : '-O2',
            #'subs/defwsr.f.o'        : '-O2',
            'fp/smvxcm.f.o'           : '-O2',
            'fp/rhogkl.f.o'           : '-O2',
            'fp/rhopos.f.o'           : '-O2',
            'fp/flocbl.f.o'           : '-O2',
            'fp/vxcnlm.f.o'           : '-O2',
            'gf/mkcpa.f.o'            : '-O2',
            'gf/mkptfp.f.o'           : '-O2',
            'mol/rhoti0.f.o'          : '-O2',
            'subs/suctrl.f.o'         : '-O3', # observed to break blm test on sandy bridge xeon with ifort 18.0.2 -march=native
        }

        if native:
            d['special']['subs/rdeq_old.f.o'] = '-O3'
            d['special']['subs/rseq.f.o'    ] = '-O3'
            d['special']['pgf/atomsr.f.o'   ] = '-O3'
            #d['special']['slatsm/dfphi.f.o' ] = '-O3'

        if knl:
            d['special']['fp/rsibl.f.o'] = '-O3' #causes ICE on archers ifort 17 with -xMIC-AVX512

        if (plt == 'intel'):
            s = xc('ifort -v')[1]
            if (b'ifort version 18.' in s):
                d['special']['nc/hmfr3c.f.o'] = '-O2'
               #d['special']['slatsm/dfphi.f.o' ] = '-O3'

    elif (plt == 'gcc'):
        ##gcc before 7 is not supported anymore anyway
        #s = xc('gfortran -v')[1]
        #gccv = re.findall(r'gcc\s+version\s+(\d+)', s)
        #if gccv < '7':
            #d['fflags'] = re.sub(r'(?ix)\-Wno\-argument\-mismatch\s?','',d['fflags'])

        if (mod[:2] == 'op'):
            d['special'] = {
                'pgf/pgfasa.f.o'          : '-O1',
                'subs/hsmq.f.o'           : '-O2', # v8.1.0+
            }

            if native or knl:
                d['special']['gf/gfg2g.f.o'] = '-O3'
                d['special']['gf/emesh.f.o'] = '-O3'
                d['special']['slatsm/fftz3.f.o'] = '-O3' # v8.1.1

    for k in d['special'].keys():
        ##no need of all this since special flags only apply to op*
        #d['special'][k] = re.sub(r'-O(\d)',
            #lambda mo: mo.group(0) if int(mo.group(1)) <= int(re.search(r'-O(\d)',d['special'][k]).group(1)) else d['special'][k], d['fflags'])
        d['special'][k] = re.sub(r'-O\d', d['special'][k], d['fflags'])
        d['special'][k] = re.sub(r'(?ix)\-x\s*host\s?','',d['special'][k])
        d['special'][k] = re.sub(r'(?ix)\-march=native\s?','',d['special'][k])
        d['special'][k] = re.sub(r'(?ix)\-march=knl\s?','',d['special'][k])
        if plt == 'cray':
            #d['special'][k] += ' -march=generic' #not supported by ifort 17 on archer knl
            d['special'][k] += ' -xSSE2'

    if (plt == 'intel' or plt == 'cray') and not nohdf5: d['special']['slatsm/h5.F90.o'] = '$(fflags) -allow nofpp_comments'

# this enables the nice runtime trace from random error code exits. shall not be treated as the other specials above.
    if 'slatsm/fmain.F90.o' not in d['special']: d['special']['slatsm/fmain.F90.o'] = '$(fflags)'
    if mod[2] == 'g' or trace: d['special']['slatsm/fmain.F90.o'] += ' -DTRACE'
    if cov: d['special']['slatsm/fmain.F90.o'] += ' -DTESTCOV=%d' % os.path.abspath('.').count('/') # this will break if not in the same directory as the output file..
    if d['special']['slatsm/fmain.F90.o'] == '$(fflags)': del d['special']['slatsm/fmain.F90.o']

    #no bounds check files
    #subs/hamfb3.f subs/roth.f subs/bloch.f subs/mksbet.f subs/hmltns.f subs/pqmix.f
    # fp/augmbl.f fp/dfrce.f fp/sopert.f

    d['qcmd'] = 'env'
    d['qhdr'] = ''

    return d

def fldict2makefl(d):
    s = '''
{header}

cc      = {cc}
cflags  = {cflags}
fc      = {fc}
fflags  = {fflags}

# The flag to set destination path for .mod files.
modflag = {modflag}

# The flag enabling OpenMP.
omp     = {omp}

# The nvidia/cuda compiler. Not necessary unless you want to play with it.
nvcc    = {nvcc}

# All general linking related flags (library and run paths, libs etc...)
ldflags = {ldflags}

# The C++ compiler to be used if necessary for glue code.
cxx     = {cxx}

# Basic KNR C compiler used only to compile ccomp. Beware the Cray architecture! Its default C compiler (cc) produces binaries which do not usually work on the compilation/login/service nodes and ccomp is needed during compilation.
knrcc   = {knrcc}

# Preprocessing flags for the custom fortran preprocessor "ccomp", used for all .f, .for and .F* files except for those under subdirectory fpgw. The compiler internal preprocessor is used for the QSGW files.
ppflags = {ppflags}

# Flags for the QSGW package in subdirectory fpgw.
gwflags = {gwflags}

# Where the 'install' target will send files.
prefix = {prefix}

# Command to use to run test jobs, NPROC and NAME will be replaced appropriately if specified. Ensure the command does not switch directory.
#qcmd = {qcmd}

# Header to include in the test job scripts, useful to load modules, set library paths etc... use \\\\n and \$$ to escape \\n and $ if needed.
#qhdr = {qhdr}

'''.format(**d)

    if 'mpirun' in d: s = s + '\n# Launcher for MPI linked binaries.\nmpirun = {mpirun}\n'.format(**d)

    s += '\n'*2

## a shoddy replacement for .format() on ancient environments:
    #for k,v in d.items(): s = s.replace('{'+k+'}',str(v))

    if len(d['special'].items()) > 0:
        sflags = sorted(set(d['special'].values()))
        #sfiles = [[] for f in sflags]
        sfiles = list(map(lambda f: [], sflags))
        for k,v in d['special'].items():
            sfiles[sflags.index(v)].append(k)

        ls = ['# Flags for routines with excessive branching which can choke a compiler but are not crucial to overall performance.']
        for i in range(len(sflags)):
            ls.append('''
lessflags{level} = {sflags}
lessfiles{level} = {sfiles}
'''.replace('{level}',str(i)).replace('{sflags}',sflags[i]).replace('{sfiles}',' '.join(sfiles[i])))  #.format(level = i, sflags = sflags[i], sfiles = ' '.join(sfiles[i])))
        s = s + '\n'.join(ls)

    return s


def shortuse(display_flags):
    #return ' usage: '+sys.argv[0]+' [%s] [dbg|opt] [openmpi|intempi] [cuda] [elpa] [ls++]'%('|'.join(plts.keys()-set(['env'])))
    return ' usage: '+sys.argv[0]+' prefix=path/to/install/to ' + display_flags
def usage(display_flags):
    s = '''# help:
#
#   key: flags in [] are optional, flags separated by "|" are incompatible with each other.
#
#  %s > flags.mk
#        path/to/lm/configure.py
#        ninja [-j nprocs] [build_targets]
#
#''' %  shortuse(display_flags)

    return s




if __name__ == '__main__':
    import os, sys

    display_flags = '[intel|gcc|cray|env] [opt|opg|dbg] [openmpi|intelmpi] [native|knl] [static] [mkl|openblas] [cuda] [scalapack] [elpa] [ls++] [cov] [trace] [as-needed] [pic]'
    all_flags = display_flags + ' [dissect] [idbc] [nohdf5] [iomp]'
    help_flags = set('-h --h --help -help usage help shortuse'.split())
    ref_flags = set(all_flags.replace('[',' ').replace(']',' ').replace('|',' ').split()) | help_flags
    pas_flags = set([a for a in sys.argv[1:] if not a.startswith('prefix=')])

    bad_flags = pas_flags - ref_flags
    if bad_flags != set():
        print('\nUnrecognised args passed: ' + ' '.join(bad_flags) + '\n',file=sys.stderr) # + '\n\nHere is how to call me:\n'

    if 'shortuse' in sys.argv:
        print(shortuse(display_flags))
        sys.exit(0)
    elif (help_flags & pas_flags) != set() or bad_flags != set():
        print(usage(display_flags),file=sys.stderr)
        sys.exit(0)


    d = get_conf(sys.argv)

    #print(usage(display_flags))
    print('''#
# the following was produced by:
#     {args}
#
'''.replace('{args}',' '.join(sys.argv))) #.format(args    = ' '.join(sys.argv)))

    print(fldict2makefl(gen_fldict(**d)))

