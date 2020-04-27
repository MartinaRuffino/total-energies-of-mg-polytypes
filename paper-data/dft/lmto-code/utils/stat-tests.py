#!/usr/bin/env python

import os, sys, re
from glob import glob as glob
from datetime import datetime

regime_ancien = False
try:
    import argparse
except ImportError:
    regime_ancien = True

if sys.version_info[0] < 3:
    from codecs import open

def loadfile(fn):
    f = open(fn,mode='r',encoding='utf8')
    s = f.read()
    f.close()
    return s

tcols = dict(
    endc = '\x1b[0m',  # stop a colour

    gry = '\x1b[30m',
    red = '\x1b[31m',
    grn = '\x1b[32m',
    yel = '\x1b[33m',
    blu = '\x1b[34m',
    vio = '\x1b[35m',
    cyn = '\x1b[36m',
    wht = '\x1b[37m',

    # bold
    gryb = '\x1b[1;30m',
    redb = '\x1b[1;31m',
    grnb = '\x1b[1;32m',
    yelb = '\x1b[1;33m',
    blub = '\x1b[1;34m',
    viob = '\x1b[1;35m',
    cynb = '\x1b[1;36m',
    whtb = '\x1b[1;37m'
)

def getstime(s):
    stime = int(datetime.now().strftime('%s'))
    try:
        stime = int(s.split('\n',1)[0].split()[1])
    except IndexError:
        pass
    return stime

def check_test(fln,nocol=False):
    log = fln
    #nam = fln[7:-4]
    nam = fln
    s = loadfile(log).strip()
    sl = s.lower()

    #stat = 'in progress'
    stat = 'running'
    stime = getstime(sl)
    wtime = ''
    if 'walltime:' in sl:
        stat = {True:'fail', False:'pass'}[
            'FAIL' in s
         or 'PASS' not in s
         or 'aborted' in sl
         or 'backtrace' in sl
         or 'command not found' in sl
         or 'error while loading shared' in sl
         or 'job returned with error status' in sl
        ]
        wtime = re.findall(r'walltime:\s*(\d*\.?\d*\s*s)',sl.rsplit('\n',1)[1])[0]
    else:
        wtime = '%d s' % (int(datetime.now().strftime('%s')) - stime)
    stat = stat.ljust(12)

    if not nocol:
        #col = {'fail':'redb', 'pass':'grnb', 'in progress':'yelb'}[stat.strip()]
        col = {'fail':'redb', 'pass':'grnb', 'running':'yelb'}[stat.strip()]
        stat = tcols[col] + stat + tcols['endc']

    #return stime,('  %s %s %7s' % (nam.ljust(32),stat.ljust(11),wtime))
    print('  %s %s %7s' % (nam.ljust(32),stat,wtime))

def quickstime(fln):
    fl = open(fln)
    sl = fl.readline()
    fl.close()
    stime = getstime(sl)
    return stime

if __name__ == '__main__':
    if regime_ancien:
        #print ('# no options available in legacy mode')
        class args_t(object): pass
        args = args_t()
        args.logfiles = '*'
        args.nots = False
        args.nocol = False
    else:
        aprs = argparse.ArgumentParser(description='Tests log files status checker.')
        aprs.add_argument('--nocol', action='store_true', help='no colours')
        aprs.add_argument('--nots', action='store_true', help='do not sort by start time')
        aprs.add_argument('logfiles', metavar='logfile', help="log file(s)", nargs='*')
        args = aprs.parse_args()

    flns = list(filter(os.path.isfile, args.logfiles))

    if flns == []:
        flns = filter(os.path.isfile, glob('checks/*.log'))
# it is nice to have the autegenerated filelist sorted by some sort of time
        if not args.nots: flns = sorted(flns, key=lambda fln: quickstime(fln))

    for fln in flns:
        check_test(fln, args.nocol or not sys.stdout.isatty())



