Massaging site data imported from as POSCAR file
================================================

It is often the case that you have site positions obtained from
another code such as VASP.  The site positions and lattice vectors can
be readily imported with the `poscar2init` and `poscar2site`
utilities.  However there are some subtleties that should be taken into account,
which this tutorial addresses.

- Find missing symmetry operations.  POSCAR files often correspond to positions
  after positions are relaxed.  With some small adjustments symmetry operations
  may appear that are not present in the raw data.

- Shorten the basis vectors.  The VASP file in this tutorial puts
  all the basis vectors so that coordinates lie between 0 and 1, with
  many close to 1.  This is allowed, but there are a number of places
  in this code which need vectors to be as short as possible, e.g. for
  Ewald sums to converge well.  Internally this is taken care of, but
  there can be complexities down the road (e.g. extra phases when
  manipulating the self-energy), that it is safer and easier to work with
  shortened basis vectors.

- Add floating orbitals when doing GW calculations on open systems.
  For open systems the augmentation spheres fill a small fraction of
  the unit cell. This is acceptable for LDA, but GW calculations
  demand a more stringent basis set, because many unoccupied states
  well above the Fermi level need to be reasonably described.  We will
  show how to improve the basis by adding floating orbitals in the
  interstitial region.  Floating orbitals do not have any augmentation
  parts, but they improve the completeness of intertitial part of the basis.


### Initial site file

Convert the POSCAR file into a site file

    poscar2site

This generates file `site`.  The first two lines should read:

    % site-data vn=3.0 xpos fast io=15 nbas=12 alat=1.889727 plat= 5.92 0.0 0.0 0.0 5.92 0.0 0.0 0.0 5.92
    C     0.8962732309     0.9999835576     0.9885225898
    ...

Note the 'xpos' in the header.  This indicates that site coordinates
are expressed in units of lattice vectors QLAT.  This POSCAR used
alat=1 Angstrom, so the site lattice vectors are rather large:
|PLAT|~5.  The LMTO codes can handle this, but various routines are
designed to operate when |PLAT|~1.  Make the change by scaling PLAT by
1/5.92, and alat by 5.92, so that the header line becomes

    % site-data vn=3.0 xpos fast io=15 nbas=12 alat=1.889727*5.92 plat= 1 0.0 0.0 0.0 1 0.0 0.0 0.0 1

Now you can run lmchk by doing:

    cp site site.mapb
    lmchk ctrl.mapb

### Finding hidden symmetry operations

There are two problems with the site file we want to fix before proceeding further.

- Shorten lattice vectors

- Find a hidden symmetry operation

lmchk displays this line:

    MKSYM:  found 1 space group operations; adding inversion generated 2 ops

But there is in fact one mirror operation.  We can kill two birds with one stone with the following:

    lmchk --shorten --fixpos:tol=2e-3 ctrl.mapb  --wsitex=site2
    cp site2.mapb site.mapb

`--fixpos:tol=2e-3` puts a "fuzz" into the site positions when checking for symmetry operations.
If a symop can be found consistent with an adjustment of the position, it will be shifted.

`--wsitex` tells lmchk to write basis vectors as multiples of QLAT.  `--wsite` will cause lmchk
to write them in Cartesian coordinates.

The output of lmchk shows that a mirror operation is now present:

    GROUPG: the following are sufficient to generate the space group:
         my

The first lines of the site file should read:

    % site-data vn=3.0 xpos fast io=62 nbas=12 alat=11.18718384 plat= 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
    #                        pos                                   vel                                    eula                   vshft  PL rlx
    C        -0.1037268  -0.0000164  -0.0114774    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111

Note e.g. that -0.1037268 = 0.8962732309-1 .  Thus the new site
positions differ from the original ones by integers. Adding integral
multiples of PLAT does not change the basis.  There is a slight
addtional shift because the basis vectors are tweaked to comport with
the mirror symmetry.

We can stop there, but inspect `site.mapb` and you will see that there
is a uniform shift of -0.0000164 along the y axis.  This is ok, but
the positions are tidier when they lie at rational multiples of the
lattice vectors.  The y axis can be shifted slightly to reflect this.
With your text editor add a constant +.0000164 to all y coordinates as
follows:

    C        -0.1037268  -0.0000164+.0000164  -0.0114774    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    N         0.1421229  -0.0000164+.0000164   0.0326819    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    H        -0.1330424  -0.0000164+.0000164  -0.1948035    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    H        -0.1776707   0.1522947+.0000164   0.0648669    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    H        -0.1776707  -0.1523276+.0000164   0.0648669    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    H         0.2211880   0.1427004+.0000164  -0.0341978    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    H         0.2211880  -0.1427333+.0000164  -0.0341978    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    H         0.1804493  -0.0000164+.0000164   0.2053792    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    Pb        0.4753318   0.4999836+.0000164   0.4779234    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    Br        0.4287771   0.4999836+.0000164  -0.0273080    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    Br        0.4334819  -0.0000164+.0000164  -0.4871517    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
    Br       -0.0324470   0.4999836+.0000164   0.4392057    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111

If you do algebra within a site file, you must eliminate the `fast` tag in the header.
With these changes run lmchk again and note site positions are a bit simpler along y.

### Improve lmto basis.

lmchk displays ratio of sphere volumes to the volume of the unit cell:

    Cell volume= 1400.11054   Sum of sphere volumes= 357.25605 (0.25516)

Only about a quarter of the unit cell is filled by augmentation
spheres.  This is acceptable for LDA, but for GW the LMTO envelope functions
are not quite complete enough.  We will fix this with some 'floating orbitals'

lmchk has an automatic generator to find "empty site" positions --- the points
farthest removed from other sites.  To find them, do:

    lmchk --findes ctrl.mapb

Inspect file poses.mapb.  21 "empty" sites were found: the sites are
ordered from largest to smallest.  We choose the first 11 sites,
because after that the sphere volumes become too small.

With your text editor, modify site.mapb by hand, using the output of
poses.mapb.  First, the header has to be modified to take into account
the extra sites being added.  We will modify it like this:

    % site-data vn=3.0 xpos io=62 nbas=12+{les?11:0} alat=11.18718384 plat= 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0

That way we can run lmchk with a variable `les`, which enables us to
include them or exclude them without modifing input files.  The
input file is set up so that

    les=0    no empty sites
    les=1    empty sites with augmentation spheres
    les=11   floating orbitals

Note that `ctrl.mapb` has this line

    ATOM=E          Z=  0  R=2.2*{les==1?1:0} LMX={lmbe}

It says that, if the empty sites are to be treated as spheres
(les==1), set the augmentation radius to 2.2.  If it is les==11, R
will be zero.  This tells the programs to treat these sites as floating orbitals
rather than as lattice positions with an augmentation sphere.

Use the `XPOS` columns extracted from poses.mapb (the site file's coordinates
are expressed as multiples of QLAT) and modify your site file to read:

    % site-data vn=3.0 xpos io=62 nbas=12+{les?11:0} alat=11.18718384 plat= 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
    #                        pos                                   vel                                    eula                   vshft  PL rlx
     C        -0.1037268  -0.0000164+.0000164  -0.0114774    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     N         0.1421229  -0.0000164+.0000164   0.0326819    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H        -0.1330424  -0.0000164+.0000164  -0.1948035    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H        -0.1776707   0.1522947+.0000164   0.0648669    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H        -0.1776707  -0.1523276+.0000164   0.0648669    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H         0.2211880   0.1427004+.0000164  -0.0341978    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H         0.2211880  -0.1427333+.0000164  -0.0341978    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H         0.1804493  -0.0000164+.0000164   0.2053792    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     Pb        0.4753318   0.4999836+.0000164   0.4779234    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     Br        0.4287771   0.4999836+.0000164  -0.0273080    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     Br        0.4334819  -0.0000164+.0000164  -0.4871517    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     Br       -0.0324470   0.4999836+.0000164   0.4392057    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.4900781   0.0000000            0.0000781    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.0729687   0.0000000            0.4213281    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.0712500   0.5000000           -0.0611719    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.2813281  -0.2256250           -0.2336719    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.2813281   0.2256250           -0.2336719    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.1100781  -0.2035156           -0.2714063    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.1100781   0.2035156           -0.2714063    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.2183594  -0.2440625            0.2259375    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.2183594   0.2440625            0.2259375    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.3161719  -0.1933594            0.2550781    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.3161719   0.1933594            0.2550781    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111

Note that poses resolves the empty sites by class, where we make them all one species.

The next steps aren't necessary, but they provide insight and
positions the empty sites a bit more uniformly in the interstitial.
With your updated `site.mapb` run `lmchk` once again

    lmchk -vles=1 ctrl.mapb

and note these line:

    Cell volume= 1400.11054   Sum of sphere volumes= 847.88067 (0.60558)

    OVMIN, 475 pairs:  fovl = 1.64497e-7   <ovlp> = 7.4%   max ovlp = 11.6%

The first line tells you that the sphere volume (using R=2.2 for the E radius)
has increased to 60%.  This will be sufficient for GW calculations.

The OVMIN line tells you that the largest overlap (using R=2.2 for the E radius)
between E and some other site is about 11%.

You can reduce this using the following:

    lmchk -vles=1 ctrl.mapb --mino~z --wsitex=site2

`--mino~z` tells `lmchk` to adjust the E positions (any site with Z=0) to reduce the overlap.
The last line now reads

    minimized: fovl = 1.5099e-9   <ovlp> = 3.4%   max ovlp = 4.1%

The E-site overlap is much reduced.

Use your text editor to merge the site positions in `site2.mapb` with
the carefully edited first line of `site.mapb`.  It is now ready
for electronic structure calculations.

Your new site file should read

    % site-data vn=3.0 xpos io=62 nbas=12+{les?11:0} alat=11.18718384 plat= 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
    #                        pos                                   vel                                    eula                   vshft  PL rlx
     C        -0.1037268   0.0000000  -0.0114774    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     N         0.1421229   0.0000000   0.0326819    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H        -0.1330424   0.0000000  -0.1948035    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H        -0.1776707   0.1523111   0.0648669    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H        -0.1776707  -0.1523112   0.0648669    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H         0.2211880   0.1427168  -0.0341978    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H         0.2211880  -0.1427169  -0.0341978    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     H         0.1804493   0.0000000   0.2053792    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     Pb        0.4753318   0.5000000   0.4779234    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     Br        0.4287771   0.5000000  -0.0273080    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     Br        0.4334819   0.0000000  -0.4871517    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     Br       -0.0324470   0.5000000   0.4392057    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.5194899   0.0000000  -0.0172236    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.0088678   0.0000001   0.4888519    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.0824243   0.5000000  -0.0328610    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.2801792  -0.2253603  -0.2427086    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.2801792   0.2253602  -0.2427086    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.1293571  -0.2578272  -0.2491955    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.1293573   0.2578272  -0.2491952    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.2175225  -0.2484617   0.2248361    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E         0.2175225   0.2484617   0.2248361    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.2682113  -0.1900164   0.2974811    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
     E        -0.2682114   0.1900164   0.2974811    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000 0.000000  0 111
