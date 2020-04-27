### Contents

clean.sh		|	removes all output files
fit.py			|	accepts an element name as argument, reads scflog files and fits EOS
init2ctrl.py		|	creates a ctrl file by running blm and lmfa
init.*			|	the structures from the Delta Codes project
kpts.lmto.summary	|	some fairly well converged k-point setups (but now quite old, perhaps not optimal)
README.txt		|
run.sh			|	bash queue script that executes the calculations

### Instructions

All the settings and options are contained in init2ctrl.py -- this is a bit of a mess because some of the options are manipulated in the main loop at the bottom of the file, while others are simply set and printed as-is, which happens towards the top of the file.  The program has an internal representation of the ctrl file as a series of sections and subsections which is used to translate simple variables (the k-mesh, lmxa, etc) into correct tokens for the ctrl file.

You need to have questaal and python in your path (with scipy and numpy): nothing in the scripts affects the path.

You can test the setup and inspect the created ctrl files, eg:

```
python init2ctrl.sh ag 
```

The run.sh file is a bash script that loops over all the elements, runs the 7 volumes and calls the fitting code.  This script can be submitted to a queuing system -- you should edit the header to suit.  If you use archer, then you should also add "aprun" in the indicated place at the top of init2ctrl.py

### Basis setup

I have been using the **extremely simple basis setup scheme** where

RSM1 = r( V = V_rmt*2 ) / sqrt(l+1)
RSM2 = r( V = V_rmt/2 ) / sqrt(l+1)

rmt is set by touching spheres (not even limited to 3.3bohr).  EH is not specified: I use -0.5Ryd which "seems reasonable".

### Discussion

For this basis setup, it actually becomes clear that our results would be in closer agreement with Wien2k if we used LFOCA=0 (I have used the usual LFOCA=1 which is our usual working method).  Clearly we use so many semi-core local orbitals that the core leaking becomes tiny, which is not usually the case, so the difference is probably indeed from the core relaxation.

For the molecular cases, one could try to creep PWEMAX upwards to some kind of convergence; I didn't attempt this because that is also not our usual method -- the 2.0Ryd illustrates a useful effect, the poorer delta values there are reasonable for LMTO.  At the same time, it would be possible to choose better basis parameters for the molecular systems, but that is not easy to automate.

### Bugs

My condition for high-lying local orbitals is currently p(l=2)>0.51, which I thought was completely ok, except that it also identifies the lower group I,II elements, so they get a high-lying d LO, too.  Actually it helps, so they can keep it.  Unfortunately, this algorithm also adds a d HLLO for chlorine, which doesn't do anything and looks a little stupid so that calculation needs to be run elsewere with the HLLO commented out in init2ctrl.py.

The calculation for He is actually pretty terrible (but they are so soft that the Delta value is essentially meaningless) -- just to get a value anyway, I increased the fitting tolerance that I use in fit.py.
