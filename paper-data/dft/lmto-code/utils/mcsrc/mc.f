C*Program mc is a stack-based calculator for matrix operations.
C It uses reverse Polish notation.
C
C*Program mc reads two-dimensional arrays from disk files and does
C various operations on them.  These include several unary operations
C (scaling, inversion, transpose, extracting a submatrix, and some
C others) and binary operations (e.g. addition, subtraction,
C multiplication,concatenation).
C
C*While the matrix operations are performed on the top array(s) on the
C stack, arrays can be assigned names and referred later.
C
C*mc reads its instructions from the command line.  Each command-line
C argument is either a reference to an array (a disk file containing an
C array or an array name), or it is some operator or other switch.
C (Operators and other switches begin with '-'; the '-' distinguishes
C them from references to arrays.)  Thus,
C   mc mat1 mat2 -+
C loads array mat1 onto the stack, then array mat2 onto the stack, and
C then adds the two.  The top stack element contains mat1+mat2, and
C mat1 and mat2 are discarded from the stack.
C In this case mat1 and mat2 are disk files containing two-dimensional
C arrays in ascii format (there is some flexibility in format of the
C disk file; see below); for this example, mat1 and mat2 must have
C the same dimensions.
C
C The unary and binary operators are described below.
C
C*Normally, after all the command-line arguments have been parsed, mc
C writes to stdout out an ascii representation of the top array on the
C stack and exits.  In this case, mc will print out in some default
C format the sum of mat1 and mat2.
C
C*switch '-a aname' copies the top array on the
C stack to some other memory location, assigns name "aname" to it,
C and discards the top array on the stack.  An array reference is made
C mc looks to see wether the name corresponds to a name it knows about;
C if so, it pushes that array onto the top of the stack.  Otherwise,
C it looks for a disk file corresponding to the references, and reads
C the file from disk.  Thus the following line, mc will:
C   mc mat1 -a first mat2 -a second first second -+
C  (1) read mat1 from disk and (2) assign it to 'first'
C  (3) read mat2 from disk and (4) assign it to 'second'
C  (... at this point, the stack is empty, but now mc will:)
C  (5) copy 'first' to the top of the stack
C  (6) copy 'second' to the top of the stack
C  (7) replace the top two elements on the stack with their sum.
C
C*One, and only one, disk file can refer to stdin.  '.' refers to stdin.
C Thus, these two lines should produce the same result:
C   mc mat1 -a first mat2 -a second first second -+
C   cat mat1 | mc . -a first mat2 -a second first second -+
C
C*Structure of disk file.  These files can be ascii or binary.
C In either case, the attributes of the array are the number of rows nr
C and columns nc.  These can be specified in one of several ways.
C The first line may (but need not) contain syntax for prescribing array
C characteristics such as nr and nc.  It uses the syntax for nr, nc:
C % rows (expression) cols (expression)
C Either nr or nc, or both may be specified. If nc is not specified, it
C counts the number of expressions in the the first line read it the
C array, and uses that number for nc.  If nr is not specified, it is
C determined by (the total number of expressions read in) / nc
C (NB: lines of ascii files are read in through a preprocessor that
C  offers much versatility that permits definition of variables,
C  conditional (if-then-else) reading of lines of the file, substitution
C  of algebraic expressions, 'do loop' constructs within the file, and
c  a few others. See routine rdfiln for documentation.)
C
C mc uses a stack-based architecture, though it can resemble an algebraic
C calculator.  It does this by accumulating an operator stack (actually
C two stacks, one for unary and one for binary) that holds operations
C until arrays are available to operate on.  Thus, the following three
C commands equivalently add arrays mat1 and mat2, since the binop '+' is
C held until two arrays are loaded onto the stack.
C   mc -+ mat1 mat2
C and
C   mc mat1 -+ mat2
C and
C   mc mat1 mat2 -+
C
C*List of unary operations and commands to invoke operation:
C (operates on array at top of stack)
C  scale:  scale array by a number
C    -s#[,#]:    scale matrix by #
C                Optional 2nd # => scale is complex
C    -snam(#1,#2)
C                scales matrix by an an element out of a named array,
C                e.g. -sa(2,3)
C
C  sub:    extract submatrix (-sub nr1,nr2,nc1,nc2)
C    -sub nr1,nr2,nc1,nc2
C                replaces top matrix on stack with subblock
C
C  subx:   purges subblock from matrix
C    -subx nr1,nr2,nc1,nc2
C                purges subblock nr1,nr2,nc1,nc2
C
C  expr:   Replace each row of matrix with a sequence of expressions
C          involving the original matrix.
C    -e# expr1 [... expr2]
C                puts "expr1 .. expr2" in columns 1..2..
C                matrix replaced by # columns of expr .   Example:
C                    -e2 x1+x2 x1-x2
C                replaces array by a two-column array whose first column
C                becomes col1+col2 of original array, and whose
C                second column becomes col10col2 of original array.
C
C  symm:   Symmetrize matrix (for square subblock only)
C    -sy         renders matrix hermitian (nr=nc block only)
C
C  include:include rows which satisfy a specified constraints, i.e.
C          exclude rows which does not satisfy the specified constraint.
C    -inc e xpr  include rows that satisfy expr<>0
C                Example: -inc 'x1>1' purges the array of any rows not
C                satisfying col1 larger than 1.
C
C  inv:    array is replaced by its inverse.
C    -i or iq    inverts matrix
C
C  LU decomposition
C    -lu         returns L for hermitian s0 where s0 = LL+
C    Example:  mc h -s1000 -lu -p  -cc -t -x h -s1000 -- -px
C
C  transpose: array is replaced by its transpose.
C    -t          transposes matrix
C
C  push:   push array onto stack
C    -p          pushes top array onto stack
C    -p+n        pushes nth array onto stack
C    -p-n        pushes nth array from bottom onto stack
C  pop:    pop array stack
C    -pop        pops stack
C
C  sort:   sorts rows in matrix
C    -sort expr  sorts rows matrix according to expr
C
C  evl:    replaces matrix by its eigenvalues
C    -evl        replaces matrix by its eigenvalues
C  evec:   replaces matrix by its eigenvectors
C    -evc        replaces matrix by its eigenvectors
C
C  v2dia:  replace matrix vector taken from diagonal, or vice-versa
C    -v2dia      expands vector into diagonal matrix, or contracts
C                diagonal elements of a matrix into a vector
C
C  assign:
C    -a[:nr=#|:nc=#] nam assign top array to nam
C  psh-assign:
C    -ap nam     same as -a, but keep top array on stack
C
C  csum:
C    -csum       sums columns
C  rsum:
C    -rsum       sums rows
C
C  cc:
C    -cc         take complex conjugate of matrix
C
C  cust:
C
C  intrp:
C    -intrp[:opts] list
C                interpolates col 2 at list of points
C                opts can be one of: [rat],[nord=#]
C                See mkdlst for list syntax.
C  diff:
C    -diff[:opts]
C                differentiates cols 2..nc wrt col1
C                opts can be one of: [rat],[nord=#]
C  int:
C    -int[:opts] x0 list
C                integrates col 2 from x0 to a list of upper-bounds
C                opts can be one of: [rat],[nord=#]
C                See mkdlst for list syntax.
C
C  real:
C    -real       take the real part of a complex matrix
C
C  split:
C    -split nam row-list column-list
C                splits array into a set of subblocks with names
C                nam11, nam12, ..., nam21, nam22, ...
C                You can use '.' for row-list or column-list, in which
C                caes it does not divide by rows (columns) but takes all of them.
C                Extract a (1:18,1:4) subarray and create four 9x2 arrays x11, x12, x21, x22:
C                mc -f9f12.6 filnam -sub 1,18,1,4 -split x 1:nr+1:9 1:4:2 x11 -show
C                Assign array to x11 and subtract it from the stack:
C                mc -f9f12.6 filnam -split x . . x11 -- -px
C
C  1: :
C    -1:n        pushes onto stack the unit matrix, dimension n
C
C  array[#1,#2] #,#,...
C                push onto stack new array of rank (#1,#2), populated
C                by #,#,... ,  which must contain #1*#2 elements.
C  rot=string :
C                specify the rotation by a sequence of rotations about
C                specified angles.  Each rotation has syntax (x,y,z)angle; thus
C                -rot='(0,0,1)pi/4,(0,1,0)pi/3,(0,0,1)pi/2' generates a rotation
C                corresponding to Euler angles alpha=pi/4 beta=pi/3 gamma=pi/2
C                Rotation is "active" in that r4z rotates x into y, and y into -x.
C  roteu=string :
C                Same as rot=string, but write out Euler angles corresponding to rotation.
C                Example: the following returns a 1x3 array (1/4,1/3,1/2)
C                mc -roteu=z:pi/4,y:pi/3,z:pi/2  -s1/pi
C  ylm~l=#[~sh|~sh2][~spin|~spin+o|~spin+oonly]
C                Given a rotation matrix on the stack, generate a rotation matrix
C                for real or spherical harmonics, for spinor rotations, or a combination.
C                Example: the matrix rotating real harmonics around z by pi/2
C                (x,y) are mapped into (y,-x), x^2-y^2 is mapped into its negative, y^2-x^2
C                mc -f9f12.6 -rot=z:pi/2 -p -ylm~l=2
C                Example: the matrix rotating spin 1/2 spinor around y by pi/2
C                mc -f9f12.6 -rot=y:pi/2 -ylm~l=0~sh~spin
C                Example: the (22) spinor block of the orbital-only part of a rotated s.h.
C                compared to the rotated s.h.
C                mc -f9f12.6 -rot=y:pi/4 -p -ylm~l=2~sh -a y \
C                -ylm~l=2~sh~spin+oonly -split a 1,nr/2+1,nr+1 1,nr/2+1,nr+1 a22 y -- -px
C
C  pwr=
C    -pwr=#      raises (col-1) matrix to a power
C
C  herm:
C    -herm       renders matrix hermitian (nr=nc block only)
C
C  rowl: create array from list of rows taken from top array
C    -rowl list  create array from list of rows of top array
C  coll: create array from list of columns taken from top array
C    -coll list  create array from list of cols of top array
C
C  rowl:mirr     reverse order of rows
C  coll:mirr     reverse order of columns
C
C  smo:
C    -smo
C
C  unx[:i] #1,#2 columns uncrossed, optionally starting at point i
C    -unx
C
C  unx2:
C   -unx2: cross formerly uncrossed lines at point of closest approach
C
C  at:
C  -at[:i] expr val:  find adjacent rows that bracket expr=val;
C                     return col with x1 interpolated betw brackets.
C
C  abs:  takes the absolute value of array, element by element
C    -abs
C
C  max:
C    -max[:g]    puts the largest element into the first column
C                Optional g returns max value of the entire array
C
C  rep: makes supermatrix out of replicas of smaller matrix
C    -rep:n1,n2
C                concatenates replicas of matrix to create (nr*n1,nc*n2) array
C
C  vassign:
C    -av[:ir,ic] nam:
C                assigns scalar variable 'nam' to element from top array on stack
C
C  shft: shift array by a constant
C    -shft=#[,#2]
C                add constant # to matrix.  Optional #2 => complex constant
C
C  roll: cyclic shift of rows and columns
C    -roll:#1[#2]
C                rows are shifted cyclically by #1 elements;
C                columns are shifted cyclically by #2 elements
C
C  nint[:#1]     For each element in s0, finds nearest integer.
C                If optional arg #1 is present, array is scaled by #1
C                before integer is taken, and scaled by 1/#1 afterwards
C                Thus -nint:2 will find nearest half-integers
C
C*List of binary operations and commands to invoke operation:
C (operates on top two arrays s0 s1)
C  cmpf:  compares numerical contents of two files
C   -cmpf[~ln=#[,#][~col=#[,#][~sep=c][~quiet][~verb][~long][~char|~tol=#][incl][excl][~nchar|~nword|~ndiff|~max]~fn1=f1~fn2=f2
C  ccat:  concatenates columns
C   -ccat        concatenate columns of s0 to s1
C                (the number of rows must be equal)
C  rcat:  concatenates rows
C   -rcat        concatenate rows of s0 to s1
C                (the number of columns must be equal)
C
C  +:  replaces top arrays on stack with their sum
C   -+           adds s1 and s0
C  -:  replaces top arrays on stack with their difference
C   --           adds s1 and -s0
C  x:  replaces top arrays on stack with their product
C   -x           multiplies s1*s0
C
C  tog: toggles top two arrays on stack
C   -tog         toggles s0,s1
C
C  gevl: replaces top arrays on stack with eigenvalues of
C        the generalized eigenvalue problem |s0 - evl s1|=0
C   -gevl        like -evl, but s1 nonorthogality.
C  gevc: replaces top arrays on stack with eigenvectors of
C        the generalized eigenvalue problem |s0 - evl s1|=0
C   -gevc        like -evc, but s1 nonorthogality.
C
C  xe: replaces top arrays on stack with element-by-element product
C   -xe          multiplies s1*s0 element by element
C  de: replaces top arrays on stack with element-by-element ratio
C   -de          divides s1/s0 element by element
C
C  cross:
C   -cross takes the cross product
C  x3:  multiplies top arrays on stack as 3D arrays.
C   -x3          multiplies s1 and s0 as 3D arrays:
C                s1=s1(n11,n21,n31); s0=s0(n10,n20,n30)
C                where n10=nr(0)/n20,n20=nc(1),n30=nc(0) and
C                      n11=n10,n21=nr(1)/n11,n31=n20
C
C  orthos:
C    -orthos     replaces s0 with s1^-1/2 s0 s1^-1/2
C
C*Other switches:
C   -nc=#:  specifies the next array loaded onto the stack must
C           have # columns.
C   -nr=#:  specifies the next array loaded onto the stack must
C           have # rows.
C   -vvar=# define a scalar variable 'var', and assign its value.
C           These may be used in any subsequent algebraic manipulation
C           including the reading of arrays which contain expressions
C           involving variables.
C   -show   shows the current stacks (data and operator).  If -show
C           is the last argument, printing of the array is suppressed.
C   -w fnam writes the top array to a file.  If -w is the last argument,
C           printing of the top array is suppressed.
C
C*Examples: file dat contains
C % cols 2
C pi pi+1 1 2
C
C and file dat2 contains
C % repeat i=1:4
C 1/{i+1} 1/{i+2} 1/{i+3} 1/{i+4}
C % end
C
C Invoking the following commands yield:
C
C
C mc dat         % rows 2 cols 2
C                   3.141593    4.141593
C                   1.000000    2.000000
C
C mc dat dat2    % rows 4 cols 4
C (prints top        0.500000    0.333333    0.250000    0.200000
C of stack)          0.333333    0.250000    0.200000    0.166667
C                    0.250000    0.200000    0.166667    0.142857
C                    0.200000    0.166667    0.142857    0.125000
C
C
C mc -sub 2,3,2,3 dat2 -+ dat     % rows 2 cols 2
C                                     3.391593    4.341593
C                                     1.200000    2.166667
C mc -sub 1,2,1-2,2 dat -sub 1,4,1,4
C puts dat2 into the upper right hand corner of a 4-by-4 array:
C                    0.000000    0.000000    3.141593    4.141593
C                    0.000000    0.000000    1.000000    2.000000
C                    0.000000    0.000000    0.000000    0.000000
C                    0.000000    0.000000    0.000000    0.000000
C
C To generate the norm of a collection of vectors of length 3 (file dat)
C   set norm = "-p -xe -csum -e1 sqrt(x1)"
C   set normalize = "-p $norm -ccat -e3 x1/x4 x2/x4 x3/x4"
C   mc dat $normalize
C For any hermitian matrix h, show h z = z e where z are evecs of h
C   mc h -p -p -evl -tog -evc -tog -v2dia -x -tog -p -evc -x --
C or equivalently,
C   mc h -a h h -evc -p -a z h -evl -v2dia -x h z -x --
C or explicitly for the 2nd eval, evec of h:
C   mc h -a h h -evl -a e h -evc -coll 2 -a z h z -x z -s'e(2,1)' -- -px
C Ditto for generalized eval problem with overlap s
C   mc h -p -p s -tog -gevl -tog s -tog -gevc -tog -v2dia\
C   -x s -tog -x -tog -p s -tog -gevc -x --
C or equivalently,
C   mc h s -p-2 -p-1 -gevc -p  -p-2 -p-1 -gevl -v2dia\
C   -x -p-2 -tog -x -p-1 -p-3 -x --
C or equivalently, one of
C   mc s -ap s h -ap h -gevc -a z s z s h -gevl -v2dia -x -x h z -x --
C   mc s -ap s h -ap h -gevc -a z s h -gevl -v2dia -a e  s z e -x -x  h z -x --
C   mc s -s100 -ap s h -s100 -ap h -gevc -a z s h -gevl -v2dia -a e  s z e -x -x  h z -x --
C Note that the same can be accomplished with -orthos:
C   mc s h -gevl s h -orthos -evl --
C Ex: multiplication of 3D arrays a3 and b3; show cplx version like real,
C     but scaled by a phase:
C   mc a3 -a:nc=2 a3 a3 b -x3 -s0,-1/2 a -real -a:nc=2 a a b -real -x3 --
C Demonstrate that two cyclic rolls (rows, cols) matches double roll:
C   mc dat.mac -a a a -w . -roll:2,0 -roll:0,1 -w . a -roll:2,1 --
C Shear a symmetry line file syml.wzns by tet:
C   mc -vtet=1.01 -vsrtet=tet^.5 -ff3.0,3f12.8,3x,3f12.8 syml.wzns\
C   -e7 x1 'x2*srtet' 'x3*srtet'  'x4/tet'  'x5*srtet'  'x6*srtet'  'x7/tet'\
C   | grep -v rows | sed s'/[.]//'
C Evaluate RMS of list of data, file dat
C   mc dat -a fn fn -s1/nr -rsum -av:1,1 bar \
C   fn -p -xe '-shft=-bar*bar' -rsum -s1/nr -e1 'sqrt(x1)'
C Permute rows, columns of hermitian matrix; show evals unchanged
C   mc h -rowl:pf=perm -coll:pf=perm -evl  h -evl -- -px:15
C Reverse order of rows and columns; put them back again
C   mc h -rowl 1,2,3 -a h h -rowl:mirr -coll:mirr -rowl:mirr -coll:mirr h --
C Formatting large numbers:
C use eg mc '-f%6;6,6d' big
C Interpolating polynomial
C mc dat.gap-tlp.cpa -intrp:rat,nord=3 0:1:.02 intrp.gap-tlp.cpa --
C Test of -unx2: should generate swap for cols 6,7, rows 11:21
C mc dat.swap -unx2 6,7 dat.swap -- -px
C Test of -at: should generate two points, x=173.325299,202.642464
C mc dat.swap -at 100 x11
C Test of -suba
C mc -f9f12.6 -1:3 -s0 -e1 i -s10 -coll 1,1  -1:2 -s0 -e1 i  -t -rowl 1,1,1 -+ -s100 estat -suba:.01 3,5,2,3 estat --
C Demonstrate -array
C mc -f18f8.2 -array'[2,3]' 6,5,4,3,2,1
C Test of looping construct: 6x6 diagonal array with elements 1,1,2,2,3,3
C mc -f18f8.2 -1:2 -a one [ ix=1:3 [ jx=1:3 one -six '?jx<>ix' -s0 '?jx>1' -ccat ] -show ] -rcat -rcat
C
C*Cast of arrays:
C Each array has an associated cast, which describes some information about
C its structure:
C    1s digit = 0 for an integer matrix (for now, not used)
C               1 for a real (double) matrix
C               2 for a complex (double) matrix
C   10s digit = 1 if matrix is hermitian
C  100s digit = 0 for normal (dense) matrix storage
C               1 for sparse matrix, compressed column storage, ie nonzero
C                 elements are grouped by column.  In the Harwell-Boeing format
C                 two integer vectors cp and ri describe their locations:
C                 cp(i) points to first entry in the compress matrix column i;
C                 thus values va(cp(i)) ... va(cp(i+1)-1), in array va
C                 contain the nonzero elements in a(*,i).
C                 ri(k) contains the row index for the kth nonzero element.
C               2 for sparse matrix row format, transposed to Harwell-Boeing.
Cl Local variables
Cl  icast  : 1s digit is cast
Cl         : 0 integer
Cl           1 double precision
Cl           2 double complex with imaginary following real
Cl           3 double complex
Cl           4 double complex with imaginary following real in columns
C*Bugs:
C -nr= and -nc= are not treated as unops.  They are passed to rdm
C and override devault values nr and nc take when a new file is read
C from disk (see rdm).
C*Updates:
C  09 Apr 19 v1.075 New -lu
C  06 Jul 18 v1.072 New -cmpf
C  03 Jul 18 v1.071 New -rot, -roteu, -ylm
C  01 Jul 18 v1.07  New -rot, -roteu, -ylm
C  30 May 18 v1.062 -split can take '.' as row or column list
C  16 May 18 v1.061 new -array and -sarray commands
C  10 Sep 17 v1.059 extension of -rsum:list; replace suba with subc; new suba
C  10 Sep 17 v1.058 new loop [..]
C  14 Apr 12 v1.057 new rowl:mirr
C  14 Apr 12 v1.056 new rowl:pf=
C  10 Mar 11 v1.055 new -nint
C  10 Mar 11 v1.054 improved -cross
C  25 Nov 09 v1.052 new -unx2 and -at
C  17 Jan 06 v1.051  -sub* can accept only 2 limits; (nr,nc) implied as others
C  17 Jun 05 v1.050: New outer product -xo
C   3 Jun 05 v1.049: New -rccat
C  14 Dec 04 v1.048: New -index, new coll:e, new max:i
C  25 Aug 04 v1.047: New -roll:
C   6 Jan 04 v1.046: New -orthos
C  11 Feb 03 v1.045: 3D multiply
C  13 Dec 02 v1.044: redid integration scheme; macro definitions
C  05 Jun 01 v1.043: -shft adds constant shift to top matrix on stack
C  04 Mar 01 v1.042: -int integrates every column, not just second
C  01 Jan 01 v1.041: -s'name(#1,#2)' : valid if 'name' is array name:
C                    scales top array on stack by (#1,#2)-th element
C                    of named array.
C  24 Dec 99 v1.040: New -avscalar-name
C  10 Nov 99 v1.039: -csum,-rsum have optional arguments
C  23 Jun 99 v1.038: mc can read in sparse matrix format.
C                    No operations implemented for now, except -px, -pxc
C  13 Apr 99 v1.037: Added -rep:nrepr,nrepc
C  24 Feb 99 v1.036: Added -r:s#, -max, -abs, -px:#
C  11 Dec 97 v1.035: Added -unx
C  11 Dec 97 v1.034: Enhanced -sub
C  11 Dec 97 v1.033: mc reads environment FMT for default format.
C  11 Dec 97 v1.032: Adapted to read multiple arrays from one file
C  11 Dec 97 v1.031: Adapted to read c*16 binary arrays
C  26 Nov 97 v1.030: Added -smo
C  25 Nov 97 v1.029: Added -tp, -rowl, -coll
C  14 Aug 97 v1.028: mesh option -int
C   4 Jun 97 v1.027  fixed and modified -w
C   4 Jun 97 v1.026  added -intrp
C  29 Jan 97 v1.025  added -a:nr=
C  29 Jan 97 v1.025  added -a:nr=
C  12 Sep 96 v1.024  added -herm
C  22 Apr 96 v1.023  formatting large-file output
C  22 Apr 96 v1.022  raise each element of a one-column matrix to a power
C  22 Mar 96 v1.021  add -1
C  22 Mar 96 v1.020  add complex matrix concatenation
C   7 Mar 96 v1.019  -rphi printout
C   5 Mar 96 v1.018  added inverse of v2dia
C  19 Feb 96 v1.017  added -split
C   6 Dec 95 v1.016  added labels (binary files only for now...)
C   6 Dec 95 v1.015  added -cross
C   9 Nov 95 v1.014  added -real
C  13 Sep 95 v1.013  -rsum and -csum work for complex arrays
C  27 Aug 95 v1.012: -sort switch is now -sort expr
C   8 Jun 95 v1.011: first pass at differentiation of columns, intgx
C  24 May 95 v1.010: load variables table with nr and nc
C   4 May 95 v1.009: -sub ok for complex arrays and extended boundaries
C  28 Apr 95 v1.008: implement trns for complex arrays
C   2 Mar 95 v1.007: -de implemented
C   1 Feb 95 v1.006: -s#,#, -xe,-cc implemented
C  11 Jan 95 v1.005: -br:nr,nc implemented
C   3 Jan 95 v1.004: eigenvalues of real nonsymmetric matrices
C  29 Nov 94 v1.003: first cut at naming arrays.
C  28 Nov 94 v1.002: first cut at complex arrays
C  23 Nov 94 v1.001.
C   8 Nov 94 First attempt.
C ----------------------------------------------------------------
      subroutine fmain
      implicit none
      real(8), parameter :: vsn=1.075d0
      character*(40) first*256,next*256,fmt,rdfmt,sortex,intrps,splitn,splitx,splity,sums,ylms*80
      character*1 f2(20),dc,dc2
      logical lmerge,ltmp,noclos,lqi,ldebug
      integer, parameter :: nsmx=128, mxlev=4, mxexcl=10
C      integer iwk(10000),nlev,iarg,n,ns,iout,
C     .  inam,i,j,i1,i2,k,ifi,jfi,marg,it(4),ib(4),ip,prmx,
C     .  bopstk(nsmx),nbops,uopstk(23),nuops,rdops,nmerge,
C     .  ilsx(nsmx),ilsy(nsmx),nlistx,nlisty,lsubc,lnohead,ms,nprec,nrep(2),nroll(2)
      integer i,i1,i2,iarg,idify,ifi,inam,iout,ip,j,jfi,k,lnohead,lsubc,marg,ms,n,
     .  nbops,nlev,nlistx,nlisty,nmerge,nprec,ns,nuops,prmx,rdops
      integer bopstk(nsmx),ib(4),ilsx(nsmx),ilsy(nsmx),it(4),iwk(10000),nrep(2),
     .  nroll(2),uopstk(23)
      integer irpt(0:mxlev),nrptarg(0:mxlev),lastskip(0:mxlev),iarg0(mxlev)
      character(len=8) :: rptvar(mxlev)
      character(len=80) :: excl(mxexcl),incl(mxexcl),trmstr
      integer vassgn(2)
      integer oiwk,owk,sz,nnz
      integer os(nsmx),nr(nsmx),nc(nsmx),osprs(2,nsmx),icast(nsmx)
      integer ncout,nclip(4)
      integer n1a,n2a,n3a,n1b,n2b,n3b
      integer iseq,nnseq,md,npseq,nnum,is1,is2,ipass,mxnum
      equivalence (f2,first)
      integer,allocatable :: iargslst(:),nseq(:),num(:,:)
      real(8),allocatable :: arr(:,:,:)
      character(1), parameter :: tab = char(9)
      character uopstr(50)*10,bopstr(17)*10,outs*80,ct*20
      double precision scale(2),cshft(2),val(3),x0,pwr,wk(1000),xx,facint,xitrp
      procedure(logical) :: cmdstr,cmdstrsyv,lsequ,a2bin,rdstrn,diffx
      procedure(integer) :: a2vec,fopng,fopnx,garg,getdig,getev,i1mach,iprint,isw,
     .  mkdlst,mkilsd,nargf,parg2,rdm,rdms,wordsw
      procedure(real(8)) :: dval,d1mach
      integer,parameter:: NULLI = -99999
      real(8), parameter :: one=1d0, zer=0d0
C ... parameters for array names
      integer mxnam,namlen
      parameter (mxnam=250,namlen=16)
      character*(namlen) symnam(mxnam),scalnm
      integer symptr(mxnam),nnam,nrn(mxnam),ncn(mxnam),icn(mxnam)
C ... Parameters for expressions
      integer ssvdef,nexpr,mxexpr,wksize,partok
      parameter (ssvdef=256, mxexpr=15,wksize=320 000 000)
      character*(ssvdef) expr(mxexpr),sincl,slist*512,strn*512,strn2*512,filel,rdcsw,doslst
      integer w(wksize)
      common /w/ w
      common /static/ iwk
      data uopstr /
     .  'scale',                !  1
     .  'sub',                  !  2
     .  'subx',                 !  3
     .  'expr',                 !  4
     .  'symm',                 !  5
     .  'include',              !  6
     .  'inv',                  !  7
     .  'transpose',            !  8
     .  'push',                 !  9
     .  'pop',                  ! 10
     .  'sort',                 ! 11
     .  'evl',                  ! 12
     .  'evec',                 ! 13
     .  'v2dia',                ! 14
     .  'assign',               ! 15
     .  'psh-assign',           ! 16
     .  'csum',                 ! 17
     .  'rsum',                 ! 18
     .  'cc',                   ! 19
     .  'cust',                 ! 20
     .  'intrp',                ! 21
     .  'diff',                 ! 22
     .  'int',                  ! 23
     .  'real',                 ! 24
     .  'split',                ! 25
     .  '1:',                   ! 26
     .  'pwr=',                 ! 27
     .  'herm',                 ! 28
     .  'rowl',                 ! 29
     .  'coll',                 ! 30
     .  'smo',                  ! 31
     .  'unx',                  ! 32
     .  'abs',                  ! 33
     .  'max',                  ! 34
     .  'rep',                  ! 35
     .  'vassign',              ! 36
     .  'shft',                 ! 37
     .  'roll',                 ! 38
     .  'unx2',                 ! 39
     .  'at',                   ! 40
     .  'nint',                 ! 41
     .  'rsum',                 ! 42
     .  'csum',                 ! 43
     .  'array',                ! 44
     .  'sarray',               ! 45
     .  'ylm',                  ! 46
     .  'cd',                   ! 47
     .  'cdi',                  ! 48
     .  'ipoly',                ! 49
     .  'lu'/                   ! 50
      data bopstr /'ccat','rcat','+','-','x','tog','gevl','gevc','xe',
     .  'de','cross','x3','orthos','index','rccat','xo','cmpf'/
      sz(i) = mod(icast(i),10)

C#ifdefC DEBUG
C      print *, 'starting mc ...'
C#endif

      call pshpr(0)
      call wkinit(wksize)
      call poppr
C     call wkprnt(1)

      call iinit(iwk,1000)
      ldebug = .false.
      iout = i1mach(2)
      ifi  = 0
      nprec = 8
C ... Number of named arrays
      nnam = 0
C ... Number of arrays on stack
      ns = 0
C ... Number of binary operations on stacks
      nbops = 0
      nuops = 0
      rdops = 10
      call defrr(os(1),1)
      nlev = 0                  ! Nesting level for repeated iteration of CL args
      irpt(0) = -1              ! Flags outside iteration loop
      irpt(1) = -1              ! Flags outside iteration loop
      nrptarg = -1; lastskip = -1 ! Counters for repeated arguments
      iarg0(1) = 0              ! starting point in repeated argument loop
      nmerge = 0                ! > 0 when multiple array merges
      nr(1) = 0
      nc(1) = 0
      sincl = ' '
      lsubc = 0
      lqi = .false.
      prmx  = 0
      noclos = .false.
      nexpr = 0
      ncout = 0
C     matrix constant scaling factor
      scale(1) = 1
      scale(2) = 0
C     matrix constant shift
      cshft(1) = 0
      cshft(2) = 0
      call ivset(nclip,1,4,-1)
      rdfmt = ' '
      goto 10
   19 continue
      outs = ' '
      call awrit0('mc: error parsing '//first//'%a ...',outs,80,iout)
      call cexit(-1,1)
   20 continue
      call mcmsg(vsn,wksize)
      call cexit(0,1)
   10 continue
      fmt = ' '
      call gtenv('FMT',fmt)
      if (fmt == ' ') then
        fmt = '(5f12.6)'
      else
        call word(fmt,1,i1,i2)
        first = '(' // fmt(i1:i2) // ')'
        fmt = first
      endif
      iarg = 1
   15 continue
      first = ' '
      if (.not. cmdstrsyv(iarg,first)) goto 999
C
C      nrec = 0
CC     if (ldebug) rdarg(8:8) = 'p'
CC     variable substitution
C      call rdfiln(0,rdarg,mxlev,loop0,nlin,list,lstsiz,
C     .  ilist,nlist,vnam,ctbl,mxchr,first,recrd,recl,nrec)

      if (ldebug) call awrit0(' mc: parsing '//first,outs,-80,i1mach(2))

      if (first(1:9) == '--version') then
        call awrit1(' mc version %d',' ',80,iout,vsn)
        call cexit(0,1)
      elseif (first(1:5) == '-quit') then
        call cexit(0,1)
      endif
      if (first(1:1) == '[') then
        nlev = nlev+1
        lastskip(nlev) = 0 ! number of arguments to skip on last iteration
        do
          if (first(lastskip(nlev)+2:lastskip(nlev)+2) /= '/') exit
          lastskip(nlev) = lastskip(nlev)+1
        enddo
        iarg = iarg+1
        ip = 1
        if (.not. cmdstrsyv(iarg,first)) goto 999
        rptvar(nlev) = 'ir'
        i = index(first(ip:),'=')
        if (i > 0) then
          rptvar(nlev) = first(ip:i-1)
          ip = i+1
        endif
        nrptarg(nlev) = mkilsd(trim(first(ip:)),-1,j)
        if (nrptarg(nlev) < 0) goto 999
        if (allocated(iargslst)) deallocate(iargslst)
        allocate(iargslst(nrptarg(nlev)))
        nrptarg(nlev) = mkilsd(trim(first(ip:)),nrptarg(nlev),iargslst)
        iarg = iarg+1
        iarg0(nlev) = iarg
        irpt(nlev) = 0          ! Flag that loop is starting
      endif

      if (lastskip(nlev) > 0 .and. irpt(nlev) == nrptarg(nlev)) then
        if (.not. cmdstr(iarg+lastskip(nlev),next)) call rx('mcx: missing "]"')
        if (next(1:1) == ']') then
          irpt(nlev) = -1
          iarg = iarg+lastskip(nlev)+1
          goto 15
        endif
      endif
      if (first(1:1) == ']' .or. irpt(nlev) == 0) then
        irpt(nlev) = irpt(nlev)+1
!        if (first(1:1) /= ']') goto 15  ! First pass
C       if (irpt(nlev) > nrptarg(nlev)) call rx('repeat list exhausted')
C       if (iarg0(nlev) == 0) call rx('mcx: encounted [ without ]')
        if (irpt(nlev) > nrptarg(nlev)) then
          irpt(nlev) = -1
          nlev = nlev-1
          if (nlev < 0) call rx("encountered ']' before '['")
          iarg = iarg+1
          goto 15
        else
          call lodsyv('i',1,dble(irpt(nlev)),ip)
          call lodsyv(trim(rptvar(nlev)),1,dble(iargslst(irpt(nlev))),ip)
          if (ldebug) call shosyv(0,0,0,6)
          iarg = iarg0(nlev)
          goto 15
        endif
      endif
C     Conditional parse next argument
      if (irpt(nlev) > 0 .and. first(1:1) == '?') then
        i = 1
        if (.not. a2bin(first,j,2,0,' ',i,-1)) goto 20
        if (j == 0) iarg = iarg+1
        iarg = iarg+1
        goto 15
      endif

      if (first(1:1) /= '-') goto 30
C ... Load variables table with nr,ns from top matrix
      if (ns > 0) then
        call lodsyv('nr',1,dble(nr(ns)),ip)
        call lodsyv('nc',1,dble(nc(ns)),ip)
      endif
      ip = 0
      if (first(1:2) == '-f') then
        fmt = '('
        call strcat(fmt,1,' ',f2(3),99,' ',n)
        call strcat(fmt,99,' ',')',1,' ',n)
      else if (first(1:2) == '-C') then
      else if (first(1:6) == '-debug') then
        ldebug = .true.
      else if (first(1:6) == '-v2dia') then
        nuops = nuops+1
        uopstk(nuops) = 14
        if (ns > 0) goto 40
C ... Variables declaration
      else if (first(1:2) == '-v') then
        first = first(1:len(first)-1) // ' '
        call wordg(first,11,'a-zA-Z0-9_*/+^-',1,i,j)
        if (first(j+1:j+1) == '=') then
          i = 2
          call parsyv(first,len(first),999,1,i)
        else
          goto 20
        endif
C ... show stack
      else if (first(1:5) == '-show') then
        call awrit5('# %i named arrays, %i on stack; pending %i unops'//
     .    ' %i bops (vsn %d)',' ',80,iout,nnam,ns,nuops,nbops,vsn)
        if (nnam > 0) then
          call awrit0('# array             nr    nc   cast',' ',80,iout)
          do  i = 1, nnam
            outs = ' '
            if (icn(i) == 1)  outs(1:10) = '   real'
            if (icn(i) == 11) outs(1:10) = '   symm'
            if (icn(i) == 2)  outs(1:10) = '   complex'
            if (icn(i) == 12) outs(1:10) = '   herm'
            call awrit2('# '//symnam(i)//'%,4i%,6i'//outs(1:10),' ',80,
     .        iout,nrn(i),ncn(i))
          enddo
        endif
        if (ns > 0) then
          call awrit0('# stack             nr    nc   cast',' ',80,iout)
          do  i = ns, 1, -1
            outs = ' '
            if (icast(i) == 1)  outs(1:10) = '   real'
            if (icast(i) == 11) outs(1:10) = '   symm'
            if (icast(i) == 2)  outs(1:10) = '   complex'
            if (icast(i) == 12) outs(1:10) = '   herm'
            call awrit3('#%,5i%,16i%,6i'//outs(1:10),' ',80,iout,
     .        i,nr(i),nc(i))
          enddo
        endif
        if (nuops > 0) then
          call awrit0('# unnop  operation',' ',80,iout)
          do  13  i = 1, nuops
   13     print '(i4,5x,a)', i, uopstr(mod(uopstk(i),100))
        endif
        if (nbops > 0) then
          call awrit0('# binop  operation',' ',80,iout)
          do  14  i = 1, nbops
   14     print '(i4,5x,a)', i, bopstr(bopstk(i))
        endif
        if (nargf() == iarg+1) call cexit(0,1)
C ... quick-read
      else if (first(1:3) == '-qr') then
        rdops = 10*(rdops/10) + 1
C ... binary-read
      elseif (garg('-br:',iarg,2,', ',2,2,marg,it,iwk) /= 0) then
        if (marg /= 2) goto 19
        rdops = 10*(rdops/10) + 3
        nr(ns+1) = iwk(1)
        nc(ns+1) = iwk(2)
C ... read modifiers
      elseif (first(1:3) == '-r:') then
        rdcsw = first(3:)
      elseif (first(1:3) == '-br') then
        rdops = 10*(rdops/10) + 2
C ... terse print-out
      else if (first(1:3) == '-px') then
        prmx = 1
        i = 4
        if (first(1:4) == '-pxc') then
          prmx = 2
          i = 5
        endif
        if (first(i:i) == ':') then
          if (.not. a2bin(first,nprec,2,0,' ',i,-1)) goto 20
        endif
C ... unop or binop
      elseif (garg('-sub',iarg,1,'xsvca ',4,0,marg,it,nclip) /= 0) then
        lsubc = it(1)
        if (lsubc >= 2 .and. lsubc <= 5) then
          scale = [1d0, 0d0]
          if (first(6:6)==':') then
            j = 6
            j = a2vec(first,len(first),j,4,', ',2,2,2,it,scale)
            if (j <= 0) goto 19
          endif
        endif
        iarg = iarg+1
        i = garg(' ',iarg,2,', ',2,4,marg,it,nclip)
        if (i /= 2 .and. i /= 4) call rx('mc: bad arg to -sub')
C       2 arguments => just (i,j) supplied; upper limit = (nr,nc)
        if (i == 2) then
          nclip(3) = nclip(2)
          nclip(2) = nr(ns)
          nclip(4) = nc(ns)
        endif
        if (lsubc == 4 .or. lsubc == 5) then
          nbops = nbops+1
          bopstk(nbops) = 17
          lsubc = lsubc - 4
        else
          nuops = nuops+1
          uopstk(nuops) = 2
          if (lsubc == 6) lsubc = 4
        endif
        if (ns > 0) goto 40
C ... unops
      else if (first(1:5) == '-sort') then
        nuops = nuops+1
        uopstk(nuops) = 11
        iarg = iarg+1
        if (.not. cmdstrsyv(iarg,sortex)) goto 20
        if (ns > 0) goto 40
      else if (first(1:3) == '-sy') then
        nuops = nuops+1
        uopstk(nuops) = 5
        if (ns > 0) goto 40
      else if (first(1:4) == '-inc') then
        iarg = iarg+1
        if (.not. cmdstrsyv(iarg,sincl)) goto 20
        nuops = nuops+1
        uopstk(nuops) = 6
        if (ns > 0) goto 40
      else if (first(1:3) == '-i ' .or. first(1:4) == '-iq ') then
        nuops = nuops+1
        uopstk(nuops) = 7
        lqi = first(1:4) == '-iq '
        if (ns > 0) goto 40
      else if (first(1:3) == '-t ') then
        nuops = nuops+1
        uopstk(nuops) = 8
        if (ns > 0) goto 40
      else if (first(1:3) == '-p ') then
        nuops = nuops+1
        uopstk(nuops) = 9
        if (ns > 0) goto 40
      elseif (garg('-p+',iarg,2,' ',1,1,k,i,n) == 1) then
        nuops = nuops+1
        uopstk(nuops) = 100*n+9
        if (ns > 0) goto 40
      elseif (garg('-p-',iarg,2,' ',1,1,k,i,n) == 1) then
        nuops = nuops+1
        uopstk(nuops) = -100*n-9
        if (ns > 0) goto 40
      else if (first(1:4) == '-pop') then
        nuops = nuops+1
        uopstk(nuops) = 10
        if (ns > 0) goto 40
      else if (first(1:4) == '-evl') then
        nuops = nuops+1
        uopstk(nuops) = 12
        if (ns > 0) goto 40
      else if (first(1:4) == '-evc') then
        nuops = nuops+1
        uopstk(nuops) = 13
        if (ns > 0) goto 40
      else if (first(1:4) == '-wap ') then
        jfi = iout
        i = max(ns,1)
        call apprm(filel,icast(i),jfi,fmt,w(os(i)),nr(i),nc(i))
        if (nargf() == iarg+1) call cexit(0,1)
      else if (lsequ(first,'-w ',3,' ',n) .or.
     .         lsequ(first,'-w:',3,' ',n) .or.
     .         lsequ(first,'-bw:',4,' ',n) .or.
     .         lsequ(first,'-bw ',4,' ',n)) then
        ltmp = lsequ(first,'-bw',3,' ',n)
        i = max(ns,1); lnohead = 0
C   ... Switches
        filel = ' '
        i2 = index(first,'w')+1
        dc = first(i2:i2)
        if (dc /= ' ' .and. wordsw(first,dc,'l','= ',i1) > 0) then
          call nwordg(first,0,dc,1,i1,i2)
          filel = first(i1+1:i2)
        endif
        if (dc /= ' ' .and. wordsw(first,dc,'nohead','',i1) > 0) then
          lnohead = 1  ! Write array without header
        endif

C        if (dc /= ' ') then
C          i1 = i2+1
C   21     call skipbl(first,len(first),i1)
C          i1 = i1+1
CC         print *, 'i1 now',i1,' ',first(i1:i1)
C          if (i1 < len(first)) then
C            if (first(i1:i1+1) == 'l=') then
C              call nwordg(first,0,dc,1,i1,i2)
C              if (i2 < 1) i2 = len(first)
C              filel = first(i1+2:i2)
C              i1 = i2+1
C              goto 21
C            endif
C            call rxs('mc: failed to parse: ',first)
C          endif
C       endif
C   ... Read name of file, or optionally label
        iarg = iarg+1
        if (.not. cmdstrsyv(iarg,first)) goto 999
        if (filel /= ' ') call awrit0('#mc write '''//first//'%a'': "'
     .    //filel//'%a"',strn,len(strn),-iout)
        if (ltmp) then
          jfi = fopng(first,-1,4)
          if (prmx /= 0) then
            print '(a)', '#mc: no sparse binary write'
          else
            call prm(1+10*lnohead,filel,icast(i),jfi,fmt,osprs(1,i),w(os(i)),nr(i),nc(i))
          endif
          call fclr(first,jfi)
        else
          if (first == '.') then
            jfi = iout
          else
            jfi = fopng(first,-1,0)
          endif
          if (prmx /= 0) then
            call xprm(nprec,icast(i),jfi*(3-2*prmx),osprs(2,i),w(os(i)),nr(i),nc(i))
            prmx = 0
          else
            call prm(0+10*lnohead,filel,icast(i),jfi,fmt,
     .        osprs(1,i),w(os(i)),nr(i),nc(i))
          endif
          if (jfi /= iout) call fclr(first,jfi)
        endif
        if (nargf() == iarg+1) call cexit(0,1)
      else if (first(1:4) == '-abs') then
        nuops = nuops+1
        uopstk(nuops) = 33
        if (ns > 0) goto 40
      else if (first(1:5) == '-nint') then
        nuops = nuops+1
        uopstk(nuops) = 41
        i = 6
        facint = 1
        if (first(i:i) == ':') then
          if (.not. a2bin(first,facint,4,0,' ',i,-1))
     .      call rxs('mc : failed to parse arguments to ',first)
        endif
        if (ns > 0) goto 40
      else if (first(1:3) == '-av') then
        vassgn(1) = 1
        vassgn(2) = 1
        i1 = 4
        if (first(4:4) == ':') then
          if (a2vec(first,len(first),i1,2,', ',2,2,2,it,vassgn) < 0)
     .      call rxs('mc : failed to parse arguments to ',first)
          i1 = i1+1
        endif
        if (.not. cmdstrsyv(iarg+1,scalnm)) goto 20
        iarg = iarg+1
        nuops = nuops+1
        uopstk(nuops) = 36
        if (ns > 0) goto 40

      else if (first(1:4)=='-unx' .or. first(1:3)=='-at') then
        nuops = nuops+1
        uopstk(nuops) = 32
        if (first(1:5) == '-unx2') uopstk(nuops) = 39
        if (first(1:3) == '-at') uopstk(nuops) = 40
        iarg = iarg+1
        wk(2) = 1
        wk(3) = 1
        if (uopstk(nuops) == 32 .and. first(5:5) == ':') then
          i = 5
          if (.not. a2bin(first,wk(3),4,0,' ',i,-1)) goto 20
        elseif (uopstk(nuops) == 39 .and. first(6:6) == ':') then
          i = 6
          if (.not. a2bin(first,wk(3),4,0,' ',i,-1)) goto 20
        elseif (uopstk(nuops) == 40 .and. first(4:4) == ':') then
          i = 4
          if (.not. a2bin(first,wk(2),4,0,' ',i,-1)) goto 20
        endif
        if (uopstk(nuops) == 40) then
          if (garg(' ',iarg,4,', ',1,4,i,it,wk) /= 1) goto 19
          iarg = iarg+1
          if (.not. cmdstrsyv(iarg,sincl)) goto 19
        else
          if (garg(' ',iarg,4,', ',2,4,i,it,wk) /= 2) goto 19
        endif
        if (ns > 0) goto 40

      else if (first(1:7) == '-array[' .or. first(1:8) == '-sarray[' .or.
     .         first(1:9) == '-sarrayr[' .or. first(1:9) == '-sarrayc[' .or. first(1:9) == '-sarrayz[') then
        j = index(first,'[')
        j = a2vec(first,len(first),j,2,',]',2,2,2,it,iwk)
        if (j /= 2) call rxs('failed to parse ',trim(first))

C       Allocate new array
        if (first(1:7) == '-array[') then
          ns = ns+1
          if (ns >= nsmx) call rxi('mc: increase nsmx, ns=',ns)
          nr(ns+1) = 0
          nc(ns+1) = 0
          icast(ns) = 1
          iarg = iarg+1
          if (.not. cmdstrsyv(iarg,next))
     .      call rx(trim(first)//' requires next argument with list of elements')
          nr(ns) = iwk(1)
          nc(ns) = iwk(2)
          j = nr(ns)*nc(ns)*sz(ns)
          call defrr(os(ns),-j)
          i = mkdlst(next,-1024*d1mach(3),0,wk)
          ltmp = i == 1 .and. i /= iwk(1)*iwk(2)
          if (ltmp) then
            k = 0
            i = a2vec(next,len_trim(next),k,4,', ',2,3,iwk(1)*iwk(2),wk,w(os(ns)))
          endif
          if (i /= iwk(1)*iwk(2)) call rx1('string after '//
     .      trim(first)//' must consist of %sx%2i elements',iwk)
          if (.not. ltmp) i = mkdlst(next,-1024*d1mach(3),i,w(os(ns)))
          goto 40
        else

          if (ns < iwk(1)*iwk(2)) call rx1('switch '//
     .      trim(first)//' requires %sx%2i arrays on stack',iwk)

          i1 = 0; if (first(8:8) == 'c' .or. first(8:8) == 'z') i1 = 10
          i2 = 0; if (first(8:8) == 'r' .or. first(8:8) == 'z') i2 = 10

          do  i = 1, iwk(1)

C           column-concatenate last iwk(2) arrays on top
            do  j = 1, iwk(2)-1
C             print *, 'i,j,ns before',j,ns,sz(ns),nr(ns),nc(ns)
              call rccatm(i1+1,ns,icast,nr,nc,os)
              call mergesn(1,ns,icast,nr,nc,os)
              icast(ns) = sz(ns)
C             print *, 'i,j,ns after',j,ns,sz(ns),nr(ns),nc(ns)
            enddo

            if (i > 1) then

              ns = ns+1
              icast(ns) = it(1)
              nr(ns)    = it(2)
              nc(ns)    = it(3)
              j = nr(ns)*nc(ns)*sz(ns)
              call defrr(os(ns),-j)
              call dcopy(j,arr,1,w(os(ns)),1)
              deallocate(arr)

C             print *, 'i,j,ns after arr copy ',j,ns,sz(ns),nr(ns),nc(ns),icast(ns)
              call rccatm(i2+2,ns,icast,nr,nc,os)
              call mergesn(1,ns,icast,nr,nc,os)
              icast(ns) = sz(ns)
C             print *, 'i,j,ns after merge',j,ns,sz(ns),nr(ns),nc(ns)

            endif

C           Copy top array to temp arr and pop from stack, for future row merge
            if (i < iwk(1)) then
C           Keep bookeeping data for ns-1 since stack will be destroyed
            it(1) = icast(ns)
            it(2) = nr(ns)
            it(3) = nc(ns)
            j = sz(ns)*nr(ns)*nc(ns)
C           print *, 'i,j,ns before copy to arr',j,ns,sz(ns),nr(ns),nc(ns),icast(ns)
            allocate(arr(j,1,1))
            call dcopy(j,w(os(ns)),1,arr,1)
            call clears(ns,nr,nc,icast)
C           print *, 'i,j,ns after copy to arr',j,ns,sz(ns),nr(ns),nc(ns),icast(ns)
            endif
          enddo

          goto 40

        endif

      elseif (garg('-a',iarg,1,' :p',1,0,i,it,val) == 1) then
C   ... look for options :nr= or :nc=
        if (it(1) == 2) then
          if (first(4:6) == 'nr=' .or. first(4:6) == 'nc=') then
            j = 6
            if (.not. a2bin(first,i,2,0,' ',j,len(first))) goto 19
            if (first(4:6) == 'nr=') then
              j = (nr(ns)*nc(ns))/i
            else
              j = i
              i = (nr(ns)*nc(ns))/j
            endif
            if (i*j /= nr(ns)*nc(ns)) call awrit4('#mc (warning) '//
     .        'size mismatch converting nr,ns=%i,%i to nr,ns=%i,%i',
     .        ' ',80,-iout,nr(ns),nc(ns),i,j)
            nr(ns) = i
            nc(ns) = j
          endif
        endif
        if (nnam>mxnam) call fexit(-1,9,'mc: too many array names',0)
        iarg = iarg+1
        if (.not. cmdstrsyv(iarg,symnam(nnam+1))) goto 20
        nuops = nuops+1
        uopstk(nuops) = 15
C   ... option -ap
        if (it(1) == 3) uopstk(nuops) = 16
        if (ns > 0) goto 40
      else if (first(1:5) == '-csum') then
        nuops = nuops+1
        uopstk(nuops) = 17
        sums = ' '
        if (first(6:6) == ':') sums = first(7:)
        if (ns > 0) goto 40
      else if (first(1:5) == '-rsum') then
        nuops = nuops+1
        uopstk(nuops) = 18
        sums = ' '
        if (first(6:6) /= ' ') sums = first(7:)
        if (ns > 0) goto 40
      else if (first(1:4) == '-ylm') then
        nuops = nuops+1
        uopstk(nuops) = 46
        ylms = first(5:)
        if (ns > 0) goto 40
      else if (first(1:4) == '-cc ') then
        nuops = nuops+1
        uopstk(nuops) = 19
        if (ns > 0) goto 40
      else if (first(1:5) == '-cdi ') then
        nuops = nuops+1
        uopstk(nuops) = 48
        if (ns > 0) goto 40
      else if (first(1:7) == '-ipoly=') then
        nuops = nuops+1
        uopstk(nuops) = 49
        k = 7
        if (a2vec(first,len_trim(first),k,4,', ',2,2,1,it,xitrp) < 1)
     .    call rx('mc failed to parse '//trim(first))
        if (ns > 0) goto 40
      else if (first(1:4) == '-cd ') then
        nuops = nuops+1
        uopstk(nuops) = 47
        if (ns > 0) goto 40
      else if (first(1:5) == '-cust ') then
        nuops = nuops+1
        uopstk(nuops) = 20
        if (ns > 0) goto 40
      else if (first(1:6) == '-split') then
        nuops = nuops+1
        uopstk(nuops) = 25
        if (.not. cmdstrsyv(iarg+1,splitn)) goto 19
        if (.not. cmdstrsyv(iarg+2,splitx)) goto 19
        if (splitx == '-' .or. splitx == '.') then
          splitx = '1,nr+1'
        endif
        if (.not. cmdstrsyv(iarg+3,splity)) goto 19
        if (splity == '-' .or. splity == '.') then
          splity = '1,nc+1'
        endif
        iarg = iarg+3
        if (ns > 0) goto 40
      else if (first(1:2) == '-1') then
        ns = ns+1
        nr(ns) = 1
        if (first(3:3) == ':') then
          j = 3
          if (.not. a2bin(first,nr(ns),2,0,' ',j,len(first))) goto 19
        endif
        nc(ns) = nr(ns)
        if (ns >= nsmx) call rxi('mc: increase nsmx, ns=',ns)
        nr(ns+1) = 0
        nc(ns+1) = 0
        icast(ns) = 11
        j = nr(ns)*nc(ns)*sz(ns)
        call defrr(os(ns),-j)
        call dcopy(nr(ns),one,zer,w(os(ns)),nr(ns)+1)
        if (ns > 0) goto 40
      else if (first(1:6) == '-roteu') then
        if (ns<1 .or. nr(ns) /= 3 .or. nc(ns) /= 3)
     .    call rx('-roteu requires that 3x3 rotation matrix resides on stack')
        call rm2eua(w(os(ns)),wk(1),wk(2),wk(3))
        call dcopy(3,wk,1,w(os(ns)),1)
        nr(ns) = 1
        goto 40
      else if (first(1:5) == '-rot=') then
        ns = ns+1; nr(ns) = 3; nc(ns) = 3; icast(ns) = 1
        j = nr(ns)*nc(ns)*sz(ns)
        call defrr(os(ns),-j)
        call a2rotm(first(6:len_trim(first)),.false.,.false.,w(os(ns)))
        goto 40
      else if (first(1:4) == '-tp ' .or.
     .         first(1:9) == '-rowl:pf=' .or.
     .         first(1:10) == '-rowl:ipf=' .or.
     .         first(1:10) == '-rowl:mirr' .or.
     .         first(1:6) == '-rowl ' .or.
     .         first(1:9) == '-coll:pf=' .or.
     .         first(1:10) == '-coll:ipf=' .or.
     .         first(1:10) == '-coll:mirr' .or.
     .         first(1:6) == '-coll ' .or.
     .         first(1:8) == '-coll:e ') then
        wk(1) = 0
        iarg = iarg+1
        call defrr(owk,1001) ! w(owk) = array size and (owk+1:) = array
        if (first(1:4) == '-tp ') then
          ns = ns+1
          nr(ns) = 1
          icast(ns) = 1
          call defrr(os(ns),300000)
          if (garg(' ',iarg,4,',:~ ',4,1000,i,iwk,w(owk)) <1) goto 19
          call expand(i,iwk,w(owk),w(os(ns)),nr(ns),nc(ns))
          j = nr(ns)*nc(ns)*sz(ns)
          call redfrr(os(ns),j)
        else
          if (first(1:6) == '-rowl ') then
            nuops = nuops+1
            uopstk(nuops) = 29
          elseif (first(1:10) == '-rowl:mirr' .or.
     .            first(1:10) == '-coll:mirr') then
            nuops = nuops+1
            uopstk(nuops) = 29
            k = nr(ns)
            if (first(1:5) == '-coll') then
              k = nc(ns)
            uopstk(nuops) = 30
            endif
            w(owk) = k
            do  i = 1, k
              w(owk+i) = k-i+1
            enddo
            iarg = iarg-1  ! no cmd line argument associated w/ list
            goto 40
          elseif (first(1:9) == '-rowl:pf=' .or.
     .            first(1:10) == '-rowl:ipf=' .or.
     .            first(1:9) == '-coll:pf=' .or.
     .            first(1:10) == '-coll:ipf=') then
            nuops = nuops+1
            uopstk(nuops) = 29
            if (first(1:5) == '-coll') uopstk(nuops) = 30
            if (ns == 0) call rx('mc: no array to permute')
            call rdperm(first,nr(ns),nc(ns),1000,w(owk),w(owk+1))
            iarg = iarg-1  ! no cmd line argument associated w/ list
            goto 40
          else
            if (first(1:8) == '-coll:e ') wk(1) = 1  ! flag :e switch
            nuops = nuops+1
            uopstk(nuops) = 30
          endif
          if (.not. cmdstrsyv(iarg,first)) goto 19
          call mkilst(first,w(owk),w(owk+1))  ! 1st elt holds #, 1: index table
          if (w(owk) < 0) goto 19
        endif
        if (ns > 0) goto 40
      else if (first(1:6) == '-intrp') then
        intrps = ' '
        if (first(7:7) == ':') intrps = first(8:)
        iarg = iarg+1
        if (.not. cmdstrsyv(iarg,slist)) goto 20
        nuops = nuops+1
        uopstk(nuops) = 21
        if (ns > 0) goto 40
      else if (first(1:5) == '-diff') then
        intrps = ' '
        if (first(6:6) == ':') intrps = first(7:)
        nuops = nuops+1
        uopstk(nuops) = 22
        if (ns > 0) goto 40
      else if (first(1:5) == '-smo ') then
        nuops = nuops+1
        uopstk(nuops) = 31
        iarg = iarg+1
        if (garg(' ',iarg,4,', ',2,4,i,it,wk) /= 4) goto 19
        if (ns > 0) goto 40

C     Shift top array by a constant amount
      else if (first(1:6) == '-shft=') then
        i = garg('-shft=',iarg,4,', ',2,2,marg,it,cshft)
        if (i < 0) goto 20
        nuops = nuops+1
        uopstk(nuops) = 37
        if (ns > 0) goto 40

C     Roll array
      elseif (garg('-roll:',iarg,2,', ',2,2,marg,it,nroll) /= 0) then
        if (marg == -1) goto 19
        if (marg == 1) nroll(2) = 0
        nuops = nuops+1
        uopstk(nuops) = 38
        if (ns > 0) goto 40

C     Scale top array by a constant factor
      else if (first(1:2) == '-s') then
C       i=1,2 if a valid scalar expression, or pair of expressions
        i = garg('-s',iarg,4,', ',2,2,marg,it,scale)
C       Look for syntax string(..) where string = array name
C       Block exits with i<0 if invalid syntax
        if (i < 0) then
          call wordg(first,1,'() ',1,i,j)
          if (first(j+1:j+1) == '(') then
            first(j+1:j+1) = ' '
            call tokmat(first(3:),symnam,nnam,namlen,' ',i,k,.false.)
          endif
          if (i >= 0) then
            iwk(1) = 1
            iwk(2) = 1
            if (first(j+2:j+2) == ')') then
            else
              k = j
              k = a2vec(first,len(first),k,2,',)',2,2,2,it,iwk)
              if (k < 0) i = -1
            endif
          endif
          if (i < 0) then
            call rxs('mc: failed to parse scaling factor: ',first)
          endif
          first(j+1:j+1) = '('
          i = i+1
          k = 4
          if (mod(icn(i),10) == 2) k = 0
          if (iwk(1) > nrn(i) .or. iwk(2) > ncn(i))
     .      call rxs('mc: array element out of range: ',first(3:))
          call dpzero(scale,2)
          call ymscop(k,1,1,nrn(i),1,iwk(1)-1,iwk(2)-1,0,0,
     .      w(symptr(i)),nrn(i)*ncn(i),scale,1)
C          if (ncn(i)*nrn(i) /= 1) call awrit0(
C     .      '#-s (warning):  using (1,1) argument of matrix '''//
C     .      first(3:)//'%a''',' ',80,iout)
        endif
        nuops = nuops+1
        uopstk(nuops) = 1
        if (ns > 0) goto 40
      elseif (first(1:4) == '-max') then
        nuops = nuops+1
        uopstk(nuops) = 34
        wk(1) = 0
        if (first(1:6) == '-max:g') wk(1) = 1
        if (first(1:6) == '-max:i') wk(1) = 2
        if (ns > 0) goto 40
      elseif (garg('-rep:',iarg,2,', ',2,2,marg,it,nrep) /= 0) then
        if (marg == -1) goto 19
        if (marg == 1) nrep(2) = 1
        nuops = nuops+1
        uopstk(nuops) = 35
        if (ns > 0) goto 40
      else if (first(1:4) == '-int') then
        intrps = ' '
        if (first(5:5) == ':') intrps = first(6:)
        if (garg(' ',iarg+1,4,' ',1,1,marg,it,x0) < 1) goto 20
        iarg = iarg+2
        if (.not. cmdstrsyv(iarg,slist)) goto 20
        nuops = nuops+1
        uopstk(nuops) = 23
        if (ns > 0) goto 40
      else if (first(1:5) == '-real') then
        nuops = nuops+1
        uopstk(nuops) = 24
        if (ns > 0) goto 40
      elseif (garg('-pwr=',iarg,4,' ',1,1,k,i,pwr) == 1) then
        nuops = nuops+1
        uopstk(nuops) = 27
        if (ns > 0) goto 40
      else if (first(1:5) == '-herm') then
        nuops = nuops+1
        uopstk(nuops) = 28
        if (ns > 0) goto 40
      else if (first(1:3) == '-lu') then
        nuops = nuops+1
        uopstk(nuops) = 50
        if (ns > 0) goto 40
      else if (first(1:2) == '-e') then
        j = 2
        if (.not. a2bin(first,i,2,0,' ',j,len(first))) i = 1
   27   nexpr = nexpr+1
        i = i-1
        if (nexpr > mxexpr) stop 'too many expressions'
        iarg = iarg+1
        if (.not. cmdstrsyv(iarg,expr(nexpr))) goto 20
        if (i > 0) goto 27
        nuops = nuops+1
        uopstk(nuops) = 4
        if (ns > 0) goto 40
      else if (first(1:4) == '-nc=') then
        i = 4
        if (.not. a2bin(first,nc(ns+1),2,0,' ',i,-1)) goto 20
      else if (first(1:4) == '-nr=') then
        i = 4
        if (.not. a2bin(first,nr(ns+1),2,0,' ',i,-1)) goto 20
      else if (first(1:5) == '-ccat') then
        nbops = nbops+1
        bopstk(nbops) = 1
        if (ns > 0) goto 40
      else if (first(1:5) == '-rcat') then
        nbops = nbops+1
        bopstk(nbops) = 2
        if (ns > 0) goto 40
      else if (first(1:3) == '-+ ') then
        nbops = nbops+1
        bopstk(nbops) = 3
        if (ns > 0) goto 40
      else if (first(1:3) == '-- ') then
        nbops = nbops+1
        bopstk(nbops) = 4
        if (ns > 0) goto 40
      else if (first(1:3) == '-x3') then
        nbops = nbops+1
        bopstk(nbops) = 12
        if (ns > 0) goto 40
      else if (first(1:3) == '-x ') then
        nbops = nbops+1
        bopstk(nbops) = 5
        if (ns > 0) goto 40
      else if (first(1:7) == '-cross ') then
        nbops = nbops+1
        bopstk(nbops) = 11
        if (ns > 0) goto 40
      else if (first(1:4) == '-tog') then
        nbops = nbops+1
        bopstk(nbops) = 6
        if (ns > 0) goto 40
      else if (first(1:5) == '-gevl') then
        nbops = nbops+1
        bopstk(nbops) = 7
        if (ns > 0) goto 40
      else if (first(1:5) == '-gevc') then
        nbops = nbops+1
        bopstk(nbops) = 8
        if (ns > 0) goto 40
      else if (first(1:7) == '-orthos') then
        nbops = nbops+1
        bopstk(nbops) = 13
        if (ns > 0) goto 40
      else if (first(1:6) == '-index') then
        nbops = nbops+1
        bopstk(nbops) = 14
        if (ns > 0) goto 40
      else if (first(1:4) == '-xe ') then
        nbops = nbops+1
        bopstk(nbops) = 9
        if (ns > 0) goto 40
      else if (first(1:4) == '-de ') then
        nbops = nbops+1
        bopstk(nbops) = 10
        if (ns > 0) goto 40
      else if (first(1:6) == '-rccat') then
        nbops = nbops+1
        bopstk(nbops) = 15
        if (ns > 0) goto 40
      else if (first(1:3) == '-xo') then
        nbops = nbops+1
        bopstk(nbops) = 16
        if (ns > 0) goto 40
      elseif (first(1:2) == '-m') then ! Handle this switch in place
        j = 0
        call chrpos(first,'=',len(first)-1,j)
        if (j. lt. len(first)-1) first(j+1:j+1) = ' '
        call macset(first(3:),j)
        if (j. lt. 0) call rxs('mc: bad macro def: ',first)
      else if (first(1:5) == '-cmpf') then ! Handle this switch in place
        dc = first(6:6)
        if (dc == ' ') call rx('-cmpf missing arguments')
        facint = 0  ! Tolerance when comparing numbers
        i = parg2(first,dc,'tol=',',',dc,4,1,1,i1,facint,facint)
        ib(1) = 1 ; ib(2) = huge(i); ib(3:4) = [1,NULLI]
        i = parg2(first,dc,'ln=',',:',dc,2,1,2,i1,ib,ib)
        if (i>0 .and. i<2) then
          ib(2) = ib(1); ib(1) = 1
        endif
        idify = NULLI; i = parg2(first,dc,'diffy=',',:',dc,2,1,1,i1,idify,idify)
        i = parg2(first,dc,'col=',',:',dc,2,1,2,i1,ib(3),ib(3))
        if (i>0 .and. i<2) then
          ib(4) = ib(3); ib(3) = 1
        endif
        call rxx(any(ib(1:3)<=0),'row and col arguments to cmpf must be > 0')
        call rxx(ib(2)<ib(1),'illegal row specification')
        call rxx(ib(4)<ib(3).and.ib(4)/=NULLI,'illegal column specification')
        j = len(strn)
        if (ib(4) > j) call rx1('cmpf column width limited to %i characters',j)
!       call pshpr(iprint())
        if (wordsw(first,dc,'verb',dc//' ',i1) > 0) call setpr(40)
        if (wordsw(first,dc,'vverb',dc//' ',i1) > 0) call setpr(50)
        if (wordsw(first,dc,'quiet',dc//' ',i1) > 0) call setpr(1)
        ltmp = wordsw(first,dc,'char',dc//' ',i1) > 0 ! .true. => compare characters
        lmerge = wordsw(first,dc,'long',dc//' ',i1) > 0 ! .true. => Use longest line
        excl(1:mxexcl) = ' '
        k = 1; i2 = 1
        do while (wordsw(first(i2:),dc,'excl','',i1) > 0 .and. k < mxexcl)
          i1 = i1+i2
          dc2 = first(i1-1:i1-1)
          do  while (k <= mxexcl)
            call nwordg(first,0,dc//dc2//' ',1,i1,i2)
            excl(k) = first(i1:i2)
            k = k+1
            if (first(i2+1:i2+1) /= dc2) exit
            i1 = i2+2
          enddo
        enddo
        incl(1:mxexcl) = ' '
        if (wordsw(first,dc,'incl=','',i1) > 0) then
          i1 = i1-5-1
          do  k = 1, mxexcl
            i = index(first(i1:),dc//'incl=')
            if (i == 0) exit
            i1 = i1+i+5
            call nwordg(first,0,dc//' ',1,i1,i2)
            incl(k) = first(i1:i2)
            i1 = i2+1
          enddo
        endif
        ct = ' '; i = wordsw(first,dc,'sep=','',i1)
        if (i > 0) then
          call nwordg(first,0,dc,1,i1,i2)
          ct = first(i1:i2)
        endif
        trmstr = ' '; i = wordsw(first,dc,'term=','',i1)
        if (i>0) then
          call nwordg(first,0,dc,1,i1,i2)
          trmstr = first(i1:i2)
        endif

        i = wordsw(first,dc,'fn1=','',i1)
        call rxx(i==0,'cmpf requires fn1=')
        call nwordg(first,0,dc,1,i1,i2)
        ifi = fopng(first(i1:i2),-1,1)
        if (idify < 0) then
          i = wordsw(first,dc,'fn2=','',i1)
          call rxx(i==0,'cmpf requires fn2=')
          call nwordg(first,0,dc,1,i1,i2)
          jfi = fopng(first(i1:i2),-1,1)
        else
          idify = idify/2
        endif

        i = 0                   ! Line counter
        iwk(1) = 0              ! number of numerical or character mismatches
        iwk(2) = 0              ! number of words
        iwk(3) = 0              ! number of characters
        iwk(4) = 0              ! longest string encountered
        val(1) = 0              ! Maximum numerical difference
        n2a = 0; n2b = 0        ! < 0 => skip reading next line, to synchronize
        n3a = 0; n3b = 0        ! Current actual line number of 1st and 2nd file
C       call snot
        do  while (i<ib(2))     ! loop until ib(2) lines read or EOF
          if (n2a < 0 .or. n2b >=0) then  ! Read unless require file 2 read only
  911       continue
            if (.not. rdstrn(ifi,strn,len(strn),.false.)) exit
            do  k = 1, len(strn)
              if (strn(k:k) == TAB) strn(k:k) = ' '
            enddo
            n3a = n3a+1
            if (strn == ' ') goto 911
            if (trmstr /= ' ') then ! Check whether to terminate diff
              if (index(strn,trim(trmstr)) > 0) goto 35 ! We are done
            endif
            do  k = 1, mxexcl
              if (excl(k) == ' ') exit
              if (index(strn,trim(excl(k))) > 0) goto 911
            enddo
            if (incl(1) /= ' ') then
              do  k = 1, mxexcl
                if (incl(k) == ' ') goto 911 ! skip if incl missing
                if (index(strn,trim(incl(k))) > 0) exit
                if (k == mxexcl) goto 911 ! skip if incl missing
              enddo
            endif
            if (idify > 0) then
              strn2 = strn(idify+1:)
              strn  = strn(1:idify-1)
            endif
          endif
          if ((n2b < 0 .or. n2a >=0) .and. idify<0) then  ! Read unless require file 1 read only
  111       continue
            if (.not. rdstrn(jfi,strn2,len(strn2),.false.)) exit
            do  k = 1, len(strn2)
              if (strn2(k:k) == TAB) strn2(k:k) = ' '
            enddo
            n3b = n3b+1
            if (strn2 == ' ') goto 111
            if (trmstr /= ' ') then
              if (index(strn,trim(trmstr)) > 0) goto 35 ! No more checks
            endif
            do  k = 1, mxexcl
              if (excl(k) == ' ') exit
              if (index(strn2,trim(excl(k))) > 0) goto 111
            enddo
            if (incl(1) /= ' ') then
              do  k = 1, mxexcl
                if (incl(k) == ' ') goto 111 ! skip if incl missing
                if (index(strn2,trim(incl(k))) > 0) exit
                if (k == mxexcl) goto 111 ! skip if incl missing
              enddo
            endif
          endif
          n2a = 0; n2b = 0
          i = i+1
          if (i < ib(1)) cycle

C         Cull columns
          if (ib(4)/=NULLI) then
            slist = strn; strn = slist(ib(3):ib(4))
            slist = strn2; strn2 = slist(ib(3):ib(4))
          endif

C         Substitute ct for blank
          if (ct /= ' ') then
            do
              ms = scan(strn,trim(ct))
              if (ms == 0) exit
              strn(ms:ms) = ' '
            enddo
            do
              ms = scan(strn2,trim(ct))
              if (ms == 0) exit
              strn2(ms:ms) = ' '
            enddo
          endif

C     ... Count words in both strings; check whether new string should be read
          call words(strn,n1a); call words(strn2,n1b)
          if (n1a == 0) n2a = -1  ! Require new line file 1
          if (n1b == 0 .and. idify<0) n2b = -1  ! Require new line file 2
          if (n1a*n1b == 0 .and. idify<0) cycle
          if (lmerge) then
            iwk(2) = iwk(2) + max(n1a,n1b) ! increment total number of words
          else
            iwk(2) = iwk(2) + min(n1a,n1b) ! increment total number of words
            if (min(n1a,n1b) == 0) cycle
          endif

          i1 = len_trim(strn); i2 = len_trim(strn2)
          iwk(3) = iwk(3) + min(i1,i2) ! increment total number of characters
          iwk(4) = max(iwk(4),min(i1,i2))
          if (lmerge) iwk(4) = max(iwk(4),i1,i2)

C     ... Compare character by character
          if (ltmp) then
            do  j = 1, min(i1,i2)
              if (strn(j:j) /= strn2(j:j)) iwk(1) = iwk(1)+1
            enddo
            if (lmerge) iwk(1) = iwk(1) + abs(i1-i2)
            cycle
          endif

C     ... Convert each number-containing word to a numerical value; compare against tol
          call words(strn,it(1))
          i2 = 0
          do  ip = 1, it(1)
            i1 = i2+1
            call nword(strn,1,i1,i2)
            j = 0; wk(ip) = 0
            ms = scan(strn(i1:i2),'0123456789') ! must contain at least one digit
            if (ms == 0) cycle
            ms = verify(strn(i1:i2),'.+-0123456789')
            if (ms == 1) cycle
            ms = verify(strn(i1:i2),'.+-0123456789eEdD')
            if (ms /= 0) cycle
            if (.not. a2bin(strn(i1:),wk(ip),4,0,' ',j,-1)) wk(ip) = 0
          enddo
          call words(strn2,it(2))
          i2 = 0
          do  ip = 1, it(2)
            i1 = i2+1
            call nword(strn2,1,i1,i2)
            j = 0; wk(500+ip) = 0
            ms = scan(strn2(i1:i2),'0123456789') ! must contain at least one digit
            if (ms == 0) cycle
            ms = verify(strn2(i1:i2),'.+-0123456789')
            if (ms == 1) cycle
            ms = verify(strn2(i1:i2),'.+-0123456789eEdD')
            if (ms /= 0) cycle
            if (.not. a2bin(strn2(i1:),wk(500+ip),4,0,' ',j,-1)) wk(500+ip) = 0
          enddo
          i1 = it(1); i2 = it(2)
          ms = iwk(1)
          x0 = 0
          val(2) = val(1)
          do  ip = 1, min(i1,i2)
            val(1) = max(val(1),abs(wk(ip)-wk(500+ip)))
            x0 = max(x0,abs(wk(ip)-wk(500+ip)))
            if (abs(wk(ip)-wk(500+ip)) > facint) iwk(1) = iwk(1)+1
          enddo
C         call snot
          if (lmerge) iwk(1) = iwk(1) + abs(i1-i2)
          n = min(n1a,n1b); if (lmerge) n = max(n1a,n1b)

C         "very verbose" mode printout
C         Don't use info because strn might contain %
          if (iwk(1)>ms .and. iprint() >= 50) then
C            call info2(50,0,0,' %,5i : '//trim(strn),n3a,2)
C            call info2(50,0,0,' %,5i : '//trim(strn2),n3b,2)
            print 444, n3a, trim(strn)
            print 444, n3b, trim(strn2)
  444       format(1x,i5,a)
          endif

          if (iwk(1)>ms .or. iprint() >= 50) then
          call info8(40,0,0,
     .      ' line %,5i  words %,3i  ndiff %,3i  count diff %,3i  cumulative%,5i max %g'//
     .      '%?;n;  glob-max %g;%j;',
     .      min(n3a,n3b),n,iwk(1)-ms,abs(i1-i2),iwk(1),x0,isw(val(1)>val(2)),val(1))
          endif

        enddo ! do  while
   35   continue ! Exit do while

C   ... Printout
        call info0(50,1,1,' ----------------------------------')
        if (val(1) <= facint) val(1) = 0
        if (wordsw(first,dc,'nchar',dc//' ',i1) > 0) then
          call info2(30,0,0,'%i',iwk(3),2)
        elseif (wordsw(first,dc,'nword',dc//' ',i1) > 0) then
          call info2(30,0,0,'%i',iwk(2),2)
        elseif (wordsw(first,dc,'ndiff',dc//' ',i1) > 0) then
          call info2(30,0,0,'%i',iwk(1),2)
        elseif (wordsw(first,dc,'max',dc//' ',i1) > 0) then
          call info2(30,0,0,'%g',val(1),2)
        else
          call info8(30,0,0,
     .      ' lines %s:%2i  cols %s:%2i  words %i  chars %i  ndiff %i%?;n;  max %g;;',
     .    [ib(1),i],[ib(3),ib(3)-1+iwk(4)],iwk(2),iwk(3),iwk(1),isw(val(1)>0),val(1),8)
        endif

C   ... Exit value
!         call fexit(min(iwk(1),255),0,'',0)
        call fexit(0,0,'',0)
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15

C --- Load s(ns+1) from existing name, or file read ---
   30 continue
      do  31  i = 1, nnam
        if (first == symnam(i)) then
          call addnws(nrn(i),ncn(i),icn(i),ns,os,nr,nc,icast)
          j = nr(ns)*nc(ns)*sz(ns)
          call dcopy(j,w(symptr(i)),1,w(os(ns)),1)
          goto 32
        endif
   31 continue
C     Check read options from rdcsw
      if (rdcsw /= ' ') then
        ct = rdcsw(1:1)
        call partk0(0,len(rdcsw),1,-1,0,len(rdcsw),-1,31,.false.)
        j = partok(rdcsw,'qr',' '//ct,ltmp,' ',0,0,0,0)
        if (ltmp) rdops = 10*(rdops/10) + 1
        j = partok(rdcsw,'nr=','= '//ct,k,' ',-1,2,0,0)
        if (j == 1) nr(ns+1) = k
        j = partok(rdcsw,'nc=','= '//ct,k,' ',-1,2,0,0)
        if (j == 1) nc(ns+1) = k
        j = partok(rdcsw,'br',' '//ct,ltmp,' ',0,0,0,0)
        if (wordsw(rdcsw,ct,'fmt=','',j) /= 0) then
          rdfmt = rdcsw(j:)
        endif
        if (ltmp) then
          rdops = 10*(rdops/10) + 2
          j = partok(rdcsw,'br,',', ,'//ct,iwk,' ',-2,2,0,0)
          if (j == 2) then
            nr(ns+1) = iwk(1)
            nc(ns+1) = iwk(2)
            rdops = 10*(rdops/10) + 3
          endif
        endif
      endif
      if (first == '.' .or. (nargf() == iarg .and. ns == 0)) then
        ifi = i1mach(1)
      else
        i = 1
        if (mod(rdops,10) == 2 .or. mod(rdops,10) == 3) i = 5
c       call fshow
        ifi = fopnx(first,72,-1,-1)
        call word(first,1,i1,i2)
        if (ifi == 0) then
          print *, 'Exit -1 mcx: attempt to open nonexistent file '''//first(i1:i2) // ''''
          call cexit(-1,1)
        endif
        ifi = fopng(first,-1,i)
      endif
      call defask(i)
      i = i/i1mach(18) - 1
      ns = ns+1
      icast(ns) = 0
      call defdr(os(ns),i)
      filel = ' '
      if (rdcsw /= ' ') then
        ct = rdcsw(1:1)
        call partk0(0,len(rdcsw),1,-1,0,len(rdcsw),-1,31,.false.)
        j = partok(rdcsw,'s=','= '//ct,k,' ',-1,2,0,0)
        if (j == 1) then
          do  55  j = 1, k
            if (mod(rdops,10) >= 2) read(ifi)
            if (mod(rdops,10) < 2) read(ifi,*)
   55     continue
        endif
        j = partok(rdcsw,'open',' '//ct,noclos,' ',0,0,0,0)
        ltmp = .false.
        j = partok(rdcsw,'spc',' '//ct,ltmp,' ',0,0,0,0)
        if (ltmp) icast(ns) = 100
        rdcsw = ' '
      endif
      if (rdfmt /= ' ') then
C        nnz = 0
CC       return the number of nonzero elements
C        i = rdms(ifi,1000+rdops,nnz,filel,w,w(os(ns)),nr(ns),nc(ns))
C        if (i < 0) call rx('cannot read nr,nc from file')
        call redfrr(os(ns),nr(ns)*nc(ns)*2)
        i = 1  ! Real for now
        rewind ifi
        call readfmt('('//trim(rdfmt)//')',mod(i,10),w(os(ns)),nr(ns),nc(ns),ifi)
      elseif (icast(ns) >= 100) then ! sparse
        nnz = 0
C       return the number of nonzero elements
        i = rdms(ifi,1000+rdops,nnz,filel,w,w(os(ns)),nr(ns),nc(ns))
        osprs(1,ns) = nnz
        rewind ifi
        call rlse(os(ns))
        call defi(osprs(2,ns),2*(nnz+1))
        call defdr(os(ns),nnz*mod(i,10))
C       Read, with the appropriately allocated arrays
        i = rdms(ifi,1000+rdops,nnz,filel,w(osprs(2,ns)),
     .    w(os(ns)),nr(ns),nc(ns))
      else
        i = rdm(ifi,1000+rdops,i,filel,w(os(ns)),nr(ns),nc(ns))
        if (i == 0) then
          read(ifi) i,j
        endif
      endif
      if (mod(i,10) == 3) then
        call redfrr(os(ns),nr(ns)*nc(ns)*2)
        call defrr(owk,nr(ns)*nc(ns)*2)
        call dcopy(nr(ns)*nc(ns)*2,w(os(ns)),1,w(owk),1)
        call zcptoy(w(owk),nr(ns),nc(ns),w(os(ns)))
        call rlse(owk)
        i = i-1
      endif
      strn = ' '
      if (filel /= ' ') call awrit0('#mc read '''//first//'%a'': "'
     .  //filel//'%a"',strn,len(strn),-iout)
      if (i < 0)
     .  call fexit(-1,9,'mc failed to parse file '//first,0)
      if (ifi /= i1mach(1) .and. .not. noclos) call fclr(first,ifi)
      noclos = .false.
      icast(ns) = i
      call rlse(os(ns))
      call defrr(os(ns),nr(ns)*nc(ns)*sz(ns))
   32 continue
      if (ns >= nsmx) call rxi('mc: increase nsmx, ns=',ns)
      nr(ns+1) = 0
      nc(ns+1) = 0
      rdops = 10
C      outs = '#mc after loading '//first
C      call awrit4('%a: np: %n:1i sp: %n:1i',outs,80,-iout,
C     .  nnam,symptr,ns,os)

C      ibpo(1) = 3
C      ibpo(2) = 1
C      ibpo(3) = 2
C      ioffo(1) = 0+1
C      ioffo(2) = 1+1
C      ioffo(3) = 3+1
C      ioffo(4) = 8+1
C      ioffn(1) = 0+1
C      ioffn(2) = 2+1
C      ioffn(3) = 4+1
C      ioffn(4) = 7+1
C      call defrr(os(ns+1),nr(ns)*nc(ns)*sz(ns))
C      call dcopy(nr(ns)*nc(ns),w(os(ns)),1,w(os(ns+1)),1)
C      call dpzero(w(os(ns)),nr(ns)*nc(ns))
C      call pmblk(3,ibpo,ioffo,ioffn,nr(ns),w(os(ns+1)),8200,
C     .  one,nr(ns)-2,nr(ns),w(os(ns)),nc(ns))

   40 continue
      iarg = iarg+1

C --- Entry point for unary operations on existing matrix ---
   42 if (nuops < 1) goto 50
      if (ldebug) call awrit1(' mc: executing uop %i: '//
     .  uopstr(uopstk(1)),' ',-80,i1mach(2),uopstk)

      call defrr(os(ns+1),nr(ns)*nc(ns)*sz(ns))

C ... Add constant cshft to s(ns) ...
      if (uopstk(1) == 37) then
        if (cshft(2) /= 0d0) call nwcast(ns,os,nr,nc,icast,2)
        call dmadd(cshft,0,0,1d0,w(os(ns)),nr(ns),1,1d0,
     .    w(os(ns)),nr(ns),1,nr(ns),nc(ns))
        if (sz(ns) == 2) then
          call rx('not ready for complex shift')
        endif
        call clrop(uopstk,nuops)
      endif

C ... Scale s(ns) ...
      if (uopstk(1) == 1) then
        if (scale(2) == 0) then
          call dscal(nr(ns)*nc(ns)*sz(ns),scale,w(os(ns)),1)
        else
C        (cludge for now: use i1mach(18) for offset to imag)
C        (If real, convert to complex)
          if (scale(2) /= 0d0) call nwcast(ns,os,nr,nc,icast,2)
          i = nr(ns)*nc(ns)
          oiwk = os(ns) + i*i1mach(18)
          call yscal(i,scale(1),scale(2),w(os(ns)),w(oiwk),1)
          icast(ns) = 2
          scale(2) = 0
        endif
        call clrop(uopstk,nuops)
      endif

C ... Matrix subblock ...
C     lsubc
C       1   exclude subblock of s0
C       2   scales  subblock of s0
C       3   copies # to subblock of s0
C       4   extract subblock of s0
      if (uopstk(1) == 2) then
        call rlse(os(ns+1))
C   ... Possible convert real to complex matrix; eliminate symmetry
        if (lsubc == 2 .or. lsubc == 3) then
          icast(ns) = mod(icast(ns),10)
          if (scale(2) /= 0) call nwcast(ns,os,nr,nc,icast,2)
          ms = ns
        else
          call stkpsh(ns,ns,os,nr,nc,icast)
          ms = ns-1
        endif
        call mclip(lsubc,scale,w(os(ms)),sz(ns),nr(ns),nc(ns),w(os(ns)),nclip)
        call ivset(nclip,1,4,-1)
        if (lsubc == 1 .or. lsubc == 4) then
          call rlse(os(ns))
          call defrr(os(ns),nr(ns)*nc(ns)*sz(ns))
          icast(ns) = mod(icast(ns),10)
          call stktog(0,ns,os,nr,nc,icast)
          call clears(ns,nr,nc,icast)
          call defrr(os(ns+1),1)
        endif
        call clrop(uopstk,nuops)
      endif

C ... Replace s(ns) by expressions ...
      if (uopstk(1) == 4) then
        if (sz(ns) /= 1) call rx('expr subst not ready for cmplx')
        call redfrr(os(ns),nr(ns)*max(nc(ns),nexpr)*sz(ns))
        call defrr(os(ns+1),nr(ns)*nexpr*sz(ns))
        call defi(oiwk,nr(ns))
        call mapdat(nexpr,expr,' ',w(oiwk),nr(ns),nc(ns),
     .    w(os(ns)),w(os(ns+1)))
        call redfrr(os(ns),nr(ns)*nc(ns)*sz(ns))
        sincl = ' '
        nexpr = 0
        call clrop(uopstk,nuops)
        goto 42
      endif

C ... Symmetrize s(ns) ...
      if (uopstk(1) == 5) then
        call dosymm(sz(ns),w(os(ns)),nr(ns),nc(ns))
        icast(ns) = icast(ns) - 10*mod(icast(ns)/10,10) + 10
        call clrop(uopstk,nuops)
      endif

C ... Replace s(ns) by rows that satisfy sincl ...
      if (uopstk(1) == 6) then
        if (sz(ns) /= 1) call rx('incl not ready for cmplx')
        call defi(oiwk,nr(ns))
        call mapdat(0,expr,sincl,w(oiwk),nr(ns),nc(ns),
     .    w(os(ns)),w(os(ns+1)))
        call rlse(oiwk)
        sincl = ' '
        call clrop(uopstk,nuops)
      endif

C ... Replace s(ns) by its inverse ...
      if (uopstk(1) == 7) then
        if (nr(ns) /= nc(ns))
     .    call fexit2(-1,1,'Exit -1 mc: inv expected nc=%i '//
     .    'but found nc=%i',nr(ns),nc(ns))
        call defi(oiwk,nr(ns))
        call defdr(owk,2*nr(ns))
        call minv(lqi,sz(ns)/=2,nr(ns),w(os(ns)),w(oiwk),w(owk),j)
        if (j /= 0) call fexit(-1,1,'Exit -1 mc: matrix singular',0)
        call rlse(oiwk)
        call clrop(uopstk,nuops)
      endif

C ... Cholesky decomposition (or inverse Cholesky) of s(ns)  ...
      if (uopstk(1) == 47 .or. uopstk(1) == 48) then
        if (nr(ns) /= nc(ns))
     .    call fexit2(-1,1,'Exit -1 mc: Cholesky expects nc=%i '//
     .    'but found nc=%i',nr(ns),nc(ns))
        if (icast(ns) < 10) call rx('mc: Cholesky decomposition requires hermitian matrix')
        if (mod(icast(ns),10) /= 1) call rx('mc: Cholesky decomposition implemented only for real matrix')
        call defdr(owk,2*nr(ns))
        ltmp = uopstk(1) == 48
        call dschd(nr(ns),nr(ns),w(os(ns)),w(wk),ltmp,j)
        call clrtri(icast(ns),nr(ns),w(os(ns)))
        if (j /= 0) call fexit(-1,1,'Exit -1 mc: matrix singular',0)
        call rlse(owk)
        call clrop(uopstk,nuops)
      endif

C ... Replace s(ns) by its transpose ...
      if (uopstk(1) == 8) then
        call dcopy(nr(ns)*nc(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
        call mtrns(sz(ns),w(os(ns+1)),w(os(ns)),nr(ns),nc(ns))
        j = nr(ns)
        nr(ns) = nc(ns)
        nc(ns) = j
        call clrop(uopstk,nuops)
      endif

C ... Push onto stack, or push and name ...
      if (mod(iabs(uopstk(1)),100) == 9 .or. uopstk(1) == 16) then
        i = iabs(uopstk(1))/100
        if (uopstk(1) >= 0) i = ns-i
        if (i < 0 .or. i > ns) then
          call fexit(-1,1,'Exit -1 mc: invalid stack element, no. %i',i)
        else
          call rlse(os(ns+1))
          call stkpsh(i,ns,os,nr,nc,icast)
          nr(ns+1) = 0
          nc(ns+1) = 0
          call defrr(os(ns+1),1)
        endif
        if (uopstk(1) == 16) then
          uopstk(1) = 15
        else
          call clrop(uopstk,nuops)
        endif
      endif

C ... Pop s off stack ...
      if (uopstk(1) == 10) then
        call clears(ns,nr,nc,icast)
        call clrop(uopstk,nuops)
      endif

C ... Sort matrix ...
      if (uopstk(1) == 11) then
        if (sz(ns) /= 1) call rx('sort not ready for cmplx')
        call dcopy(nr(ns)*nc(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
        call defi(oiwk,nr(ns))
C       Put binary representation of expr 'sortex' into wk
        if (sz(ns) /= 1) call rx('expr subst not ready for cmplx')

C#ifdefC DVHTEST
C        print '(''# testing'')'
C        call dvheap(nr(ns),nc(ns),w(os(ns)),w(oiwk),0d0,0)
C#endif
        call dcopy(nr(ns)*nc(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
        nr(ns+1) = nr(ns)
        nc(ns+1) = nc(ns)
        call mapdat(1,sortex,' ',w(oiwk),nr(ns+1),nc(ns+1),
     .    w(os(ns)),w(1))
        call msort(nr(ns),nc(ns),icast(ns),w(oiwk),w(os(ns)),
     .    w(os(ns+1)),w(os(ns)))
        call clrop(uopstk,nuops)
      endif

C ... Replace s by its eigenvalues or eigenvectors ...
      if (uopstk(1) == 12 .or. uopstk(1) == 13) then
        if (nr(ns) /= nc(ns))
     .    call fexit2(-1,1,'Exit -1 mc: evl expected nr=nc '//
     .    'but found nr=%i and nc=%i',nr(ns),nc(ns))
C i=0 for evals, 1 for evecs
        i = 0
        if (uopstk(1) == 13) i = 1
C If nonsymmetric real, convert to complex
        if (icast(ns) == 1) call nwcast(ns,os,nr,nc,icast,2)
        call defcc(os(ns+1),nr(ns))
        call defcc(os(ns+2),nr(ns)*nc(ns))
        call defi(owk,nr(ns))
        icast(ns+1) = getev(nr(ns),w(os(ns)),w(os(ns)),icast(ns),
     .    .false.,i,w(os(ns+1)),w(os(ns+2)),w(owk))
        call rlse(owk)
        if (uopstk(1) == 12) then
          i = sz(ns+1)*nr(ns)
          nc(ns) = 1
        else
          os(ns+1) = os(ns+2)
          i = sz(ns+1)*nr(ns)**2
        endif
        call dcopy(i,w(os(ns+1)),1,w(os(ns)),1)
        call rlse(os(ns))
        call defrr(os(ns),i)
        call defrr(os(ns+1),1)
        icast(ns) = icast(ns+1)
        call clrop(uopstk,nuops)
      endif

C ... Expand vector into diagonal matrix or vice-versa ...
      if (uopstk(1) == 14) then
C       Expand diagonal part into matrix
        if (nc(ns) == 1) then
          call rlse(os(ns))
          i = nr(ns)**2*sz(ns)
          call redfrr(os(ns),i)
          call defrr(os(ns+1),nr(ns)*sz(ns))
          call dcopy(nr(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
          call dpzero(w(os(ns)),i)
          call v2dia(sz(ns),nr(ns),w(os(ns+1)),0,w(os(ns)))
          nc(ns) = nr(ns)
          icast(ns) = mod(icast(ns),10) + 10
C       Extract diagonal part from matrix
        elseif (nc(ns) == nr(ns)) then
          call v2dia(sz(ns),nr(ns),w(os(ns)),1,w(os(ns)))
          call redfrr(os(ns),nr(ns)*sz(ns))
          nc(ns) = 1
          icast(ns) = icast(ns) - mod(icast(ns),100) + mod(icast(ns),10)
          call defrr(os(ns+1),1)
        else
          call fexit(-1,1,'Exit -1 mc: v2dia expected nc=1 '//
     .      'but found nc=%i',nc(ns))
        endif
        call clrop(uopstk,nuops)
      endif

C ... Assign name to top stack element ...
      if (uopstk(1) == 15) then
        call rlse(os(ns+1))
C   ... See if synmam already exists; if so use inam = nnam
        inam = nnam+1
        do  i = 1, nnam
          if (symnam(nnam+1) == symnam(i)) inam = i
        enddo
C   ... Name did not exist already; create a new nam
        if (inam == nnam+1) nnam = nnam+1
        call nams(ns,os,nr,nc,icast,symptr(inam),nrn(inam),ncn(inam),
     .    icn(inam))
        call clrop(uopstk,nuops)
        call defrr(os(ns+1),1)
      endif

C ... Sum columns or rows ...
      if (uopstk(1) == 17 .or. uopstk(1) == 18) then
        if (uopstk(1) == 17) then
          call mtrns(sz(ns),w(os(ns)),w(os(ns+1)),nr(ns),nc(ns))
          i = nc(ns)
          j = nr(ns)
          nc(ns) = 1
        else
          call dcopy(nr(ns)*nc(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
          i = nr(ns)
          j = nc(ns)
          nr(ns) = 1
        endif

        if (sums(1:1) /= ' ') then
          nr(ns) = i; nc(ns) = j ! undo the change just made
          do  ipass = 1, 2 ! First pass => just count number of rows (colums)

          md = 0 ! panel index
          mxnum = 0 ! maximum size of list
          npseq = 0 ! number of panels obtained through seq=, needed for wordg below
          iseq = 0  ! Sequence index, 0 if not in sequence
          do ! until no more lists
            call wordg(sums,1,'; ',md-npseq+1,is1,is2)
            doslst = ' '
            if (is2 < is1) exit
            doslst = sums(is1:is2)

            if (iseq == 0 .and. doslst(1:4) == 'seq=') then
              if (nnum == 0) call rx('mcx: sequence declared before list')
              nnseq = mkilsd(doslst(5:),-1,iwk)
              if (nnseq <= 0) call rx('mcx: bad sequence')
              allocate(nseq(nnseq))
              nnseq = mkilsd(doslst(5:),nnseq,nseq)
              iseq = 1
              do  while (iseq > 0)
                md = md+1
                if (ipass == 2) then
                  num(1:nnum,md) = num(1:nnum,md-1) + nseq(iseq) - num(1,md-1)
                endif
                iseq = mod(iseq+1,nnseq+1)
                if (iseq > 0) then
                  npseq = npseq + 1 ! Used in finding the next string, see wordg above
                else
                  deallocate(nseq)
                endif
              enddo
            else
              nnum = mkilsd(doslst,-1,iwk)
              if (nnum <= 0) call rx('mcx: bad list')
              md = md+1
              mxnum = max(mxnum,nnum)
              if (ipass == 1) cycle
              nnum = mkilsd(doslst,nnum,num(1,md))
            endif
          enddo
          if (ipass == 1) then
            if (uopstk(1) == 17) call rx('not ready for summing columns')
            allocate(arr(md,nc(ns),sz(ns)),num(mxnum,md))
            call iinit(num,size(num))
          endif
          enddo
          call sumrx(sz(ns),w(os(ns)),nr(ns),nc(ns),arr,md,nc(ns),num,mxnum)
          nr(ns) = md
          icast(ns) = mod(icast(ns),10)
          call rlse(os(ns))
          call defrr(os(ns),nr(ns)*nc(ns)*sz(ns))
          call dcopy(size(arr),arr,1,w(os(ns)),1)
          deallocate(arr,num)
        else
          iwk(1) = 1
          iwk(2) = i
          call sumr(sz(ns),w(os(ns+1)),iwk(1),iwk(2),i,j,w(os(ns)))
        endif
        call defrr(os(ns+1),1)
        call clrop(uopstk,nuops)
      endif

C ... Replace matrix by its complex congugate ...
      if (uopstk(1) == 19) then
        if (sz(ns) == 2) then
          call dscal(nr(ns)*nc(ns)*2,-one,w(os(ns)),1)
          call dscal(nr(ns)*nc(ns),-one,w(os(ns)),1)
        else
          call awrit0('#cc (warning): no cc of real matrix',' ',80,iout)
        endif
        call clrop(uopstk,nuops)
      endif

C ... Customize ...
      if (uopstk(1) == 20) then
C       call cust(ns,osprs,os,nr,nc,icast)
        call clrop(uopstk,nuops)
        call defrr(os(ns+1),1)
      endif

C ... Interpolate to a spec'd list of points
      if (uopstk(1) == 21) then
        call itrpsw(intrps,it)
        i = mkdlst(slist,-1024*d1mach(3),0,w(os(ns+1)))
        call defrr(os(ns+1),i*nc(ns)*icast(ns))
        nr(ns+1) = mkdlst(slist,-1024*d1mach(3),i,w(os(ns+1)))
        nc(ns+1) = nc(ns)
        icast(ns+1) = icast(ns)
        call intrp(w(os(ns)),nr(ns),nc(ns),icast(ns),w(os(ns+1)),
     .    nr(ns+1),it)
        call stktog(0,ns+1,os,nr,nc,icast)
C       call intrp(slist,it,ns,os,nr,nc,icast)
        call clrop(uopstk,nuops)
      endif

C ... Differentiate wrt col 1 ...
      if (uopstk(1) == 22) then
        call itrpsw(intrps,it)
        if (.not. diffx(it,w(os(ns)),nr(ns),nc(ns),icast(ns)))
     .    call rx('mc failed differentiate file')
        call clrop(uopstk,nuops)
      endif

C ... Integrate col2 wrt col1 from start to a list of upper limits
      if (uopstk(1) == 23) then
        call itrpsw(intrps,it)
        i = mkdlst(slist,-1024*d1mach(3),0,w(os(ns+1)))
        nc(ns+1) = nc(ns)
        icast(ns+1) = 1
        if (it(1) == 2) i = nr(ns)
        call defrr(os(ns+1),-i*nc(ns+1))
        nr(ns+1) = mkdlst(slist,-1024*d1mach(3),i,w(os(ns+1)))
        call intgx(intrps,x0,it,i,w(os(ns+1)),
     .    w(os(ns)),nr(ns),nc(ns),icast(ns))
        nr(ns+1) = i
        call stktog(0,ns+1,os,nr,nc,icast)
        call clrop(uopstk,nuops)
      endif

C ... Real part of a complex matrix
      if (uopstk(1) == 24) then
        if (mod(icast(ns),10) == 2) icast(ns) = icast(ns)-1
        i = nr(ns)*nc(ns)*sz(ns)
        call redfrr(os(ns),i)
        call defrr(os(ns+1),1)
        call clrop(uopstk,nuops)
      endif

C ... Force a matrix hermitian or symmetric
      if (uopstk(1) == 28) then
        call rxx(nr(ns)/=nc(ns),
     .    'nonsquare matrix cannot be made hermitian')
C       copy s(ns) to s(ns+1), replace s by its transpose
        call dcopy(nr(ns)*nc(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
        call mtrns(sz(ns),w(os(ns+1)),w(os(ns)),nr(ns),nc(ns))
C       complex congugate of transpose matrix
        if (sz(ns) == 2) then
          call dscal(nr(ns)*nc(ns)*2,-one,w(os(ns)),1)
          call dscal(nr(ns)*nc(ns),-one,w(os(ns)),1)
        endif
C       add original to cc-tranpose
        call daxpy(nr(ns)*nc(ns)*sz(ns),one,w(os(ns+1)),1,w(os(ns)),1)
C       and scale for average
        call dscal(nr(ns)*nc(ns)*sz(ns),one/2,w(os(ns)),1)
        i = nr(ns)*nc(ns)*sz(ns)
        call redfrr(os(ns),i)
        call defrr(os(ns+1),1)
        icast(ns) = icast(ns) + 10*(1-getdig(icast(ns),1,10))
        call clrop(uopstk,nuops)
      endif

C ... Create a new array out of rows and columns of top stack
C     wk(1) -> flag for coll:e
      if (uopstk(1) == 29 .or. uopstk(1) == 30) then
        nr(ns+1) = nr(ns)
        nc(ns+1) = w(owk)
        if (uopstk(1) == 29) then
          nr(ns+1) = w(owk)
          nc(ns+1) = nc(ns)
        endif
        icast(ns+1) = sz(ns)
        call redfrr(os(ns+1),nr(ns+1)*nc(ns+1)*sz(ns+1))
        call crowl(uopstk(1)==29,wk(1)/=0,w(owk),w(owk+1),sz(ns),
     .    w(os(ns)),nr(ns),nc(ns),w(os(ns+1)),nr(ns+1),nc(ns+1))
        call defps2(owk,os(ns+1))
        call rlse(owk)
        call stktog(0,ns+1,os,nr,nc,icast)
        call clrop(uopstk,nuops)
      endif

C ... Replace list of points with a smoothed Gaussian
      if (uopstk(1) ==  31) then
        if (wk(4) == 0) call rx('smooth: dx is zero')
        nr(ns+1) = nint((wk(3)-wk(2))/wk(4)+1)
        nc(ns+1) = 2
        icast(ns+1) = icast(ns)
        call redfrr(os(ns+1),nr(ns+1)*nc(ns+1)*sz(ns+1))
        call smooth(wk,w(os(ns)),nr(ns),nc(ns),w(os(ns+1)),nr(ns+1),
     .    nc(ns+1))
        call stktog(0,ns+1,os,nr,nc,icast)
        call clears(ns,nr,nc,icast)
        call clrop(uopstk,nuops)
      endif

      if (uopstk(1) ==  32) then
        call unx(wk,w(os(ns)),nr(ns),nc(ns),sz(ns))
        call clrop(uopstk,nuops)
      endif

      if (uopstk(1) ==  39) then
        call unx2(wk,w(os(ns)),nr(ns),nc(ns),sz(ns))
        call clrop(uopstk,nuops)
      endif

      if (uopstk(1) ==  40) then
        call rxx(sz(ns)/=1,'-at not allowed w/ cmplex array')
        call findx(0,wk,sincl,w(os(ns)),nr(ns),nc(ns),w,nr(ns+1))
        nc(ns+1) = 2
        icast(ns+1) = icast(ns)
        call redfrr(os(ns+1),nr(ns+1)*nc(ns+1)*sz(ns+1))
        call findx(1,wk,sincl,w(os(ns)),nr(ns),nc(ns),
     .    w(os(ns+1)),nr(ns+1))
        ns = ns+1
        call defrr(os(ns+1),1)
        call clrop(uopstk,nuops)
      endif

      if (uopstk(1) ==  33) then
        call abss(w(os(ns)),nr(ns),nc(ns),icast(ns))
        call clrop(uopstk,nuops)
      endif

C ... max value in an array
      if (uopstk(1) ==  34) then
        call maxs(isw(wk(1)==1)+10*isw(wk(1)==2),
     .    w(os(ns)),nr(ns),nc(ns),sz(ns))
        if (wk(1) == 0 .or. wk(1) == 2) then
          nc(ns) = 1
          call redfrr(os(ns),nr(ns)*nc(ns)*sz(ns))
        elseif (wk(1) == 1) then
          nr(ns) = 1
          nc(ns) = 3
          icast(ns) = 1
          call redfrr(os(ns),nr(ns)*nc(ns)*sz(ns))
        endif
        call defrr(os(ns+1),1)
        call clrop(uopstk,nuops)
      endif

C ... replicate array
      if (uopstk(1) ==  35) then
        if (nrep(1)*nrep(2) > 0) then
          nr(ns+1) = nr(ns)
          nc(ns+1) = nc(ns)
          icast(ns+1) = icast(ns)
          nr(ns) = nrep(1)*nr(ns)
          nc(ns) = nrep(2)*nc(ns)
          icast(ns) = sz(ns)
          call redfrr(os(ns),nr(ns)*nc(ns)*sz(ns))
          call defrr(os(ns+1),nr(ns+1)*nc(ns+1)*sz(ns+1))
          call dcopy(nr(ns)*nc(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
          call mcrep(icast,nr(ns+1),nrep(1),nc(ns+1),nrep(2),
     .      w(os(ns+1)),w(os(ns)))
        endif
        call clrop(uopstk,nuops)
      endif

C ... roll array
      if (uopstk(1) ==  38) then
        if (nroll(1)/=0 .or. nroll(2) /= 0) then
          nr(ns+1) = nr(ns)
          nc(ns+1) = nc(ns)
          icast(ns+1) = icast(ns)
          icast(ns) = sz(ns)
          call defrr(os(ns+1),nr(ns+1)*nc(ns+1)*sz(ns+1))
          call dcopy(nr(ns)*nc(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
          call mcroll(icast,nr(ns+1),nroll(1),nc(ns+1),nroll(2),
     .      w(os(ns+1)),w(os(ns)))
          nroll(1) = 0
          nroll(2) = 0
        endif
        call clrop(uopstk,nuops)
      endif

C ... assign a variable to top array on stack
      if (uopstk(1) ==  36) then
        if (sz(ns) /= 1) call rx('var. assign not ready for complex')
        xx = 0
        if (vassgn(1) <= nr(ns) .and. vassgn(2) <= nc(ns))
     .  xx = dval(w(os(ns)),vassgn(1)+nr(ns)*(vassgn(2)-1))
        call lodsyv(scalnm,1,xx,k)
        call clrop(uopstk,nuops)
      endif

C ... Split matrices at top of stack into subblocks
      if (uopstk(1) == 25) then
        nlistx = mkilsd(splitx,-1,ilsx)
        if (nlistx<0 .or. nlistx>size(ilsx)) goto 19
        nlisty = mkilsd(splity,-1,ilsy)
        if (nlisty<0 .or. nlisty>size(ilsy)) goto 19
        nlistx = mkilsd(splitx,nlistx,ilsx)
        nlisty = mkilsd(splity,nlisty,ilsy)

        call rlse(os(ns+1))
        do  251  i = 2, nlistx
        do  251  j = 2, nlisty
C     ... Set the clip limits
          nclip(1) = ilsx(i-1)
          nclip(2) = ilsx(i)-1
          nclip(3) = ilsy(j-1)
          nclip(4) = ilsy(j)-1
C     ... Push stack, replace top by subblock
          call stkpsh(ns,ns,os,nr,nc,icast)
          call mclip(4,0d0,w(os(ns-1)),sz(ns),nr(ns),nc(ns),w(os(ns)),
     .      nclip)
          call redfrr(os(ns),nr(ns)*nc(ns)*sz(ns))
          icast(ns) = mod(icast(ns),10)
          call ivset(nclip,1,4,-1)
C     ... Name it, push off top
          nnam = nnam+1
          call awrit2('%x'//splitn//'%a%i%i',symnam(nnam),
     .      namlen,0,i-1,j-1)
          call nams(ns,os,nr,nc,icast,symptr(nnam),nrn(nnam),ncn(nnam),
     .      icn(nnam))
  251   continue
        call clrop(uopstk,nuops)
        call defrr(os(ns+1),1)
      endif

C ... Raise matrix to a power
      if (uopstk(1) == 27) then
        call dcopy(nr(ns)*nc(ns)*sz(ns),w(os(ns)),1,w(os(ns+1)),1)
        call pwrs(nr(ns),nc(ns),icast(ns),pwr,w(os(ns+1)),w(os(ns)))
        call clrop(uopstk,nuops)
      endif

C ... Nearest integer
      if (uopstk(1) == 41) then
        call tonint(w(os(ns)),nr(ns),nc(ns),icast(ns),facint)
        call clrop(uopstk,nuops)
      endif

C ... Rotation matrix of spherical harmonics
      if (uopstk(1) == 46) then
        call rlse(os(ns+1))
C       call wkprnt(1); call wkinfo

        if (ns<1 .or. nr(ns) /= 3 .or. nc(ns) /= 3)
     .    call rx('-ylm requires that 3x3 rotation matrix resides on stack')
        dc = ylms(1:1)
        k = wordsw(ylms,dc,'l=','',j) + 2
        if (dc == ' ' .or. k < 3) call rx('missing or improper switches following -ylm')
        if (a2vec(ylms,len_trim(ylms),k,2,', '//dc,3,2,1,it,j) < 1)
     .      call rx('mc failed to parse '//trim(ylms))
        ip = (j+1)**2

C       call wkprnt(1); call wkinfo

        if (wordsw(ylms,dc,'sh','',k) > 0) then
          i = 1; if (wordsw(ylms,dc,'sh2','',k) > 0) i = 2
        endif

C       Rotation matrix for spin 1/2
        if (wordsw(ylms,dc,'spin','',k) > 0) then
          call rm2eua(w(os(ns)),wk(1),wk(2),wk(3))
          call addnws(ip*2,ip*2,2,ns,os,nr,nc,icast)
          ms = 1; if (wordsw(ylms,dc,'spin+o','',k)>0) ms = ms+2
          if (wordsw(ylms,dc,'spin+oonly','',k)>0) ms = 2
          call makerots(ms+10*i,j+1,wk,w(os(ns)))
          call ztoyy(w(os(ns)),nr(ns),nc(ns),nr(ns),nc(ns),1,0)
          call stktog(0,ns,os,nr,nc,icast)
          call clears(ns,nr,nc,icast); call rlse(os(ns+1))
        else
          call addnws(ip,ip,1,ns,os,nr,nc,icast)
          call ylmrtg(ip,w(os(ns-1)),w(os(ns)))
          call stktog(0,ns,os,nr,nc,icast)
          call clears(ns,nr,nc,icast); call rlse(os(ns+1))
C         call wkprnt(1); call wkinfo

          if (wordsw(ylms,dc,'sh','',k) > 0) then
            i = 2-i
            call addnws(ip,ip,2,ns,os,nr,nc,icast)
            call s2sph(0+100*i,j+1,j+1,w(os(ns-1)),ip,ip,ip,ip,w(os(ns)))
            call stktog(0,ns,os,nr,nc,icast)
            call clears(ns,nr,nc,icast); call rlse(os(ns+1))
C           call wkprnt(1); call wkinfo
          endif

        endif

        call clrop(uopstk,nuops)
        goto 42
      endif

C ... Lagrange interpolating polynomial
      if (uopstk(1) == 49) then
        if (nc(ns) /= 1) call rx('mcx requires a single column for Lagrange interpolation')
        call ipoly(nr(ns),xitrp,w(os(ns)),w(os(ns+1)))
        call clrop(uopstk,nuops)
      endif

C ... Replace s0 with LU decomposition of it
      if (uopstk(1) == 50) then
        if (nr(ns) /= nc(ns)) call rx('mc: sought LU decomposition of nonsquare matrix')
        i = icast(ns)
        call nwcast(ns,os,nr,nc,icast,2)
        call ztoyy(w(os(ns)),nr(ns),nc(ns),nr(ns),nc(ns),0,1)
        call ludcmp(w(os(ns)),nr(ns),icast(ns))
        call ztoyy(w(os(ns)),nr(ns),nc(ns),nr(ns),nc(ns),1,0)
        if (icast(ns)-i == 1) call nwcast(ns,os,nr,nc,icast,1)
        icast(ns) = icast(ns) - 10*getdig(icast(ns),1,10)
        call clrop(uopstk,nuops)
        goto 42
      endif

C ... Clean up for unary operators ...
      nr(ns+1) = 0
      nc(ns+1) = 0
      call rlse(os(ns+1))
      goto 42

C --- Entry point for binary operations on top 2 matrices ---
   50 continue
      lmerge = .false.
      if (ns < 2 .or. nbops < 1) goto 60
      if (bopstk(nbops) == 15 .and. ns < 4) goto 15

      if (ldebug) call awrit1(' mc: executing binop %i: '//
     .  bopstr(bopstk(1)),' ',-80,i1mach(2),bopstk)

C ... Row or column concatenation ...
      if (bopstk(nbops) == 1 .or. bopstk(nbops) == 2) then
C       call wkdbg
        call rccatm(bopstk(nbops),ns,icast,nr,nc,os)
        nbops = nbops-1
        lmerge = .true.

C ... Addition and subtraction ...
      elseif (bopstk(nbops) == 3 .or. (bopstk(nbops) == 4)) then
        if (nr(ns-1) /= nr(ns))
     .    call fexit2(-1,1,'Exit -1 mc: add expected nr=%i '//
     .    'but found nr=%i',nr(ns-1),nr(ns))
        if (nc(ns-1) /= nc(ns))
     .    call fexit2(-1,1,'Exit -1 mc: add expected nc=%i '//
     .    'but found nc=%i',nc(ns-1),nc(ns))
        if (sz(ns-1) /= sz(ns)) then
          call awrit0('#mc (warning): cast mismatch ...',' ',80,iout)
          call sscast(ns,os,nr,nc,icast)
        endif
        i = min(sz(ns-1),sz(ns))
        call defrr(os(ns+1),nr(ns)*nc(ns)*i)
        val(1) = one
        if (bopstk(nbops) == 4) val(1) = -one
        call dmadd(w(os(ns-1)),nr(ns-1),1,one,w(os(ns)),nr(ns),1,val,
     .    w(os(ns+1)),nr(ns-1),1,nr(ns-1),nc(ns-1)*i)
        nr(ns+1) = nr(ns-1)
        nc(ns+1) = nc(ns-1)
        icast(ns+1) = i
        if (icast(ns-1) == icast(ns)) icast(ns+1) = icast(ns)
        nbops = nbops-1
        lmerge = .true.

C ... Multiplication ...
      elseif (bopstk(nbops) == 5) then
        if (nc(ns-1) /= nr(ns))
     .    call fexit2(-1,1,'Exit -1 mc: mpy expected nr=%i '//
     .    'but found nr=%i',nc(ns-1),nr(ns))
        if (sz(ns-1) /= sz(ns)) then
          call awrit0('#mc (warning): cast mismatch ...',' ',80,iout)
          call sscast(ns,os,nr,nc,icast)
        endif
        call defrr(os(ns+1),nr(ns-1)*nc(ns)*sz(ns))
        if (sz(ns) == 1) then
          call dmpy(w(os(ns-1)),nr(ns-1),1,w(os(ns)),nr(ns),1,
     .      w(os(ns+1)),nr(ns-1),1,nr(ns-1),nc(ns),nc(ns-1))
          icast(ns+1) = 1
        else
          nr(ns+1) = nr(ns-1)
          nc(ns+1) = nc(ns)
          call zmpy(w(os(ns-1)),nr(ns-1),1,nr(ns-1)*nc(ns-1),
     .      w(os(ns)),nr(ns),1,nr(ns)*nc(ns),w(os(ns+1)),nr(ns+1),1,
     .      nr(ns+1)*nc(ns+1),nr(ns-1),nc(ns),nc(ns-1))
          icast(ns+1) = 2
        endif
        nr(ns+1) = nr(ns-1)
        nc(ns+1) = nc(ns)
        nbops = nbops-1
        lmerge = .true.

C ... Outer product
      elseif (bopstk(nbops) == 16) then
        stop 'not ready for xo'

C ... 3D Multiplication ...
C     Treats arrays (ns-1) ns as 3D arrays A and B, dimensioned
C       B(n1b=nr(ns)/n2b,n2b=nc(ns-1),n3b=nc(ns))
C       A(n1a=n1b,n2a=nr(ns-1)/nr(ns),n3a=n2b)
C     and creates C(nr=n1a*n2a,nc=n3b):
C       and creates C(n1a,n2a,nc) = C(nr=n1a*n2a,nc=n3b) (2D form):
C       C(i,j,k) = sum_m A(i,j,m) B(i,m,k)
      elseif (bopstk(nbops) == 12) then

        n1b = nr(ns)/nc(ns-1)
        n2b = nc(ns-1)
        n3b = nc(ns)
        n1a = n1b
        n2a = nr(ns-1)/n1a
        n3a = nc(ns-1)
        if (n2a*n1a /= nr(ns-1))
     .    call fexit2(-1,1,'Exit -1 mc: mpy3 nr(ns-1)=%i '//
     .    'not an integer multiple of nr(ns)=%i',nr(ns-1),nr(ns))
        if (sz(ns-1) /= sz(ns)) then
          call awrit0('#mc (warning): cast mismatch ...',' ',80,iout)
          call sscast(ns,os,nr,nc,icast)
          call rx('mpy3 not ready for size mismatch')
        endif
        call defrr(os(ns+1),n1a*n2a*n3b*sz(ns))
        if (sz(ns) == 1) then
          call mpy3(n1a,n2a,n3a,n1b,n2b,n3b,
     .      w(os(ns-1)),w(os(ns)),w(os(ns+1)))
          icast(ns+1) = 1
        else
          call ztoyy(w(os(ns-1)),nr(ns-1),nc(ns-1),nr(ns-1),nc(ns-1),
     .      0,1)
          call ztoyy(w(os(ns)),nr(ns),nc(ns),nr(ns),nc(ns),0,1)
          call zmpy3(n1a,n2a,n3a,n1b,n2b,n3b,
     .      w(os(ns-1)),w(os(ns)),w(os(ns+1)))
          icast(ns+1) = 2
          call ztoyy(w(os(ns-1)),nr(ns-1),nc(ns-1),nr(ns-1),nc(ns-1),
     .      1,0)
          call ztoyy(w(os(ns)),nr(ns),nc(ns),nr(ns),nc(ns),1,0)
          nr(ns+1) = n1a*n2a
          nc(ns+1) = n3b
          call ztoyy(w(os(ns+1)),nr(ns+1),nc(ns+1),nr(ns+1),nc(ns+1),
     .      1,0)
        endif
        nr(ns+1) = n1a*n2a
        nc(ns+1) = n3b
        nbops = nbops-1
        lmerge = .true.

C ... Toggle top elements on stack ...
      elseif (mod(bopstk(nbops),100) == 6) then
        i = iabs(bopstk(nbops))/100
        call stktog(i,ns,os,nr,nc,icast)
        nbops = nbops-1

C ... Gen. eval problem: top two matrices replaced by its evals or evecs
      elseif (bopstk(nbops) == 7 .or. bopstk(nbops) == 8) then
        if (nr(ns) /= nc(ns))
     .    call fexit2(-1,1,'Exit -1 mc: evl expected nr=nc '//
     .    'but found nr=%i and nc=%i',nr(ns),nc(ns))
        if (nr(ns-1) /= nc(ns-1))
     .    call fexit2(-1,1,'Exit -1 mc: evl expected nr=nc '//
     .    'but found nr=%i and nc=%i',nr(ns-1),nc(ns-1))
        call defrr(os(ns+1),nr(ns)*sz(ns))
        call defrr(os(ns+2),nr(ns)*nc(ns)*sz(ns))
        i = 0
        if (bopstk(nbops) == 8) i = 1
        call defi(owk,nr(ns))
        icast(ns+1) = getev(nr(ns),w(os(ns)),w(os(ns-1)),icast(ns),
     .    .true.,i,w(os(ns+1)),w(os(ns+2)),w(owk))
        call rlse(owk)
        if (bopstk(nbops) == 7) then
          i = sz(ns+1)*nr(ns)
          nc(ns) = 1
        else
          os(ns+1) = os(ns+2)
          i = sz(ns+1)*nr(ns)**2
        endif
        nr(ns+1) = nr(ns)
        nc(ns+1) = nc(ns)
        nbops = nbops-1
        lmerge = .true.

C ... Overwrite s0 with S1^-1/2  S0  S1^-1/2
C     bug: should check thst S1 is hermitian is
      elseif (bopstk(nbops) == 13) then
        if (icast(ns) /= 12 .or. icast(ns-1) /= 12)
     .    call fexit(-1,1,'Exit -1 mc: orthos implemented only '//
     .    'for hermitian s0 and s1',0)
        if (nr(ns) /= nc(ns))
     .    call fexit2(-1,1,'Exit -1 mc: orthos expected nr=nc '//
     .    'but found nr=%i and nc=%i',nr(ns),nc(ns))
        if (nr(ns-1) /= nc(ns-1))
     .    call fexit2(-1,1,'Exit -1 mc: evl expected nr=nc '//
     .    'but found nr=%i and nc=%i',nr(ns-1),nc(ns-1))
        call defrr(owk,2*nr(ns))
        call diagno(nr(ns),w(os(ns)),w(os(ns-1)),w(owk),.false.,3,0,0,
     .    0d0,i,0d0,0d0)
        call rlse(owk)
        call stktog(i,ns,os,nr,nc,icast)
        call clears(ns,nr,nc,icast)
        nbops = nbops-1

C ... Overwrite s0(i) with s1(s0(i))
      elseif (bopstk(nbops) == 14) then
        if (sz(ns) /= 1) call rx('index matrix cannot be complex')
        nr(ns+1) = nr(ns)
        nc(ns+1) = nc(ns-1)
        icast(ns+1) = mod(icast(ns-1),10)
        call defrr(os(ns+1),nr(ns+1)*nc(ns+1)*sz(ns+1))
        call s0mps1(w(os(ns-1)),nr(ns-1),nc(ns-1),icast(ns+1),w(os(ns)),
     .    nr(ns),.false.,w(os(ns+1)))
        nbops = nbops-1
        lmerge = .true.

C ... Multiply or divide arrays element by element, or cross product
      elseif (bopstk(nbops) >= 9 .and. bopstk(nbops) <= 11) then
        if (nr(ns-1) /= nr(ns) .and. nr(ns-1) /= 1)
     .    call fexit2(-1,1,'Exit -1 mc: op expected nr=%i '//
     .    'but found nr=%i',nr(ns-1),nr(ns))
        if (nc(ns-1) /= nc(ns))
     .    call fexit2(-1,1,'Exit -1 mc: op expected nc=%i '//
     .    'but found nc=%i',nc(ns-1),nc(ns))
        if (sz(ns-1) /= sz(ns)) then
          call awrit0('#mc (warning): cast mismatch ...',' ',80,iout)
          call sscast(ns,os,nr,nc,icast)
        endif
        i = sz(ns)
        call defrr(os(ns+1),nr(ns)*nc(ns)*sz(ns))
        nr(ns+1) = nr(ns)
        nc(ns+1) = nc(ns)
        icast(ns+1) = i
        if (bopstk(nbops) == 11) then
          call rxx(sz(ns)/=1,'no cross product on cmplex vectors')
          call crss(w(os(ns-1)),w(os(ns)),nr(ns-1),nr(ns+1),nc(ns+1),
     .      w(os(ns+1)))
        else
          call cmprd(w(os(ns-1)),w(os(ns)),nr(ns+1),nc(ns+1),
     .      one,i==1,.false.,.false.,bopstk(nbops)==10,w(os(ns+1)))
        endif
        nbops = nbops-1
        lmerge = .true.

C ... Copy subblock of one matrix into subblock of another
      elseif (bopstk(nbops) == 17) then
        call mscop(lsubc,ns,os,nr,nc,icast,nclip,scale)
        call ivset(nclip,1,4,-1)
        nbops = nbops-1

C ... Row + column concatenation ...
      elseif (bopstk(nbops) == 15) then
C       Keep bookeeping data for ns-1 since it will be overwritten
        i1 = icast(ns-1)
        i2 = os(ns-1)
        i = nr(ns-1)
        j = nc(ns-1)
C       call wkdbg
C       Concatenate columns of s3, s2
        call rccatm(1,ns-2,icast,nr,nc,os)
C       copy new bookeeping data into ns+2; restore ns-1 data
C       Array pointer in ns-2 no longer needed :
C       use for temp storage and restore ns-1 pointers ...
        icast(ns-2) = icast(ns-1)
        os(ns-2) = os(ns-1)
        nr(ns-2) = nr(ns-1)
        nc(ns-2) = nc(ns-1)
        icast(ns-1) = i1
        os(ns-1) = i2
        nr(ns-1) = i
        nc(ns-1) = j
C       Concatenate columns of s1, s0
        call rccatm(1,ns,icast,nr,nc,os)
C       Copy ns+1 info to ns-1; shrink stack by 2 elements
        icast(ns-1) = icast(ns+1)
        os(ns-1)    = os(ns+1)
        nr(ns-1)    = nr(ns+1)
        nc(ns-1)    = nc(ns+1)
        ns = ns-2
C       Concatenate rows of newly created matrices
        call rccatm(2,ns+1,icast,nr,nc,os)
C       Copy ns+2 info to ns+1
        icast(ns+1) = icast(ns+2)
        os(ns+1)    = os(ns+2)
        nr(ns+1)    = nr(ns+2)
        nc(ns+1)    = nc(ns+2)

        nbops = nbops-1
        lmerge = .true.
      endif

C ... Copies s(ns+1) to s(ns-1) and reduces stack size by 1 ...
      if (lmerge) then
        call mergesn(1,ns,icast,nr,nc,os)
C        call rlse(os(ns-1))
C        call defrr(os(ns-1),nr(ns+1)*nc(ns+1)*sz(ns+1))
C        call dcopy(nr(ns+1)*nc(ns+1)*sz(ns+1),
C     .    w(os(ns+1)),1,w(os(ns-1)),1)
C        nr(ns-1) = nr(ns+1)
C        nc(ns-1) = nc(ns+1)
C        icast(ns-1) = icast(ns+1)
C        call clears(ns,nr,nc,icast)
      endif

C --- Cleanup operators ---
   60 continue
      if (ns > 1 .and. nbops > 0) goto 50
C      if (ns > 1 .and. nbops > 0) then
CC       Special case : requires 4 matrices on stack
C        if (bopstk(nbops) == 15 .and. ns < 4) goto 15
CC       Else evaluate next binop
C        goto 50
C      endif
      goto 15

C --- Write top of stack on exit ---
  999 continue
      if (ns == 0) goto 30
      if (prmx /= 0) then
        call xprm(nprec,icast(ns),iout*(3-2*prmx),osprs(1,ns),w(os(ns)),nr(ns),nc(ns))
        prmx = 0
      else
        call prm(0,' ',icast(ns),iout,fmt,osprs(1,ns),w(os(ns)),nr(ns),nc(ns))
      endif

      end

      subroutine readfmt(rdfmt,i,s,nr,nc,ifi)
      character*(*) rdfmt
      integer ifi
      double precision s(nr,nc)

      read(ifi,trim(rdfmt)) s
  333 format(1p,4e20.13)

      end
C      subroutine snot
C      end
