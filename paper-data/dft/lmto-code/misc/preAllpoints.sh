#! /bin/sh
# This script prepare several files for the irrk calculation
# allpoints points : required by lmgf
# allpoints_m : required by  matlab code
#
# Usage:
#       lmf fe --rs=1,0 --pr80 -vnit=0 --rpos| tee out.tmp    # prepare the output file
#       preAllpoints.sh  # text filter
#

iallqStartM2=` grep -n "BZMESH: qp mapping" out.tmp | cut -f1 -d:`
iirrqStartM2=`grep -n "irreducible QP from" out.tmp | cut -f1 -d:`
iirrqFinishP1=`grep -n "TETIRR: sorting"    out.tmp | cut -f1 -d:`

iallqStart=$(( $iallqStartM2+2))
iallqFinish=$(( $iirrqStartM2-1))
iirrqStart=$(( $iirrqStartM2+2))
iirrqFinish=$(( $iirrqFinishP1-1))

sed -n "$iallqStart, $iallqFinish p" out.tmp  > allpoints
sed -n "$iirrqStart, $iirrqFinish p" out.tmp  > points
awk '{$7=""}1' allpoints  | tr -s '\( ,\), \,' ' '  > allpoints_m

# iallq
# echo $iy
# echo 'ix=' $ix  'iy=' $iy 'iz=' $iz
# awk '{print $7}' allpoints
# cat allpoints |tr -s '\( ,\), \,' ' '
# cut -f1-6  allpoints
# awk '{for(i=1;i<=NF;i++)if(i!=7)printf$i OFS;print""}' allpoints  | tr -s '\( ,\), \,' ' ' | tr -s '\)' '    '
