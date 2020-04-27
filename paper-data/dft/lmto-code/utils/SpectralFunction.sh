#!/bin/bash

#------------------------------------------------------
# A script to generate the spectral function plots.
#------------------------------------------------------

version="1.1 beta"

export GDFONTPATH=/usr/share/fonts/bitstream-vera
# User may have to edit this path

#----------------------------------------------------
function show_help() #Help and usage
{

detailed=false
if [ ! -z $1 ] ; then
  if [ $1 == "H" ] ; then detailed=true;fi
fi

printf "
A script to visualize and save the spectral function.\n
USAGE : `basename $0` -<OPTIONS>

OPTIONS:

    -v | --viz | --visualize : To visualize the spectral function.
Note that visualization is NOT enabled by default, but saving is."
if $detailed ; then printf "
User will be given option to save at the end."
fi
printf "

Presently the script can save output in three formats: PNG, SVG and EPS.\n\
Use the following options to use them:
   -p : PNG format (General purpose. This is the DEFAULT option)
   -e : EPS format (Useful for LaTeX)
   -s : SVG format (In case you need better control and quality)
In order to save the output in more than one format use, for example :
   `basename $0` -se
"
if $detailed ; then  printf "
If the script is run without any options following is assumed:
   `basename $0` -p
producing only PNG output.
"
fi

printf "
There are additional options:
   -n : To view and save the negative part of the spectral function."
if $detailed ; then  printf "
Ideally Spectral Function should not have any negative parts, but sometimes
numerical procedure may introduce it." ; fi

printf "

   -f : --factor : Shrink the colorbar range by a factor N.
                   Number N must be follow the -f flag."
if ! $detailed ; then printf " (See --info for details.)
      "
else
printf "

   This option is useful when the details of the spectral function are
   washed out in overall range. By reducing the range the details at the
   bottom of the colorbar are reveled. As an example if, the colorbar
   range is at 0:1000 then '-f 4' option will shrink it to 0:250. NOTE
   that only colorbar range is changed; meaning values above 250 are not
   discarded but are represented by the same color.

";fi

printf "
   -w : Draw the spectral function over an energy range e1:e2.
        argument e1:e2 must follow the -w flag."
if ! $detailed ; then printf " (See --info for details.)
      "
else
printf "
   e1 and e2 are the lower and upper energy bounds for the figure;
   each must be a real number.
";fi

printf "
   -h| --help : This help message.
   -i| --info : Extended help and details of procedure.
   --faq      : Troubleshoot. Possible problems and their solutions.
   --version  : Version information
"
if $detailed ; then
printf "
Spectral function is read from the file 'spf.<system>'. The name of the
system is extracted from 'ctrl.<system>' file. In case the script fails to
access ctrl file, it will ask user for the correct system name.

Note that the script will process specfun.<system> file and produce
'specfun.pm3d' file, which will be later used by gnuplot. The script will also
generate 'gnu.plt' file, which is executed by gnuplot. If desired, user may
tweak the file manually to change the appearance of a particular output.

Note on the output: The script will produce a few output images depending upon
the flags.  It will use a 'base name' for the output files. Typically the base
name is same as your system name. The script will extract the system name from
the extension of 'ctrl' file. Optionally user can provide a different base name
(without whitespaces). At the end the output files will look like
'<basename>_UP.png' and '<basename>_DN.png', etc.

In case the negative part is asked for it will be stored in the separate files:
'<basename>_UP_negative.png ' and '<basename>_DN_negative.png', etc.
"
fi

printf "
Use --info for extended help and --faq for troubleshoot.\n
Bugs/Comments?: Bhalchandra Pujari, pujari@unl.edu\n
"
      exit 0
}


#--------------------------------------------------------------------

function faq()  # Troubleshoot
{
      printf "
FAQ:

* What are the prerequisite for this script?
  - On a standard Linux platform you are likely to have everything installed.
  Along with a standard Bourne Shell, it requires Gnuplot with support to PNG,
  SVG and EPS terminals. It also uses programs like AWK, SED.

* What are the units?
  - Spectral function is in the arbitrary units. That is why color bar has no
  units. The energy (vertical axis) is in electron volts (eV). As is the case
  in standard band-structures, horizontal axis has no units but is scaled as
  per the distances among the symmetry points in brillouin zone.

* I get a message 'Could not find/open font when opening font Vera, Using
 default'.What does it mean?
  - This is a non-critical issue. In the PNG mode Gnuplot searches for 'Vera'
  font. If not found it switches to default font. Vera is nicer than the
  default.  You can install Vera font to get rid of this message. In case you
  have it installed and continue getting the error, change 'export
  GDFONTPATH=<path>' statement at the beginning of the script to point to the
  directory where the font is located. Note that your EPS and SVG outputs are
  not affected by this error.

* What is the 'specfun.pm3d' file? Is it OK to delete it?
  - The raw data from specfun.<system> is processed in order to visualize it in
  gnuplot. The data in pm3d file should be same except for some additional
  blank lines. YES, it is safe to remove .pm3d file. More information on pm3d
  format can be found in Gnuplot manual.

* Most of the plot appears in the same color. What do I do?
  - Sometime a few maxima of the spectral function are so high in value that
  in comparison, other values do not stand out in the plot. In that case it is
  useful to shrink the colorbar range to reveal the details. '-f' flag would do
  just that. Note that -f must be followed by a number, which is used as a
  factor to shrink the colorbar. Also note that although the values on the
  colorbar change, it does not mean the higher values are discarded. It simply
  means they are shown in the same color as that of maxima of colorbar.

* What does -f flag do?
  - Please check the previous question.

* I don't like the colors of the plot. How can I change them?
  - There are two ways. First: After you run the script, a 'gnu.plt' file would
  have been generated. You can change the definition of 'palette' in the file
  and run 'gnuplot gnu.plt', which will change the appearance.  Second: You can
  change the script itself to change the definition of the palette. Note
  though, that the second method is permanent and will change the appearance of
  all future plots. Gnuplot manual will give you more info on how to define
  the palettes.

* What is SVG? What is so special about it?
  - SVG is Scalable Vector Graphics format. It allows more flexibility on the
  post-processing of the output image. For example you can change the axis
  properties or resize the plots without loosing the image quality. SVG files
  can be handled using programs like 'inkscape'. It can be saved in any other
  graphics format.


Use `basename $0` --info for detailed help.

"


exit 0

}

#--------------------------------------------------------------------

function show_error()
{

      detailed=false
            if [ ! -z $2 ] ; then
                  if [ $2 == "U" ] ; then detailed=true;fi
            fi
      printf "\nError! $1 \n"

      if $detailed ; then
      printf "Please use -h for usage and --faq for troubleshoot.\n"
      fi


      exit 1


}

#--------------------------------------------------------------------

function test_command()
{
      for a in $* ; do
            which $a >/dev/null 2>/dev/null
            if [ $? != 0 ] ; then
                  show_error "Command not found: $a \nIs $a in your PATH?\n"
            fi
      done

}

#--------------------------------------------------------------------
function show_version()
{
      printf "`basename $0`: Version - $version\nBug reports to: pujari@unl.edu\n\n"; exit 0
}

isnumber() #test if number
{
      test "$2" && printf '%f' "$2" > /dev/null 2>&1;

      if [ $? -ne 0 ] ; then
            show_error "'$2' is not a number. $1 flag must be associated with a number.\n" "U"
      fi
}

isnumberf() #test if number
{
echo $1

      test "$1" && printf '%f' "$1" > /dev/null 2>&1;

      if [ $? -ne 0 ] ; then
            show_error "'$cbfactor' is not a number. -f flag must be associated with a number.\n" "U"
      fi
}

#--------------------------------------------------------------------
#-----MAIN PROGRAM----------------------

#-Check for installed packages.
test_command gnuplot awk cat sed tr


printf "\n"
pngout=false
epsout=false
svgout=false
skipsvg=true
skipeps=true
skippng=false
negative=false
showhelp=false
showinfo=false
showversion=false
showfaq=false
defaulton=false
visual=false
cbfactor="1.0"


while getopts "pPsSeEnNhHiIvV-:f:w:" flag ; do    # Until you run out of parameters . . .
  case "$flag" in
   p| P)
      pngout=true
      skippng=false
      ;;
   s| S)
      svgout=true
      skipsvg=false
      ;;
   e| E)
      epsout=true
      skipeps=false
      ;;
   n| N)
      negative=true
      ;;
   h| H)
      showhelp=true
      ;;
   i| I)
      showinfo=true
      ;;
   v| V)
     visual=true
      ;;
   f| F)
      cbfactor=$OPTARG
      ;;
   w| W)
      window=$OPTARG
      ;;
   -)
       case "${OPTARG}" in
            factor)
               cbfactor=$OPTARG
               ;;
            [vV]is|[vV]isualise|[vV]iz)
               visual=true
               ;;
            [fF][aA][qQ])
               showfaq=true
               ;;
            [vV]ersion)
               showversion=true
               ;;
            [hH]elp)
               showhelp=true
               ;;
            [iI]nfo)
               showinfo=true
               ;;
            *)
              show_error "Invalid argument --$OPTARG" "U"
       esac
       ;;
   *)
      show_error "Invalid argument -$OPTARG" "U"
  esac
done
  shift $((OPTIND-1))

      [ "$flag" = "--" ] && shift


isnumber -f $cbfactor
if [ `echo "$cbfactor >= 1 "| bc -lq` -ne 1 ] ; then
      show_error "Factor must be greater than 1\n"
fi

if ! [ -z ${window+x} ]; then
isnumber -w `echo $window | awk -F : '{print $1}'`
isnumber -w `echo $window | awk -F : '{print $2}'`
window='['$window'][][]'
fi

if $showinfo ; then show_help "H" ; fi
if $showversion ; then show_version ; fi
if $showhelp ; then show_help ; fi
if $showfaq ; then faq ; fi

if [ `ls ctrl.* | wc -l` -ne 1 ] ; then
    printf "Oops! No unique ctrl file found. Please enter the system name: \n"
    read base
 else
    ctrlname=`ls ctrl*`
    base=`echo ${ctrlname#*.}`
fi

output="$base"
specfile="spf.${base}"
if [ ! -s $specfile ] ; then show_error "File $specfile does not exist" "U" ;fi
#printf "Converting spectral function to gnuplot's format\n\n"
awk '/^[[:blank:]]*#/ {next}; NF < 3 {next}; $1 != prev {printf "\n"; prev=$1}; {print}' $specfile > specfun.pm3d
if [ $? != 0 ] ; then
   show_error "Oops! Something went wrong! Have you finished calculating spectral function?\n\
If you think the $specfile is fine, consider reporting this issue to developers."
fi

if ! $pngout && ! $svgout && ! $epsout ; then
      defaulton=true
fi

while : ; do  #infinite while-do loop

if  $visual   ; then  #only visualization
      printf "Visualizing the spectral function. \n"
      keepread="n"
else
      keepread="y"
      switch=true

      if $defaulton ; then
      pngout=true
      printf "No output format specified, using default.\n"
      fi

      printf "Output format(s):  "
      if $pngout ; then printf "PNG " ; fi
      if $epsout ; then printf "EPS " ; fi
      if $svgout ; then printf "SVG " ; fi
      printf "\n"
      if $negative; then
        printf "Both positive and negative parts of spectral function will be saved.\n"
      fi
      printf "\n"
fi


#When $visual switch is on producing output is skipped with keepread=n
while [ ${keepread} == "y" ] ; do #keep reading user input till the flag changes
  if $switch ; then
  printf "The basename for the output files will be:\n\n   $output \n\nSounds good? (y/n) "
  read ans
  else
  ans="n"
#  keepread="n"
  fi


  if [ -z  ${ans} ] ; then
    keepread="n"
  elif [ $ans == "n" ] ; then
    printf 'Please enter the basename for output files: (Do not include the extension like '.eps')\n'
    read output
    output=`echo $output | tr -d ' '`
    if [ -z  ${output} ] ; then
     show_error "No input provided. Exiting."
    fi
   else         #anything other than ans='n' is considered 'y'
    keepread="n"
  fi


  upfile="${output}_UP"
  dnfile="${output}_DN"
  up_negative_file="${output}_UP_negative"
  dn_negative_file="${output}_DN_negative"


  if $svgout ; then
    if [[ -f "$upfile.svg" ]]; then
       printf "\n NOTE:SVG output already exists. Skipping SVG.\n"
       skipsvg=true
#       keepread="y"
     else
        keepread="n"
     fi
  fi
  if $pngout ; then
    if [[ -f "$upfile.png" ]]; then
       printf "\n NOTE: PNG output already exists. Skipping PNG.\n"
       skippng=true
#       keepread="y"
     else
       keepread="n"
     fi
  fi

  if $epsout ; then
    if [[ -f "$upfile.eps" ]]; then
       printf "\n NOTE: EPS output already exists. Skipping EPS.\n"
       skipeps=true
#       keepread="y"
     else
       keepread="n"
     fi
  fi

  if  $skipsvg && $skippng && $skipeps  ; then #if all the formats already present
    printf "\nOops! Looks like all the output files are already present.\n\
You can either suggest a new name  or I can overwrite them.\n\
Do you want me to overwrite?(y/[n]/q): "
    read overwrt
    if [ -z  ${overwrt} ] ; then
      keepread="y"
    elif [ $overwrt == "y" ] ; then
      keepread="n"
      printf "WARNING! Output files will be overwritten\n"
    elif [ $overwrt == "q" ] ; then
      printf "Quitting\n"
      exit 0
    else
      keepread="y"
    fi
  fi
      switch=false
      skipeps=false
      skippng=false
      skipsvg=false

done



gnufile="gnu.plt"
rm -f $gnufile



#set symmetry lines
for pts in `head -n1 $specfile | sed -e 's/#//'` ; do
  echo "set arrow from  graph 0,first ${pts}, graph 0 to graph 1,first ${pts},0 nohead front lt -1" >> $gnufile
done

#find max-min

upmax=`sort -n -k3 $specfile | tail -n1 | awk '{print $3}'`
dnmax=`sort -n -k4 $specfile | tail -n1 | awk '{print $4}'`


upcbrange=`echo "(${upmax}/3.0)/${cbfactor}" | bc -lq`
dncbrange=`echo "(${dnmax}/3.0)/${cbfactor}" | bc -lq`

if  $visual  ; then
  hashcharcter="#"
  hashcharcter2=" "
  pauseline="pause -1"
else
  hashcharcter=" "
  hashcharcter2="#"
  pauseline=" "
  printf "Please wait while I generate the output files..\n\n"
fi


cat >> $gnufile << EOF
set pm3d at b
unset surf
set view 180,90
unset ytics
#set palette defined (0 "white", .1 "blue", 1 "red");
set palette rgbformulae -10,-13,-26
set yzeroax lt 2
#set zr [0:]

EOF

if $pngout || $defaulton || $visual ; then
cat >> $gnufile << EOF
set cbrange resto
set cbrange[0:$upcbrange]
$hashcharcter set terminal png enhanced  nocrop  font "Vera" size 800,600
$hashcharcter set out "${upfile}.png"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:3
$hashcharcter unset out
$hashcharcter2 print "Now showing the UP channel"
$hashcharcter2 print "Close the window and hit enter to view Down channel"
$pauseline
set cbrange resto
set cbrange[0:$dncbrange]
$hashcharcter set terminal png enhanced   nocrop font "Vera" size 800,600
$hashcharcter set out "${dnfile}.png"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:4
$hashcharcter unset out
$hashcharcter2 print "Close the window and hit enter to proceed."
$pauseline
EOF
fi

if $svgout && ! $visual ; then
cat >> $gnufile << EOF
set cbrange resto
set cbrange[0:$upcbrange]
set term svg fsize 22
set out "${upfile}.svg"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:3
unset out
set cbrange resto
set cbrange[0:$dncbrange]
set term svg fsize 22
set out "${dnfile}.svg"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:4
unset out
EOF
fi

if $epsout && ! $visual ; then
cat >> $gnufile << EOF
set cbrange resto
set cbrange[0:$upcbrange]
set term pos enh col eps 24
set out "${upfile}.eps"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:3
unset out
set cbrange resto
set cbrange[0:$dncbrange]
set term pos enh col eps 24
set out "${dnfile}.eps"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:4
unset out
EOF
fi


if $negative ; then
#find max-min on negative side

upmax=`sort -nr -k3 $specfile | tail -n1 | awk '{print $3}'`
dnmax=`sort -nr -k4 $specfile | tail -n1 | awk '{print $4}'`


upcbrange=`echo "(${upmax})/${cbfactor}" | bc -lq`
dncbrange=`echo "(${dnmax})/${cbfactor}" | bc -lq`
cat >> $gnufile << EOF
set auto
set zr [:0]
EOF

if $pngout || $defaulton || $visual  ; then
cat >> $gnufile << EOF
set cbrange resto
set cbrange[$upcbrange:0]
$hashcharcter set terminal png enhanced  nocrop  font "Vera" size 800,600
$hashcharcter set out "${up_negative_file}.png"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:3
$hashcharcter unset out
$hashcharcter2 print "Now showing the negative part of spectracl funnction in UP channel"
$hashcharcter2 print "Close the window and hit enter to view Down channel"
$pauseline
set cbrange resto
set cbrange[$dncbrange:0]
$hashcharcter set terminal png enhanced  nocrop  font "Vera" size 800,600
$hashcharcter set out "${dn_negative_file}.png"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:4
$hashcharcter unset out
$hashcharcter2 print "Close the window and hit enter to proceed."
$pauseline
EOF
fi

if $svgout && ! $visual ; then
cat >> $gnufile << EOF
set cbrange resto
set cbrange[$upcbrange:0]
set term svg fsize 22
set out "${up_negative_file}.svg"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:3
unset out
set cbrange resto
set cbrange[$dncbrange:0]
set term svg fsize 22
set out "${dn_negative_file}.svg"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:4
unset out
EOF
fi

if $epsout && ! $visual ; then
cat >> $gnufile << EOF
set cbrange resto
set cbrange[$upcbrange:0]
set term pos enh col eps 24
set out "${up_negative_file}.eps"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:3
unset out
set cbrange resto
set cbrange[$dncbrange:0]
set term pos enh col eps 24
set out "${dn_negative_file}.eps"
sp $window "specfun.pm3d" u (\$1*13.60569806):2:4
unset out
EOF
fi
fi # if-loop for $negative

gnuplot $gnufile

if [ $? -ne 0 ]; then
  show_error "Looks like something went wrong while plotting."
fi

if  $visual  ; then
   printf "Visualizing done. Do you want to save the output? (y,[n]) \n"
   read vis_ans
   if [ -z  ${vis_ans} ] ; then
     printf "No input.\n"
     break
   elif [ $vis_ans == "y" ] ; then
     printf "Proceeding to save the output\n"
     visual=false
   else
     printf "Plots are not saved.\n"
     break
   fi
else
   break
fi


done

if [ $? -ne 0 ]; then
  show_error "Looks like something went wrong. This shouldn't have happened.\n\
Please report the issue with maximum possible information. Exit status: $?"
fi

printf "\nEverything done! Now exiting.\n"

exit 0

#---------------------------------------------------------
# A script by Bhalchandra Pujari. Feb 2013.
# Major modification on March 1st 2013
# Bugs? Report it to : pujari@unl.edu
#---------------------------------------------------------
