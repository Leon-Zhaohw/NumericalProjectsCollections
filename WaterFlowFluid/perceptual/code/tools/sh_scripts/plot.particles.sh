#!/bin/bash -e
# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2015 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL)
# http://www.gnu.org/licenses
#
# Gnuplot-based rendering script for 2D particles
#
# ----------------------------------------------------------------------------

W=600				# width
H=600				# height
TRANS=0				# translation
SCALE=1.0			# scaling
RX=64				# x domain [0,RX]
RY=64				# y domain [0,RY]
C_DIM=1				# num of elements for color data
VEC=0				# draw vectors
TCUT=0				# remove type?
RT=999				# then, remove this type
O_NAME="particles.png"		# output
C_NAME="particles.pT.sep"	# this data will decide color
T_NAME="particles.pT.sep"	# type data
V_NAME="particles.pV.sep"	# vector data
RANGE_A=0			# minimum value in color range
RANGE_B=1			# maximum value in color range
VSCALE=1			# scaling factor for vector data
C_SET="MATLAB"
BG="white"
FG="black"

function echooption() {
    echo -e -n "\t\t[\033[1m$1\033[0m ($2)]\n"
}
function usage() {
    echo "Usage: $(basename $0) [options] space-delimiter-position-data-file"
    echo -e "\toptions:"
    echooption "-a value" "min value in color range"
    echooption "-b value" "max value in color range"
    echooption "-c filepath" "filepath for color data"
    echooption "-t filepath" "filepath for type data; use with -r option"
    echooption "-v filepath" "filepath for vector (arrow) data"
    echooption "-o filepath" "output PNG filepath"
    echooption "-r value" "remove this type; valid only when -t option is given (this will be used with sed command, so it can be a regexp string)"
    echooption "-d 2/3" "color data type; 2D/3D vector, length is used"
    echooption "-e value" "scaling factor for vector data"
    echooption "-f white/grey/black/#1234AB/..." "foreground color; default is black"
    echooption "-g white/grey/black/#1234AB/..." "background color; default is white"
    echooption "-s MATLAB/BLUES/YGB" "color preset; default is MATLAB"
    echooption "-h int" "image height"
    echooption "-w int" "image width"
    echooption "-x value" "domain x"
    echooption "-y value" "domain y"
    echooption "-S value" "scale data"
    echooption "-T value" "translate data"
}

if [ $# -eq 0 ]; then usage; exit $?; fi

while getopts ":a:b:c:d:e:f:g:h:o:r:s:t:v:w:x:y:S:T:" optname; do
    case "$optname" in
	"a") RANGE_A=$OPTARG ;;
	"b") RANGE_B=$OPTARG ;;
	"c") C_NAME=$OPTARG ;;
	"d") C_DIM=$OPTARG ;;
	"e") VSCALE=$OPTARG ;;
	"f") FG=$OPTARG ;;
	"g") BG=$OPTARG ;;
	"h") H=$OPTARG ;;
	"o") O_NAME=$OPTARG ;;
	"r") TCUT=1; RT=$OPTARG ;;
	"s") C_SET=$OPTARG ;;
	"t") T_NAME=$OPTARG ;;
	"v") VEC=1; V_NAME=$OPTARG ;;
	"w") W=$OPTARG ;;
	"x") RX=$OPTARG ;;
	"y") RY=$OPTARG ;;
	"S") SCALE=$OPTARG ;;
	"T") TRANS=$OPTARG ;;
	*)
	    echo "Unknown option: -$OPTARG"
	    exit $?
	    ;;
    esac
done
shift $(($OPTIND - 1))

if [[ ( $RANGE_B < $RANGE_A ) ]]; then
    C_SET=$C_SET"_REV"
fi

I_NAME=$@

[[ -e $I_NAME ]] || (echo "Can't find file: $I_NAME" 1>&2 && exit 1)
[[ -e $C_NAME ]] || (echo "Can't find file: $C_NAME" 1>&2 && exit 1)
[[ $TCUT -eq 1 ]] && ([[ -e $T_NAME ]] || (echo "Can't find file: $T_NAME" 1>&2 && exit 1))
[[ $VEC -eq 1 ]]  && ([[ -e $V_NAME ]] || (echo "Can't find file: $V_NAME" 1>&2 && exit 1))

[[ ${I_NAME##*.} == "gz" ]] && I_NAME="<(zcat $I_NAME)"
[[ ${C_NAME##*.} == "gz" ]] && C_NAME="<(zcat $C_NAME)"
[[ ${T_NAME##*.} == "gz" ]] && T_NAME="<(zcat $T_NAME)"
[[ ${V_NAME##*.} == "gz" ]] && V_NAME="<(zcat $V_NAME)"

gnuplot <<PLOT
unset key

set output "$O_NAME"
set terminal png size $W, $H background rgb '$BG'

# change a color of border.
set border lc rgb '$FG'

# change text colors of	 tics
set xtics textcolor rgb '$FG'
set ytics textcolor rgb '$FG'

set macros
MATLAB=	     "defined (0 '#000090', 1 '#000FFF', 2 '#0090FF', 3 '#0FFFEE', 4 '#90FF70', 5 '#FFEE00', 6 '#FF7000', 7 '#EE0000', 8 '#7F0000')"
MATLAB_REV = "defined (0 '#7F0000', 1 '#EE0000', 2 '#FF7000', 3 '#FFEE00', 4 '#90FF70', 5 '#0FFFEE', 6 '#0090FF', 7 '#000FFF', 8 '#000090')"
BLUES=	     "defined (0 '#F7FBFF', 1 '#DEEBF7', 2 '#C6DBEF', 3 '#9ECAE1', 4 '#6BAED6', 5 '#4292C6', 6 '#2171B5', 7 '#084594')"
BLUES_REV=   "defined (0 '#084594', 1 '#2171B5', 2 '#4292C6', 3 '#6BAED6', 4 '#9ECAE1', 5 '#C6DBEF', 6 '#DEEBF7', 7 '#F7FBFF')"
YGB=	     "defined (0 '#FFFFD9', 1 '#EDF8B1', 2 '#C7E9B4', 3 '#7FCDBB', 4 '#41B6C4', 5 '#1D91C0', 6 '#225EA8', 7 '#0C2C84')"
YGB_REV=     "defined (0 '#0C2C84', 1 '#225EA8', 2 '#1D91C0', 3 '#41B6C4', 4 '#7FCDBB', 5 '#C7E9B4', 6 '#EDF8B1', 7 '#FFFFD9')"

datacut =  ($TCUT==0) ? "'< exec bash -c \"cut -d \\\  -f 1-2 $I_NAME | paste - $C_NAME\"'" : "'< exec bash -c \"cut -d \\\  -f 1-2 $I_NAME | paste $T_NAME - $C_NAME | sed \"/^$RT/d\" - | cut -f 2- -\"'"
vdatacut = ($TCUT==0) ? "'< exec bash -c \"cut -d \\\  -f 1-2 $I_NAME | paste - $V_NAME\"'" : "'< exec bash -c \"cut -d \\\  -f 1-2 $I_NAME | paste $T_NAME - $V_NAME | sed \"/^$RT/d\" - | cut -f 2- -\"'"
pointplot = "@datacut using (\$1*$SCALE+$TRANS):(\$2*$SCALE+$TRANS):($C_DIM==1 ? \$3 : ($C_DIM==2 ? sqrt(\$3*\$3+\$4*\$4) : sqrt(\$3*\$3+\$4*\$4+\$5*\$5))) with points pt 5 ps 0.2 palette"
plotstr = ($VEC==0) ? pointplot : "@pointplot, @vdatacut using (\$1*$SCALE+$TRANS):(\$2*$SCALE+$TRANS):(\$3*$VSCALE):(\$4*$VSCALE) with vectors head size 0.08,20,60 filled lw 1"

unset title
unset xlabel
unset ylabel
unset colorbox
unset key

# define grid
set style line 12 lc rgb '$FG' lt 0 lw 1
set grid back ls 12

set xrange [0:$RX]
set yrange [0:$RY]
set cbtics scale 0
set cbrange [$RANGE_A:$RANGE_B]
set palette @$C_SET
plot @plotstr

quit
PLOT
