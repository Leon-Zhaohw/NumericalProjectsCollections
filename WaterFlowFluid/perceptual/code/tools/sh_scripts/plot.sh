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
# Frame data plotting script
#
# ----------------------------------------------------------------------------

function echooption() {
    echo -e -n "\t\t[\033[1m$1\033[0m ($2)]\n"
}
function usage() {
    echo "Usage: $(basename $0) [options] datafiles ..."
    echo -e "\toptions:"
    echooption "-t" "total; default is average"
    echooption "-a" "magnitude for scalar or length square for vector"
    echooption "-d filename" "filename for particle type data; it should be in the same directory with datafiles"
    echooption "-p type#" "used with -d option only; select data type using this number, other types are excluded"
    echooption "-g graph-title" "used only for pdf file format"
    echooption "-l data-label" "used for matlab and pdf file formats"
    echooption "-c 1/2/3" "pick a component"
    echooption "-v 1/2/3" "data dimension"
    echooption "-w weight" "scaling value; multiplied to each data element"
    echooption "-o filename.{pdf,txt,m}" "output filename; supported types are pdf (gnuplot), text, and matlab"
    echooption "-s" "open the output file after plotting"
}
function clip() {
    if [ ! $PTYPE_FILE == "" ]; then
	typefile="$(dirname $1)/$PTYPE_FILE"
	datafile="$1"
	[[ -e $typefile ]] || (echo "Can't find file: $typefile" 1>&2 && exit $?)
	[[ ${typefile##*.} == "gz" ]] && typefile="<(zcat $typefile)"
	[[ ${datafile##*.} == "gz" ]] && datafile="<(zcat $datafile)"
	eval paste $typefile $datafile | sed '/^[^'${SELECT_PTYPE}']/d' | cut -f 2-
    else
	([[ ${1##*.} == "gz" ]] && zcat $1) || cat $1 || exit $?
    fi
}
function select_scalar() {
    ls -1 "$@" > /dev/null 2>&1; if [ ! $? -eq 0 ]; then echo "Can't find file: $@"; exit 1; fi
    for i in "$@"; do
	clip $i | awk -v abs=$ABS -v avg=$AVG '{ total += ((abs==1) ? sqrt(('$WEIGHT'*$'$COMPONENT')^2) : '$WEIGHT'*$'$COMPONENT');	      count++ } END { print (avg==1) ? total/count : total }' - >> $TMPFILE || exit $?
    done
}
function length_vector() {
    ls -1 "$@" > /dev/null 2>&1; if [ ! $? -eq 0 ]; then echo "Can't find file: $@"; exit 1; fi
    for i in "$@"; do
	clip $i | awk -v sqr=$ABS -v avg=$AVG '{ sum=0; for (k=1; k<='$DIM'; k++) sum+=('$WEIGHT'*$k)^2; total+=((sqr==1) ? sum : sqrt(sum)); count++ } END { print (avg==1) ? total/count : total }' - >> $TMPFILE || exit $?
    done
}

if [ $# -eq 0 ]; then usage; exit $?; fi

TMPDIR="/tmp/"
RANDOMSTR=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32`
OUTFILE="/tmp/plot.sh.pdf"
TMPFILE="$TMPDIR/$(basename $OUTFILE).$RANDOMSTR.tmp"

PTYPE_FILE=""
SELECT_PTYPE=2

COMPONENT=1
ABS=0
AVG=1
DIM=1
SHOW=0
WEIGHT=1
VALUE_TYPE="Average"
TITLE=''
LABEL=''
while getopts ":ac:d:g:hl:o:p:stv:w:" optname; do
    case "$optname" in
	"a") ABS=1 ;;
	"c") COMPONENT=$OPTARG ;;
	"d") PTYPE_FILE=$OPTARG ;;
	"g") TITLE="$OPTARG" ;;
	"h") usage; exit $? ;;
	"l") LABEL="$OPTARG" ;;
	"o") OUTFILE="$OPTARG"; TMPFILE="$TMPDIR/$(basename $OUTFILE).$RANDOMSTR.tmp" ;;
	"p") SELECT_PTYPE=$OPTARG ;;
	"s") SHOW=1 ;;
	"t") VALUE_TYPE="Total"; AVG=0 ;;
	"v") DIM=$OPTARG ;;
	"w") WEIGHT=$OPTARG ;;
	*)
	    echo "Unknown option: -$OPTARG"
	    exit $?
	    ;;
    esac
done
shift $(($OPTIND - 1))

echo -n > $TMPFILE

if [ $DIM -eq 1 ]; then
    echo -n "$VALUE_TYPE scalar value (Scaling: $WEIGHT) ... "		&& [ $ABS -eq 1 ] && echo -n "(magnitude) "	&& VALUE_TYPE="$VALUE_TYPE magnitude"
    select_scalar "$@"
else
    echo -n "$VALUE_TYPE ${DIM}D vector length (Scaling: $WEIGHT) ... " && [ $ABS -eq 1 ] && echo -n "(length square) " && VALUE_TYPE="$VALUE_TYPE length square"
    length_vector "$@"
fi

EXT="${OUTFILE##*.}"
if [ $EXT == "m" ] || [ $EXT == "M" ]; then
    VNAME=${LABEL//[[:blank:]]/}
    echo "d_$VNAME = [" > "$OUTFILE"
    awk '{ print $1";" }' $TMPFILE >> "$OUTFILE"
    echo "];" >> "$OUTFILE"
elif [ $EXT == "txt" ] || [ $EXT == "TXT" ]; then
    mv $TMPFILE $OUTFILE
elif [ $EXT == "pdf" ] || [ $EXT == "pdf" ]; then
    gnuplot -persist <<PLOT
unset key
set terminal pdf font "Times,10"
set output "$OUTFILE"
set title "$TITLE"
set xlabel "Frame"
set ylabel "$VALUE_TYPE ($LABEL)"
plot '$TMPFILE' with lines
quit
PLOT

    if [ "$?" = "0" ] && [ $SHOW -eq 1 ]; then
	PDFVIEWER=${PDFVIEW:=`which xdg-open || which evince || which qpdfview`}
	if [ "$PDFVIEWER" == "" ]; then
	    echo "Cannot find any PDF viewer set using PDFVIEW"
	else
	    $PDFVIEWER "$OUTFILE"
	fi
    fi
fi

rm -f "$TMPFILE"
echo "saved to" "$OUTFILE"
