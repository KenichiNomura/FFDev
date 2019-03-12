#!/bin/sh

#if [ "x$QC" == "x" ]; then
#   t=`dirname  $0`
#   export QC=`(cd $t >/dev/null 2>&1; pwd -P)`
#fi
export QC=/usr/usc/qchem/5.1.1
#echo "QC = $QC"

function qc { 
   cd $QC
}
function get_shellvar { 
   TT="${1} "; grep "${TT}" $QC/config/shellvar.txt | tail -n 1 | awk '{print $2}'
}
function add_path { 
   TT="$1"; export PATH=`echo ${TT}:$PATH| sed -e "s|:${TT}:|:|g"` 
}

if [ "x${QCPLATFORM}" == "x" ]; then
   export QCPLATFORM=$(get_shellvar  QCPLATFORM)
fi
#echo "QCPLATFORM = $QCPLATFORM"
if [ "x${QCSCRATCH}" == "x" ]; then
   export QCSCRATCH=$SCRATCHDIR
fi
if [ "x${QCLOCALSCR}" == "x" ]; then
   export QCLOCALSCR=$(get_shellvar QCLOCALSCR)
fi
if [ "x${QCLOCALSCR}" == "x" ]; then
   unset QCLOCALSCR
fi
if [ "x${QCAUX}" == "x" ]; then
   export QCAUX=$(get_shellvar QCAUX)
echo $QCAUX
fi
if [ "x${QCMPI}" == "x" ]; then
   export QCMPI=$(get_shellvar QCMPI)
fi
if [ "x${QCRSH}" == "x" ]; then
   export QCRSH=$(get_shellvar QCRSH)
fi

export QCPROG=$QC/exe/qcprog.exe
if [ -e $QC/exe/qcprog.exe_s ] ; then
   export QCPROG_S=$QC/exe/qcprog.exe_s
else
   export QCPROG_S=$QCPROG
fi
add_path "${QC}/bin"
add_path "${QC}/exe"
#echo "$PATH"

if [ ! -e $QCAUX/license/qchem.license.dat ]; then 
     echo Cannot find license file: $QCAUX/license/qchem.license.dat not found 
fi

