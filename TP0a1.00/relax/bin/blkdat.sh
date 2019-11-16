#!/bin/bash
# KAL - get X from input
if [ $# -ne 1 ] ; then
   echo "You must input experiment name (ex 01.0)"
   exit
fi

# Check that environment is set
if [ -z ${BASEDIR} ] ; then
   echo "Environment not set "
   exit
fi

set -x
X=$1
#
# --- create the relaxation blkdat.input from the expt version.
#
#
echo "Levitus (NOAA World Ocean Atlas 1994) Climatology"              > blkdat.input
echo "  00	  'month ' = month of climatology (01 to 12)"            >> blkdat.input
egrep "'iversn'|'iexpt '|'yrflag'" ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'idm   '|'jdm   '"          ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'itest '|'jtest '|'kdm   '" ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'nhybrd'|'nsigma'"          ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'dp00. '|'ds00. '"   ${BASEDIR}/expt_${X}/blkdat.input | egrep -v "'dp00i '"    >> blkdat.input
egrep "'thflag'|'thbase'|'sigma '" ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'thkmin'"                   ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
