#!/bin/sh 
# 

MINPARAMS=1
if [ $# -lt "$MINPARAMS" ]
then
  echo
  echo "Usage: DoRecoDiff.sh 'CMSSW_X_Y_Z'"
  exit
fi  


CUR_DIR=`pwd`
cd /afs/cern.ch/cms/sw/slc4_ia32_gcc345/cms/cmssw/$1/src
eval `scramv1 runtime -sh`
cd $CUR_DIR
$CUR_DIR/ConfdbMakeOfflineConfig.py >& CMSSWold_cff.py
eval `scramv1 runtime -sh`
export THE_OLD_RELEASE_NAME=$1
$CUR_DIR/ConfdbDiffOfflineConfig.py >& thediff.twikiformat.txt
exit
