#!/usr/bin/env bash

source ~/.bashrc
export TOKENPATH=/home/anewman/tokens/
BASE="/home/anewman/GPS/Kivu"
PROC="$BASE/KivuGNSS.processing"
BIN="$PROC/bin"
DATA="$BASE/KivuGNSS"
# set to run on linux network at GT
pushd $DATA  # doesn't seem to work with ~
   # above processing within a single script now
   conda run -n Maps $BIN/get_rates2.py
   #plot vector maps
   conda run -n Maps $BIN/Kivu_GNSS_Map.py
   conda run -n Maps $BIN/Kivu_GNSS_Map_common.py
   # overlay images on map   
   $BIN/overlays.sh
   # create small copies of files
   $BIN/small_plots.sh plots/*png
   # organize plots
   \mv plots/*_TS.png plots/TS/
   \mv plots/*_TS_sm.png plots/TS/small
   \mv plots/*_TS_Common.png plots/TS/Common/
   \mv plots/*_TS_Common_sm.png plots/TS/Common/small
   \mv plots/*UKF.png plots/UKF/
   \mv plots/*UKF_sm.png plots/UKF/small

   # send it all to github
   $BIN/gitpush.sh
   # cp files to webdir
   cp -a ./plots/* ~/html/research/KivuGNSS/plots 
popd
