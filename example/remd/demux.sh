#!/bin/bash
set -e

PREFIX="grest"
REPLEX=500
SAVE_INTERVAL=500
MD_LOG="rep1/${PREFIX}.log"
EXTRA=$(echo $REPLEX $SAVE_INTERVAL | awk '{print $1 / $2 - 1}')

mdtbx cmd demux.pl $MD_LOG $EXTRA

mdtbx cmd gmx trjcat -f rep*/${PREFIX}.xtc -demux replica_index.xvg -o demux_trj.xtc

echo done
