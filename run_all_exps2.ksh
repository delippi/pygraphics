#!/bin/ksh

startTime=2018032800
exppref="fv3gfs_dl2rw_"
exps="fv3gfs_dl2rw_NODA fv3gfs_dl2rw_NATURE"
exps="fv3gfs_dl2rw_NATURE"
tmpdir='ptmp'
fhr=0
fhrmax=120

for exp in $exps; do
    while [[ $fhr -le $fhrmax ]] ; do
        echo $fhr
        python /gpfs/hps3/emc/meso/save/Donald.E.Lippi/pygraphics/threaded_namv4_2d.py $exp $fhr $startTime $tmpdir
        ((fhr=$fhr+3))
    done
done
