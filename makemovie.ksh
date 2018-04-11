#!/bin/ksh

exps="NATURE NODA"
fhr=0
fhrmax=120
fhrmax=120
typeset -Z3 fhr

cd /gpfs/hps3/emc/meso/save/Donald.E.Lippi/pygraphics/figs/201803280



while [[ $fhr -le $fhrmax ]]; do
    echo $fhr
    toppng=NATURE_REFC_CONUS_f${fhr}_CONUSNEST.png
    botpng=NODA_REFC_CONUS_f${fhr}_CONUSNEST.png
    montage -tile 1x2 -geometry +4+4 $toppng $botpng ./movie/REFC_CONUS_f${fhr}_CONUSNEST.png
    ((fhr=$fhr+3))
done
    cd movie
    convert -delay 100 -loop 0 *.png REFC_CONUS_movie.gif

