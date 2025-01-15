ANALYSISFILE=./mdanalysis.py

TOTALFOLDER=/path/to/demos/400/E0_6e-4
XYZFOLDER=$TOTALFOLDER/xyznve
python $ANALYSISFILE $XYZFOLDER

PACFOLDER=$TOTALFOLDER/pac

if [ ! -d $PACFOLDER ]; then
    mkdir $PACFOLDER
fi

COORDFOLDER=$TOTALFOLDER/coord

if [ ! -d $COORDFOLDER ]; then
    mkdir $COORDFOLDER
fi

sleep 3s

mv $XYZFOLDER/simu_*.xc.xyz.pac.txt $PACFOLDER
mv $XYZFOLDER/simu_*.xc.xyz.scoord.txt $COORDFOLDER
