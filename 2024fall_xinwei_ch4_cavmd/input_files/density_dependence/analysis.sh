#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=4G
##SBATCH --mem-per-cpu=1024M
#SBATCH --job-name=pac
#SBATCH --partition=taoeli
#SBATCH --time=7-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

nmol=400
e0=0
amp=6e-3

ANALYSISFILE=/path/to/mdanalysis.py

WASTEFOLDER=/path/to/waste
if [ ! -d $WASTEFOLDER ]; then
    mkdir $WASTEFOLDER
fi

TOTALFOLDER=/path/to
XYZFOLDER=$TOTALFOLDER/xyznve
python $ANALYSISFILE $XYZFOLDER

DACFOLDER=$TOTALFOLDER/dac

if [ ! -d $DACFOLDER ]; then
    mkdir $DACFOLDER
fi

PACFOLDER=$TOTALFOLDER/pac

if [ ! -d $PACFOLDER ]; then
    mkdir $PACFOLDER
fi

COORDFOLDER=$TOTALFOLDER/coord

if [ ! -d $COORDFOLDER ]; then
    mkdir $COORDFOLDER
fi

AACFOLDER=$TOTALFOLDER/aacvst

if [ ! -d $AACFOLDER ]; then
    mkdir $AACFOLDER
fi

sleep 3s

mv $XYZFOLDER/simu_*.xc.xyz.dac.txt $DACFOLDER
mv $XYZFOLDER/simu_*.xc.xyz.pac.txt $PACFOLDER
mv $XYZFOLDER/simu_*.xc.xyz.scoord.txt $COORDFOLDER
mv $XYZFOLDER/simu_noneq_"$amp"_"$traj".xc.xyz.aac_*.txt $AACFOLDER
