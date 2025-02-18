#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --job-name=nve
#SBATCH --partition=taoeli
#SBATCH --time=7-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

amp=6e-3 # In equilibrium simulation this parameter does not work
nmol=100
E0=0e-4
traj=1
echo "running nonequilibrate NVE for nmol=$nmol, E0=$E0, traj=$traj"

HOMEFOLDER=/path/to/ml_100_ch4/E0_$E0
ORIGINFOLDER=$HOMEFOLDER/running_input_$E0
CHECKPOINTFOLDER=$HOMEFOLDER/checkpoint
INPUTFOLDER=$HOMEFOLDER/nve_running_input_"$E0"
MLFOLDER=/path/to/ml_100_ch4/ml_potential

if [ ! -d $INPUTFOLDER ]; then

    mkdir $INPUTFOLDER
    cp $ORIGINFOLDER/* $INPUTFOLDER
    cp $MLFOLDER/* $INPUTFOLDER
    cd $INPUTFOLDER
    
    sleep 10s

    if [ -f 'photon_params.json' ]; then

        echo "remove useless photon_params.json"    
        rm photon_params.json

    fi

    python -u $INPUTFOLDER/create_photon_params.py
    sleep 10s

else

    echo "Inputfolder for E0=$E0 exists"
    sleep 10s

fi

cd $HOMEFOLDER

# run 1 noneq NVE simulation
NVEFOLDER=$HOMEFOLDER/nve_"$traj"
XYZFOLDER=$HOMEFOLDER/xyznve
OUTFOLDER=$HOMEFOLDER/out
if [ ! -d $XYZFOLDER ]; then
    mkdir $XYZFOLDER
fi
if [ ! -d $OUTFOLDER ]; then
    mkdir $OUTFOLDER
fi
if [ ! -d $NVEFOLDER ]; then

    mkdir $NVEFOLDER
    cp $INPUTFOLDER/* $NVEFOLDER
    sleep 10s
    cd $NVEFOLDER
    sed -i "s/mesitylene-pimd.1/ch4.E0_$E0.traj_noneq_"$amp"_$traj/" input_traj.xml
    sed -i "s/mesitylene-pimd.1/ch4.E0_$E0.traj_noneq_"$amp"_$traj/" in.lmp
    sed -i "s/RESTART/simu_nve_"$traj".checkpoint/" input_traj.xml
    sed -i "s/'simu'/'simu_nve_6e-4_$traj'/g" input_traj.xml

fi

NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/simu_nve_"$traj".checkpoint
cp $NEWNVECHECKPOINTFILE $NVEFOLDER
cd $NVEFOLDER

rm -rf /tmp/ipi_ch4.E0_$E0.traj_noneq_"$amp"_$traj
i-pi input_traj.xml & > output.output &
sleep 60s
lmp < in.lmp
sleep 60s
            
mv $NVEFOLDER/simu_nve_6e-4_$traj.xc.xyz $XYZFOLDER
mv $NVEFOLDER/simu_nve_6e-4_$traj.out $OUTFOLDER 
