#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --job-name=noneqnve
#SBATCH --partition=taoeli
#SBATCH --time=7-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

amp=6e-3
nmol=100
E0=0e-4
traj=1
echo "running nonequilibrate NVE for nmol=$nmol, E0=$E0, traj=$traj"

HOMEFOLDER=/path/to/ml_100_ch4/E0_$E0
ORIGINFOLDER=$HOMEFOLDER/running_input_$E0
CHECKPOINTFOLDER=$HOMEFOLDER/checkpoint
INPUTFOLDER=$HOMEFOLDER/noneqnve_running_input_"$E0"_"$amp"
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

    sed -i "s/'''/ /g" $INPUTFOLDER/create_photon_params.py
    sed -i "s/\[0.0006/\[$amp/" $INPUTFOLDER/create_photon_params.py
    python -u $INPUTFOLDER/create_photon_params.py
    sleep 10s

else

    sleep 20s
    echo "Inputfolder for E0=$E0 exists"

fi

cd $HOMEFOLDER

# run 1 noneq NVE simulation
NVEFOLDER=$HOMEFOLDER/noneqnve_"$traj"_"$amp"
XYZFOLDER=$HOMEFOLDER/noneqxyznve_"$amp"
OUTFOLDER=$HOMEFOLDER/noneqout_"$amp"
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
    sed -i "s/'simu'/'simu_noneqnve_6e-4_"$traj"_"$amp"'/g" input_traj.xml

fi

NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/simu_nve_"$traj".checkpoint
cp $NEWNVECHECKPOINTFILE $NVEFOLDER
cd $NVEFOLDER

rm -rf /tmp/ipi_ch4.E0_$E0.traj_noneq_"$amp"_$traj
i-pi input_traj.xml & > output.output &
sleep 60s
lmp < in.lmp
sleep 60s
            
mv $NVEFOLDER/simu_noneqnve_6e-4_"$traj"_"$amp".xc.xyz $XYZFOLDER
mv $NVEFOLDER/simu_noneqnve_6e-4_"$traj"_"$amp".out $OUTFOLDER 
