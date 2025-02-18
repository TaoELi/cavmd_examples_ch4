#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --job-name=nvt
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
CHECKPOINTFOLDER=$HOMEFOLDER/checkpoint_old
NEWCHECKPOINTFOLDER=$HOMEFOLDER/checkpoint

if [ ! -d $NEWCHECKPOINTFOLDER ]; then
    mkdir $NEWCHECKPOINTFOLDER
fi

INPUTFOLDER=$HOMEFOLDER/nvt_running_input_"$E0"
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

# run 1 noneq NVT simulation
NVTFOLDER=$HOMEFOLDER/nvt_"$traj"

if [ ! -d $NVTFOLDER ]; then

    mkdir $NVTFOLDER
    cp $INPUTFOLDER/* $NVTFOLDER
    sleep 10s
    cd $NVTFOLDER
    sed -i "s/mesitylene-pimd.1/ch4.E0_$E0.traj_noneq_"$amp"_$traj/" input_eq.xml
    sed -i "s/mesitylene-pimd.1/ch4.E0_$E0.traj_noneq_"$amp"_$traj/" in.lmp
    sed -i "s/RESTART/init_$(($traj-1)).checkpoint/" input_eq.xml
    sed -i "s/'simu'/'simu_nvt_$traj'/g" input_eq.xml

fi

NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_$(($traj-1)).checkpoint
cp $NEWNVECHECKPOINTFILE $NVTFOLDER
cd $NVTFOLDER

rm -rf /tmp/ipi_ch4.E0_$E0.traj_noneq_"$amp"_$traj
i-pi input_eq.xml & > output.output &
sleep 60s
lmp < in.lmp
sleep 60s
            
cp simu_nvt_$traj.checkpoint init_nvt_$(($traj-1)).checkpoint
mv $NVTFOLDER/init_nvt_$(($traj-1)).checkpoint $NEWCHECKPOINTFOLDER
