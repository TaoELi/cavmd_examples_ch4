#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
##SBATCH --mem-per-cpu=1024M
#SBATCH --job-name=cavityloss
#SBATCH --partition=standard
#SBATCH --time=1-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

amp=6e-3
nmol=400
E0=6e-4
traj=1
FREQIDX=0

echo "running nonequilibrate NVE for nmol=$nmol, freqidx=$FREQIDX, traj=$traj"

TOTALFOLDER=/path/to/fix_frequencies
HOMEFOLDER=$TOTALFOLDER/freqidx_"$FREQIDX"

if [ ! -d $HOMEFOLDER ]; then
    mkdir $HOMEFOLDER
fi

ORIGINFOLDER=$TOTALFOLDER/E0_$E0/running_input_$E0
CHECKPOINTFOLDER=$TOTALFOLDER/E0_$E0/checkpoint
INPUTFOLDER=$HOMEFOLDER/noneq_cavity_loss_running_input_"$FREQIDX"
 
if [ "$traj" = 1 ]; then
    echo "generating new input files"
else
    echo "sleep 30s to make sure first job generate all input files"
    sleep 30s
fi

if [ ! -d $INPUTFOLDER ]; then

    mkdir $INPUTFOLDER
    cp $ORIGINFOLDER/* $INPUTFOLDER
    cd $INPUTFOLDER

    if [ -f 'photon_params.json' ]; then

        echo "remove useless photon_params.json"    
        rm photon_params.json

    fi

    sed -i "s/'''/ /g" $INPUTFOLDER/create_photon_params.py
    sed -i "s/freqidx=0/freqidx=$FREQIDX/" $INPUTFOLDER/create_photon_params.py
    sed -i "s/\[0.0006/\[$amp/" $INPUTFOLDER/create_photon_params.py
    python -u $INPUTFOLDER/create_photon_params.py
    sleep 10s

else

    echo "Inputfolder for E0=$E0 exists"
    sleep 10s

fi

cd $HOMEFOLDER

# run 1 noneq NVE simulation
NVEXYZ=$HOMEFOLDER/noneq_gaussian_xyznve_"$FREQIDX"
mkdir $NVEXYZ
OUTFOLDER=$HOMEFOLDER/noneq_gaussian_out_"$FREQIDX"
mkdir $OUTFOLDER

NVEFOLDER=$HOMEFOLDER/noneq_gaussian_nve_"$FREQIDX"_"$traj"

if [ ! -d $NVEFOLDER ]; then

    mkdir $NVEFOLDER
    cp $INPUTFOLDER/* $NVEFOLDER
    cd $NVEFOLDER
    sed -i "s/mesitylene-pimd.1/ch4.E0_$E0.cavity_loss_traj_$traj/" input_traj.xml
    sed -i "s/mesitylene-pimd.1/ch4.E0_$E0.cavity_loss_traj_$traj/" in.lmp
    sed -i "s/RESTART/init_$(($traj-1)).checkpoint/" input_traj.xml
    sed -i "s/'simu'/'simu_noneq_"$amp"_$traj'/g" input_traj.xml
    sed -i "s/500/$(($nmol*5))/g" input_traj.xml
    sed -i "s/501/$(($nmol*5+1))/g" input_traj.xml

fi

NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_$(($traj-1)).checkpoint
cp $NEWNVECHECKPOINTFILE $NVEFOLDER
cd $NVEFOLDER

rm -rf /tmp/ipi_ch4.E0_$E0.cavity_loss_traj_$traj
i-pi input_traj.xml & > output.output &
sleep 60s
lmp < in.lmp
sleep 60s
            
cd $HOMEFOLDER
mv $NVEFOLDER/simu_noneq_"$amp"_$traj.xc.xyz $NVEXYZ
mv $NVEFOLDER/simu_noneq_"$amp"_$traj.out $OUTFOLDER
