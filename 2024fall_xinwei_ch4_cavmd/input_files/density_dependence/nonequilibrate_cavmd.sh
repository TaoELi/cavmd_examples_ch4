#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
##SBATCH --mem-per-cpu=1024M
#SBATCH --job-name=noneqnve
#SBATCH --partition=standard
#SBATCH --time=1-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

amp=6e-3
nmol=400
E0=0e-4
traj=1
name=1
echo "run single nonequilibrium NVE simulation for nmol=$nmol, E0=$E0, traj=$traj, amp=$amp, density=$name"

HOMEFOLDER=/path/to/density_dependence/n400v"$name"/E0_$E0
ORIGINFOLDER=$HOMEFOLDER/running_input_$E0
CHECKPOINTFOLDER=$HOMEFOLDER/checkpoint
INPUTFOLDER=$HOMEFOLDER/noneq_"$amp"_running_input_"$E0"

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
    sed -i "s/\[0.0006/\[$amp/" $INPUTFOLDER/create_photon_params.py
    OLDPATH=/work/taoeli/users/astolfo/project/ch4_cavity/origin_input/freq.txt
    NEWPATH=$INPUTFOLDER/freq.txt
    sed -i "s|$OLDPATH|$NEWPATH|" $INPUTFOLDER/create_photon_params.py
    sqrtN=$(echo "print(($nmol//100)**(-0.5))" | python)
    sed -i "s|abs(freq\[e0,0\]|abs(freq\[e0,0\] * $sqrtN|" $INPUTFOLDER/create_photon_params.py
    newE0=$(echo "print($E0*($nmol//100)**(-0.5))" | python)
    sed -i "s/nmol=100/nmol=$nmol/" $INPUTFOLDER/create_photon_params.py
    python -u $INPUTFOLDER/create_photon_params.py
    sleep 10s

else
    echo "Inputfolder for E0=$E0 exists"
    sleep 10s

fi

cd $HOMEFOLDER

# run 1 noneq NVE simulation
NVEXYZ=$HOMEFOLDER/noneqxyznve_"$amp"
mkdir $NVEXYZ
OUTFOLDER=$HOMEFOLDER/noneqout_"$amp"
mkdir $OUTFOLDER

NVEFOLDER=$HOMEFOLDER/eqnve_"$amp"_"$traj"

if [ ! -d $NVEFOLDER ]; then

    mkdir $NVEFOLDER
    cp $INPUTFOLDER/* $NVEFOLDER
    cd $NVEFOLDER
    sed -i "s/mesitylene-pimd.1/ch4.E0_$E0.traj_noneq_"$amp"_$traj/" input_traj.xml
    sed -i "s/mesitylene-pimd.1/ch4.E0_$E0.traj_noneq_"$amp"_$traj/" in.lmp
    sed -i "s/RESTART/init_$(($traj-1)).checkpoint/" input_traj.xml
    sed -i "s/'simu'/'simu_noneq_"$amp"_$traj'/g" input_traj.xml
    sed -i "s/500/$(($nmol*5))/g" input_traj.xml
    sed -i "s/501/$(($nmol*5+1))/g" input_traj.xml

fi

NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_$(($traj-1)).checkpoint
cp $NEWNVECHECKPOINTFILE $NVEFOLDER
cd $NVEFOLDER

rm -rf /tmp/ipi_ch4.E0_$E0.traj_noneq_"$amp"_$traj
i-pi input_traj.xml & > output.output &
sleep 60s
lmp < in.lmp
sleep 60s
            
cd $HOMEFOLDER
mv $NVEFOLDER/simu_noneq_"$amp"_$traj.xc.xyz $NVEXYZ
mv $NVEFOLDER/simu_noneq_"$amp"_$traj.out $OUTFOLDER