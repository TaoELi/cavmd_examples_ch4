#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
##SBATCH --mem-per-cpu=1024M
#SBATCH --job-name=test
#SBATCH --partition=taoeli
#SBATCH --time=7-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

nmol=400
E0=0e-4
name=1
TOTALFOLDER=/path/to/density_dependence/n400v$name

if [ ! -d $TOTALFOLDER ]; then
    mkdir $TOTALFOLDER
fi

HOMEFOLDER=$TOTALFOLDER/E0_$E0
if [ ! -d $HOMEFOLDER ]; then
    mkdir $HOMEFOLDER
fi

ORIGINFOLDER=$TOTALFOLDER/origin_input
INPUTFOLDER=$HOMEFOLDER/running_input_$E0
mkdir $INPUTFOLDER
cp $ORIGINFOLDER/* $INPUTFOLDER
cd $INPUTFOLDER
newE0=$(echo "print($E0*($nmol//100)**(-0.5))" | python)
sed -i "s/0.0000/$newE0/" $INPUTFOLDER/create_photon_params.py
python $INPUTFOLDER/create_photon_params.py

# run NVT equilibrate
CHECKPOINTFOLDER=$HOMEFOLDER/checkpoint
if [ ! -d $CHECKPOINTFOLDER ]; then
    mkdir $CHECKPOINTFOLDER
fi

NVTPATTERN="<step>300000</step>"
NVTFOLDER=$HOMEFOLDER/nvt

if [ -f "$CHECKPOINTFOLDER/init_0.checkpoint" ]; then
    
    if grep -q $NVTPATTERN "$CHECKPOINTFOLDER/init_0.checkpoint"; then
        echo 'NVT equilibrate for 150 ps has finished'
    fi
    
else

    if [ ! -d "$NVTFOLDER" ]; then

        mkdir $NVTFOLDER
        cp $INPUTFOLDER/* $NVTFOLDER
        cd $NVTFOLDER
        sed -i "s/mesitylene-pimd.1/ch4."$name"."$nmol".E0_"$E0".eq/" input_eq.xml
        sed -i "s/mesitylene-pimd.1/ch4."$name"."$nmol".E0_"$E0".eq/" in.lmp

    fi

    cd $NVTFOLDER
    echo 'NVT checkpoint file does not exist, run NVT equilibrate 150 ps'
    rm -rf /tmp/ipi_ch4."$name"."$nmol".E0_"$E0".eq
    i-pi input_eq.xml & > nvteq.output &
    sleep 60s
    lmp < in.lmp
    sleep 60s
            
    cd $HOMEFOLDER
    mv $NVTFOLDER/init_0.checkpoint $CHECKPOINTFOLDER

fi

# run 40 NVE simulation
NVEPATTERN="<step>40000</step>"
NVEXYZ=$HOMEFOLDER/xyznve
mkdir $NVEXYZ
OUTFOLDER=$HOMEFOLDER/out
mkdir $OUTFOLDER

for traj in {1..40}
do  
    
    NVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_$traj.checkpoint
    NVEFOLDER=$HOMEFOLDER/nve_$traj
    
    if [ -f "$NVECHECKPOINTFILE" ]; then
    
        if grep -q $NVEPATTERN "$NVECHECKPOINTFILE"; then
    	      echo "Found checkpoint for $traj-th NVE 20 ps simulation finished, skip $traj-th sequential job"
        fi

    else

        if [ ! -d "$NVEFOLDER" ]; then

            mkdir $NVEFOLDER
            cp $INPUTFOLDER/* $NVEFOLDER
            cd $NVEFOLDER
            sed -i "s/mesitylene-pimd.1/ch4."$name"."$nmol".E0_"$E0".traj_"$traj"/" input_traj.xml
            sed -i "s/mesitylene-pimd.1/ch4."$name"."$nmol".E0_"$E0".traj_"$traj"/" in.lmp
            sed -i "s/RESTART/init_$(($traj-1)).checkpoint/" input_traj.xml
            sed -i "s/'simu'/'simu_$traj'/g" input_traj.xml
            sed -i "s/500/$(($nmol*5))/g" input_traj.xml
            sed -i "s/501/$(($nmol*5+1))/g" input_traj.xml
                        
        fi
        
        NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_$(($traj-1)).checkpoint
        cp $NEWNVECHECKPOINTFILE $NVEFOLDER
        cd $NVEFOLDER
        echo "$traj-th NVE checkpoint file does not exist, run NVE simulation 20 ps"
        rm -rf /tmp/ipi_ch4."$name"."$nmol".E0_"$E0".traj_"$traj"
        i-pi input_traj.xml & > output.output &
        sleep 60s
        lmp < in.lmp
        sleep 60s
                        
        cp simu_$traj.checkpoint init_$traj.checkpoint
        cd $HOMEFOLDER
        mv $NVEFOLDER/init_$traj.checkpoint $CHECKPOINTFOLDER
        mv $NVEFOLDER/simu_$traj.xc.xyz $NVEXYZ
        mv $NVEFOLDER/simu_$traj.out $OUTFOLDER    
    
    fi
    
done
