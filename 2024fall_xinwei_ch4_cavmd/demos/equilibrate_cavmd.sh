nmol=400
E0=6e-4
name=ch4
CURRENTFOLDER=$(pwd)
TOTALFOLDER=$CURRENTFOLDER/"$nmol"

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

# run 1 NVE simulation
NVEPATTERN="<step>10000</step>"
NVEXYZ=$HOMEFOLDER/xyznve
mkdir $NVEXYZ
OUTFOLDER=$HOMEFOLDER/out
mkdir $OUTFOLDER

CHECKPOINTFOLDER=$TOTALFOLDER/checkpoint
NVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_1.checkpoint
NVEFOLDER=$HOMEFOLDER/nve_1

if [ -f "$NVECHECKPOINTFILE" ]; then

    if grep -q $NVEPATTERN "$NVECHECKPOINTFILE"; then
            echo "Found checkpoint for 1-th NVE 5 ps simulation finished, skip 1-th sequential job"
    fi

else

    if [ ! -d "$NVEFOLDER" ]; then

        mkdir $NVEFOLDER
        cp $INPUTFOLDER/* $NVEFOLDER
        cd $NVEFOLDER
        sed -i "s/mesitylene-pimd.1/ch4."$name"."$nmol".E0_"$E0".traj_"1"/" input_traj.xml
        sed -i "s/mesitylene-pimd.1/ch4."$name"."$nmol".E0_"$E0".traj_"1"/" in.lmp
        sed -i "s/RESTART/init_$((1-1)).checkpoint/" input_traj.xml
        sed -i "s/'simu'/'simu_1'/g" input_traj.xml
        sed -i "s/500/$(($nmol*5))/g" input_traj.xml
        sed -i "s/501/$(($nmol*5+1))/g" input_traj.xml
                    
    fi
    
    NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_$((1-1)).checkpoint
    cp $NEWNVECHECKPOINTFILE $NVEFOLDER
    cd $NVEFOLDER
    echo "1-th NVE checkpoint file does not exist, run NVE simulation 20 ps"
    rm -rf /tmp/ipi_ch4."$name"."$nmol".E0_"$E0".traj_"1"
    i-pi input_traj.xml & > output.output &
    sleep 10s
    lmp < in.lmp
    sleep 10s
                    
    cp simu_1.checkpoint init_1.checkpoint
    cd $HOMEFOLDER
    mv $NVEFOLDER/init_1.checkpoint $CHECKPOINTFOLDER
    mv $NVEFOLDER/simu_1.xc.xyz $NVEXYZ
    mv $NVEFOLDER/simu_1.out $OUTFOLDER    

fi
