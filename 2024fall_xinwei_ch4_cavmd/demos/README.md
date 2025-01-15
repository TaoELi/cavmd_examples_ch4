Equilibrium and nonequilibrium simulation with cw pulse
====================================

Overall, the goal of this demo is to obtain the average vibrational energy dynamics per molecule in each symmetry coordinate. In practice, we first run a 150-ps NVT simulation to reach equilibrium. Starting from the thermally equilibrated geometry, we simulate one 20-ps NVE trajectory. Then, we perform an additional nonequilibrium simulation to study the polariton relaxation and energy transfer dynamics under the NVE ensemble. A total of 400 methane molecules are expliclity simulated at 110 K under periodic boundary conditions. Because only a single NVE trajectory is run, this demo calculation does not produce the converged data. To run a converged simulation similar to that of **example.png** (as shown in the paper), 40 20-ps NVE trajectories are required.

Requirement
-----------------

This simulation requires the use of [LAMMPS](https://www.lammps.org/) and [cavity-md-ipi](https://github.com/TaoELi/cavity-md-ipi).

Usage
-----

**Equilibrium part:**

Lines 1-2 in equilibrate_cavmd.sh set the number of methane molecules $N_{\rm{simu}}=400$ and the effective light-matter coupling strength $\tilde{\varepsilon} = 6 \times 10 ^{-4} / \sqrt{N_{\rm{simu}}/100} = 3 \times 10 ^{-4}$ a.u. 

```bash
nmol=400
E0=6e-4
```

Change the /path/to in line 4 in equilibrate_cavmd.sh to the place you want to run this simulation.

```bash
TOTALFOLDER=/path/to/demos/"$nmol"
```

Run NVT simulations to reach equilibrium geometry and get equilibrium NVE trajectories using the provided equilibrate_cavmd.sh script.

```bash
sh equilibrate_cavmd.sh
```

This takes about 6 hours and a half to finish. Final thermally equilibrared geometry init_0.checkpoint and other 40 NVE final geometries are provided in folder /data/checkpoint_400_E0_6e-4/. 

**To skip the long NVT simulation**, create a new folder /path/to/demos/400/E0_6e-4/checkpoint, copy the init_0.checkpoint file to the new folder and rerun the equilibrate_cavmd.sh script. The script will automatically skip the long NVT simulation. The equilibrium NVE simulation will take about 1 hour to finish.

```bash
mkdir /path/to/demos/400/E0_6e-4/checkpoint
cp /data/checkpoint_400_E0_6e-4/init_0.checkpoint /path/to/demos/400/E0_6e-4/checkpoint
sh equilibrate_cavmd.sh
```

**Nonequilibrium part**:

Lines 1-4 in nonequilibrate_cavmd.sh set the cw pulse amplitude $E_0 = 6 \times 10 ^{-3}$ a.u., the number of methane molecules $N_{\rm{simu}}=400$ , the effective light-matter coupling strength $\tilde{\varepsilon} = 6 \times 10 ^{-4} / \sqrt{N_{\rm{simu}}/100} = 3 \times 10 ^{-4}$ a.u. and the **$traj**-th simulation is performed.

```bash
amp=6e-3
nmol=400
E0=6e-4
traj=1
```

Change the /path/to in line 7 in nonequilibrate_cavmd.sh to the place you want to run this simulation.

```bash
HOMEFOLDER=/path/to/demos/"$nmol"/E0_$E0
```

Get nonequilibrium NVE trajectories using the provided nonequilibrate_cavmd.sh script 

```bash
sh nonequilibrate_cavmd.sh
```

This takes about 1 hour to finish.

**Dynamics in each symmetry coordinate**:

Change the /path/to in lines 3-4 in analysis.sh to the place you ran this simulation and the foldername where the NVE trajectories are stored either equilibrium or nonequilibrium. To get average vibrational energy dynamics in each symmetry coordinate, you need to run this script both the equilibrium and nonequilibrium.

```bash
OTALFOLDER=/path/to
XYZFOLDER=$TOTALFOLDER/xyznve
```

You can also change line 7 and line 13 the change the folders where the FFT result for symmetry coordinates and square symmetry coordinate trajectories $[v_{\sigma\lambda}(t)]^2$ are.

```bash
PACFOLDER=$TOTALFOLDER/pac
COORDFOLDER=$TOTALFOLDER/coord
```

Run the analysis script

```bash
sh analysis.sh
```

This takes about 1 minute to finish.

**Plotting**:

Change the lines 52-53 to your square symmetry coordinate trajectories path, data1 is nonequilibrium trajectory and data2 is equilibrium trajectory. In /demos/data/ folder, Both 1 trajectory results (in coord_1 and noneqcoord_1) and averaged 40 trajectories results (in coord_40 and noneqcoord_40) are provided.
```python
data1 = np.loadtxt(f'./data/noneqcoord_1/simu_noneq_6e-3_1.xc.xyz.scoord.txt')
data2 = np.loadtxt(f'./data/coord_1/simu_1.xc.xyz.scoord.txt')
```
Run the python script plot.py

```bash
python plot.py
```
