Equilibrium and nonequilibrium simulation with cw pulse
====================================

Briefly, the goal is to get the average vibrational energy dynamics in each symmetry coordinate. In practice, starting from the thermally equilibrated geometry, one 5-ps NVE trajectory is simulated. And we also perform an additional nonequilibrium simulation to study the polariton relaxation and energy transfer dynamics under the NVE ensemble. This example runs on a periodic box of 400 methane molecules at 110 K. This example does not produce useful data, to run a converged simulation similar to that of example.png, forty 20-ps NVE trajectories both equilibrium and nonequilibrium are necessary.

Requirement
-----------------

This simulation requires the use of [LAMMPS](https://www.lammps.org/) and [cavity-md-ipi](https://github.com/TaoELi/cavity-md-ipi).

Usage
-----

**Equilibrium part:**

The line 1-2 in equilibrate_cavmd.sh set the number of methane molecules $N_{\rm{simu}}=400$ and effective light-matter coupling strength $\tilde{\varepsilon} = 6 \times 10 ^{-4} / \sqrt{N_{\rm{simu}}/100} = 3 \times 10 ^{-4}$ a.u. 

```bash
nmol=400
E0=6e-4
```

Run equilibrate_cavmd.sh script to start equilibrium NVE simulation.

```bash
sh equilibrate_cavmd.sh
```

This takes about 15 minutes to finish. 

**Nonequilibrium part**:

The line 1-4 in nonequilibrate_cavmd.sh set the cw pulse amplitude $E_0 = 6 \times 10 ^{-3}$ a.u., the number of methane molecules $N_{\rm{simu}}=400$ , effective light-matter coupling strength $\tilde{\varepsilon} = 6 \times 10 ^{-4} / \sqrt{N_{\rm{simu}}/100} = 3 \times 10 ^{-4}$ a.u. and traj-th simulation is performed.

```bash
amp=6e-3
nmol=400
E0=6e-4
traj=1
```

Get nonequilibrium NVE trajectories using the provided nonequilibrate_cavmd.sh script 

```bash
sh nonequilibrate_cavmd.sh
```

This takes about 15 minutes to finish.

**Dynamics in each symmetry coordinate**:

Run the analysis script to get the FFT result for symmetry coordinates and square symmetry coordinate trajectories $[v_{\sigma\lambda}(t)]^2$.

```bash
sh analysis.sh
```

This takes about 1 minute to finish.

**Plotting**:

In plot.py file, data1 is nonequilibrium trajectory and data2 is equilibrium trajectory. In /demos/data/ folder, the averaged 40 trajectories results (in coord_40 and noneqcoord_40) are also provided.

```python
data1 = np.loadtxt(f'./400/E0_6e-4/noneqcoord/simu_noneq_6e-3_1.xc.xyz.scoord.txt')[:2501,:]
data2 = np.loadtxt(f'./400/E0_6e-4/coord/simu_1.xc.xyz.scoord.txt')[:2501,:]
```

Run the python script plot.py and get output figure.png

```bash
python plot.py
```
