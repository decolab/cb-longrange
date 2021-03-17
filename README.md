# Long_Range
The code in this directory is the code that was used for Deco et al. ms
“Rare long-range cortical connections enhance information processing”.

1. System requirements
The code in this directory works for Matlab2018b and newer. Due to the
massive computations for computing information processing for brain
imaging data, we have included the SLURM files to be uploaded and run on
the cluster.

2. Installation guide
Move all of the matlab files and dependencies into a folder and then use
the two main programs:
1. Ring (and demo program):
run_demo_ring_hopf_G.m (

2. fMRI whole-brain model:
slurm.sbatch_hopf_turbu_longrange.m (to be run on a SLURM compatible
cluster)
read_hopf_SClong.m (to read the output of the cluster results)

3. Demo
Load “run_demo_ring_hopf_G.m” into Matlab
Run this program
Depending the G_range (in line 55) and the computing power of computer,
this will take a couple of hours.
The expected output is all the variables shown in the ms. As an example,
the code will plot the turbulence as a function global coupling (G) for
short range and long range connections.

4. Instructions for use
You will need a license for Matlab to run the ring results. Due to the
computational demands for the neuroimaging results, you will need a
SLURM compatible cluster.
The two main programs will generate and reproduce all of the
quantitative results in the manuscript.
