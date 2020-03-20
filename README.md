# mutect2-nf

GATK4 mutect2 test for missing variants in NGS composite run

# setup

Clone this repository from GitHub:

```
git clone --recursive https://github.com/varshini712/mutect2-nf.git
cd mutect2-nf
```
# usage

Make sure you have your container image files and reference files created/specified in your working dir.
Use Makefile in this repo to submit job to slurm

Installs nextflow
```
make install
```
Submits sbatch job to slurm
```
make submit
```

