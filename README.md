SpaceHASTEN version 0.1
Developed by Tuomo Kalliokoski, Orion Pharma <tuomo.kalliokoski at orionpharma.com>
2024-09-10

# Introduction

This is a tool that allows you screen BiosolveIT's .space-files using molecular docking tool Glide.
You can get large number of good scoring compounds from those vast chemical spaces of billions of compounds just by docking few million structures.
Hardware requirements: few hundred CPU cores, 1 reasonable GPU (2024) and few hundred gigabytes of disk space.
Only Linux is supported (tested on Ubuntu 22.04.4 and Rocky Linux 8.8).

SpaceHASTEN requires following commercial software:

* [Schrödinger Suite (version 2023-4): phase/ligprep, glide, python API (run)][https://www.schrodinger.com/release-download/]
* [SpaceLight (version 1.3.0)][https://www.biosolveit.de/download/?product=spacelight]
* [FTrees (version 6.11.0)][https://www.biosolveit.de/download/?product=ftrees]

Depending where you work, anaconda3 [may or may not be free for you][https://www.anaconda.com/pricing]:
* [miniconda3][https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh]

In addition, these free tools must be installed:

* [slurm workload manager (21.08.5 and 23.02.4 tested)][https://https://slurm.schedmd.com/documentation.html]
* [chemprop (version 1.7.1)][https://github.com/chemprop/chemprop/archive/refs/tags/v1.7.1.tar.gz]

Following commands in Sep 2024 can be used to install chemprop. This assumes that you have miniconda/anacond installed with mamba,
modify the included conda_activation_example.sh`(the same script required when installing SpaceHASTEN as well):

```
source conda_activation_example.sh
tar xzf v1.7.1.tar.gz
cd chemprop-1.7.1
conda create -y -n chemprop python=3.8
conda activate chemprop
mamba install -y pytorch=2.3.1 torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
mamba env update -f environment.yml
pip install -e .
```

In addition to software, you'll need some chemical spaces to search.
[Download them from BiosolveIT][https://www.biosolveit.de/chemical-spaces/]

Enamine provides diverse sets of Enamine REAL compounds that are good sources of seed molecules for Enamine REALSpace searches.
[Download Enamine REAL lead-like subset][https://enamine.net/compound-collections/real-compounds/real-database-subsets]

For more information, please see the manuscript describing how the method works. FIXME: link to manuscript.

# Installation and useage

Please check the requirements once more that make sure that you have all needed software and hardware.
If you still think that you have all pieces in place, follow these instructions:

* Check out this repository to some temporary directory: `git clone https://github.com/TuomoKalliokoski/SpaceHASTEN`
* Go to that directory and run the installer: `cd SpaceHASTEN ; python3 install_spacehasten.py`
* Start verify script as suggested by the installer to see that everything is running smoothly.
