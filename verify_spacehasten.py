# SpaceHASTEN: check environment sanity
#
# Copyright (c) 2024 Orion Corporation
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software 
# without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
#
import os
import cfg
import time
import glob
#
# these are here to check that the required software is available
import pandas
import rdkit
import tqdm
import numpy

def write_confgen_file(c):
    w = open("test.inp","wt")
    w.write("[SET:ORIGINAL_LIGANDS]\n")
    w.write("    VARCLASS   Structures\n")
    w.write("    FILES   examples.smi,\n")
    w.write("\n")
    w.write("[STAGE:LIGPREP]\n")
    w.write("    STAGECLASS   ligprep.LigPrepStage\n")
    w.write("    INPUTS   ORIGINAL_LIGANDS,\n")
    w.write("    OUTPUTS   LIGPREP_OUT,\n")
    w.write("    RECOMBINE   YES\n")
    w.write("    RETITLE   YES\n")
    w.write("    MIXLIGS   YES\n")
    w.write("    SKIP_BAD_LIGANDS   YES\n")
    w.write("    UNIQUEFIELD   s_m_title\n")
    w.write("    OUTCOMPOUNDFIELD   s_m_title\n")
    w.write("    USE_EPIK   YES\n")
    w.write("    METAL_BINDING   NO\n")
    w.write("    PH   7.0\n")
    w.write("    PHT   2.0\n")
    w.write("    NRINGCONFS   1\n")
    w.write("    COMBINEOUTS   NO\n")
    w.write("    STEREO_SOURCE   parities\n")
    w.write("    NUM_STEREOISOMERS   32\n")
    w.write("    REGULARIZE   NO\n")
    w.write("\n")
    w.write("[STAGE:POSTLIGPREP]\n")
    w.write("    STAGECLASS   ligprep.PostLigPrepStage\n")
    w.write("    INPUTS   LIGPREP_OUT,\n")
    w.write("    OUTPUTS   POSTLIGPREP_OUT,\n")
    w.write("    UNIQUEFIELD   s_m_title\n")
    w.write("    OUTVARIANTFIELD   s_phase_variant\n")
    w.write("    PRESERVE_NJOBS   YES\n")
    w.write("    LIMIT_STEREOISOMERS   YES\n")
    w.write("    MAXSTEREO   4\n")
    w.write("    REMOVE_PENALIZED_STATES   YES\n")
    w.write("\n")
    w.write("[STAGE:MANAGE]\n")
    w.write("    STAGECLASS   phase.DBManageStage\n")
    w.write("    INPUTS   POSTLIGPREP_OUT,\n")
    w.write("    OUTPUTS   DATABASE,\n")
    w.write("    DATABASE  test.phdb\n")
    w.write("    NEW   YES\n")
    w.write("    MULTIPLE_CONFS   NO\n")
    w.write("    CONSIDER_STEREO   NO\n")
    w.write("    GENERATE_PROPS   NO\n")
    w.write("    CREATE_SUBSET   NO\n")
    w.write("    SKIP_DUPLICATES   NO\n")
    w.write("\n")
    w.write("[STAGE:CONFSITES]\n")
    w.write("    STAGECLASS   phase.DBConfSitesStage\n")
    w.write("    INPUTS   DATABASE,\n")
    w.write("    CONFS   auto\n")
    w.write("    MAX_CONFS   1\n")
    w.write("    GENERATE_PROPS   YES\n")
    w.write("\n")
    w.write("[USEROUTS]\n")
    w.write("    USEROUTS   DATABASE,\n")
    w.close()

def write_train_slurm(c):
    """
    Write slurm script for training

    :param control_dir: Work directory
    :param args: the args
    """
    jobname="verifytrain"
    w = open("submit_"+jobname+".sh","wt")
    w.write("#!/bin/bash\n")
    w.write("#SBATCH -J "+jobname+"\n")
    w.write("#SBATCH -o /dev/null\n")
    w.write("#SBATCH -e /dev/null\n")
    w.write("#SBATCH "+c.SLURM_GPU_PARAMETER+"\n")
    if c.SLURM_GPU_EXCLUSIVE == "1":
        w.write("#SBATCH --exclusive\n")
    w.write("cd $HOME/SPACEHASTEN/VERIFY\n")
    w.write(c.SLURM_PREPARE_ANACONDA + "\n")
    w.write(c.SLURM_ACTIVATE_CHEMPROP + "\n")
    w.write("chemprop_train --dataset_type regression --target_columns docking_score --data_path example.csv --save_dir verifytrain_model --batch_size 250 --no_cache_mol\n")
    w.write("touch jobdone-verifytrain\n")
    w.close()

def write_dock_slurm(c):
    jobname="verifydock_test"
    cpuname=jobname+"_cpu1"
    personal_scratch=c.SCRATCH_DEFAULT+"/"+os.getenv("USER")+"/"+cpuname
    w = open("submit_"+jobname+".sh","wt")
    w.write("#!/bin/bash\n")
    w.write("#SBATCH -J "+jobname+"\n")
    w.write("#SBATCH -o /dev/null\n")
    w.write("#SBATCH -e /dev/null\n")
    w.write("#SBATCH -p "+c.SLURM_PARTITION+"\n")
    w.write("#SBATCH -n 1\n")
    w.write("export SCHRODINGER_FEATURE_FLAGS=\"\"\n")
    w.write("curdir=$(pwd)\n")
    w.write("rm -fr "+personal_scratch+"\n")
    w.write("mkdir -p "+personal_scratch+"/"+cpuname+"\n")
    w.write("cp "+c.SPACEHASTEN_DIRECTORY+"/examples.smi "+personal_scratch+"/\n")
    w.write("cp "+c.SPACEHASTEN_DIRECTORY+"/grid-test_dock.zip "+personal_scratch+"/\n")
    w.write("cp test.inp "+personal_scratch+"/\n")
    w.write("cp test.in "+personal_scratch+"/\n")
    w.write("cd "+personal_scratch+"\n")
    w.write("$SCHRODINGER/pipeline -prog phase_db test.inp -OVERWRITE -WAIT -HOST localhost:1 -NJOBS 1\n")
    w.write("$SCHRODINGER/phase_database $(pwd)/test.phdb export -omae $(pwd)/test -get 1 -limit 99999999 -WAIT\n")
    w.write("rm -fr $(pwd)/test.phdb\n")
    w.write("$SCHRODINGER/glide test.in -OVERWRITE -WAIT -NJOBS 1 -HOST localhost:1\n")
    w.write("tar --exclude=results-test.tar.gz -czf results-test.tar.gz .\n")
    w.write("mv results-test.tar.gz $curdir\n")
    w.write("cd $curdir\n")
    w.write("rm -fr "+personal_scratch+"\n")
    w.write("touch jobdone-test-CPU1\n")
    w.close()

def check_slurm(c):
    curdir = os.getcwd()
    dock_dir = os.getenv("HOME")+"/SPACEHASTEN/VERIFY"
    os.chdir(dock_dir)
    print("Looking for the sbatch:")
    sbatch_ok = os.system("which sbatch")
    if sbatch_ok != 0:
        print("The sbatch command from slurm is not available.")
        exit()
    print("Running a simple job in "+c.SLURM_PARTITION+" slurm partition...")
    slurm_ok = os.system("sbatch --partition="+c.SLURM_PARTITION+" --wrap='echo Hello World!'")
    if slurm_ok != 0:
        print("Unable to run stuff in slurm partition "+c.SLURM_PARTITION+".")
        exit()
    print("Running a docking job in slurm partition "+c.SLURM_PARTITION+"...")
    write_confgen_file(c)
    w = open("test.in","wt")
    for line in open(c.SPACEHASTEN_DIRECTORY+"/test_dock.in"): w.write(line)
    w.write("LIGANDFILE   test_1.maegz\n")
    w.close()
    write_dock_slurm(c)

    os.system("rm -f jobdone-test-CPU1")
    print("Running docking via slurm at "+dock_dir+" ...")
    os.system("sbatch submit_verifydock_test.sh")
    os.chdir(curdir)
    jobs_left = 1 - len(glob.glob(dock_dir+"/jobdone-test-CPU*"))
    while jobs_left>0:
        time.sleep(5)
        jobs_left = 1 - len(glob.glob(dock_dir+"/jobdone-test-CPU*"))
    os.chdir(curdir)
    os.system("rm -fr "+c.SCRATCH_DEFAULT+"/"+os.getenv("USER")+"/verifydock_test_results_cpu1")
    os.system("mkdir -p "+c.SCRATCH_DEFAULT+"/"+os.getenv("USER")+"/verifydock_test_results_cpu1")
    os.system("tar xzf "+dock_dir+"/results-test.tar.gz -C "+c.SCRATCH_DEFAULT+"/"+os.getenv("USER")+"/verifydock_test_results_cpu1/")
    if not os.path.exists(c.SCRATCH_DEFAULT+"/"+os.getenv("USER")+"/verifydock_test_results_cpu1/test.csv"):
        print("Docking did not produce output file. Check that your computing nodes can see SpaceHASTEN installation nodes.)")
        print("Also, check that you do not have licensing issues.")
        exit()
    lines_in_output = len(open(c.SCRATCH_DEFAULT+"/"+os.getenv("USER")+"/verifydock_test_results_cpu1/test.csv").readlines())
    if lines_in_output < 2:
        print("Docking output is empty (test.csv).")
        exit()
    else:
        print("Docking output:",lines_in_output)
    print("Docking job finished successfully.")
    os.system("rm -fr "+c.SCRATCH_DEFAULT+"/"+os.getenv("USER")+"/verifydock_test_results_cpu1")
    
def check_training(c):
    curdir = os.getcwd()
    train_dir = os.getenv("HOME")+"/SPACEHASTEN/VERIFY"
    os.system("cp "+c.SPACEHASTEN_DIRECTORY+"/example.csv "+train_dir+"/")
    os.chdir(train_dir)
    os.system("rm -f jobdone-verifytrain")
    print("Testing chemprop training via slurm ...")
    write_train_slurm(c)
    os.system("sbatch submit_verifytrain.sh")
    os.chdir(curdir)
    jobs_left = 1 - len(glob.glob(train_dir+"/jobdone-verifytrain*"))
    while jobs_left>0:
        time.sleep(5)
        jobs_left = 1 - len(glob.glob(train_dir+"/jobdone-verifytrain*"))
    if not os.path.exists(train_dir+"/verifytrain_model/test_scores.csv"):
        print("Model not created.")
        exit()
    model_specs = open(train_dir+"/verifytrain_model/test_scores.csv").readlines()[-1]
    print("Model specs:",model_specs)
    cuda = True
    for line in open(train_dir+"/verifytrain_model/verbose.log"):
        if "'cuda': False," in line:
            cuda = False
    if not cuda:
        print("**********************************************************")
        print("WARNING!!!!!! CUDA not available, training will be slower.")
        print("**********************************************************")
    else:
        print("Training was running on CUDA.")
    print("Model trained succesfully.")

def check_biosolveit(c):
    print("Using the default chemical space:",c.SPACES_FILE_DEFAULT)
    if not os.path.exists(c.SPACES_FILE_DEFAULT):
        print("This file is missing.")
        exit()
    verdir = c.SCRATCH_DEFAULT+"/"+os.getenv("USER")+"/verify_biosolveit"
    os.system("rm -fr "+verdir)
    os.mkdir(verdir)
    
    print("Running locally 1 Spacelight....")
    spacelight_ok = os.system(c.EXE_SPACELIGHT_DEFAULT + " -i " + c.SPACEHASTEN_DIRECTORY + "/example.smi -s "+c.SPACES_FILE_DEFAULT +" -o "+verdir+"/spacelightresult.csv --max-nof-results 100 --min-similarity-threshold 0.5 --thread-count 1")
    if spacelight_ok != 0 or not os.path.exists(verdir+"/spacelightresult_1.csv"):
        print("Spacelight failed.")
        exit()
    print("SpaceLight output:",len(open(verdir+"/spacelightresult_1.csv").readlines()))
    print("Running locally 1 FTrees...")
    ftrees_ok = os.system(c.EXE_FTREES_DEFAULT +  " -i " + c.SPACEHASTEN_DIRECTORY + "/example.smi -s "+c.SPACES_FILE_DEFAULT+" -o "+verdir+"/ftreesresult.csv --max-nof-results 100 --min-similarity-threshold 0.9 --thread-count 1")
    if ftrees_ok != 0 or not os.path.exists(verdir+"/ftreesresult_1.csv"):
        print("Ftrees failed.")
        exit()
    print("FTrees output:",len(open(verdir+"/ftreesresult_1.csv").readlines()))
    os.system("rm -fr "+verdir)
    print("SpaceLight and FTrees tested succesfully.")

def check_pigz():
    print("Looking for pigz:")
    pigz_ok = os.system("which pigz")
    if pigz_ok != 0:
        print("pigz not installed.")
        exit()

print("This script checks that all bits and pieces required to run SpaceHASTEN are in place.")

if os.path.exists(os.getenv("HOME")+"/SPACEHASTEN/VERIFY"):
    print("ERROR: please remove previous verification directory $HOME/SPACEHASTEN/VERIFY before re-run!")
    exit()

c = cfg.SpaceHASTENConfiguration()
print("SpaceHASTEN directory:",c.SPACEHASTEN_DIRECTORY)
print("Creating $HOME/SPACEHASTEN/VERIFY directory that should be visible to all computing nodes as well...")
os.system("mkdir -p $HOME/SPACEHASTEN/VERIFY")

check_pigz()
check_slurm(c)
check_training(c)
check_biosolveit(c)

print("All checks passed. You are ready to run SpaceHASTEN ("+c.SPACEHASTEN_DIRECTORY+"/spacehasten)")
