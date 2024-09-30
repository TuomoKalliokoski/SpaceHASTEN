# SpaceHASTEN: configuration
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
import sys
import configparser
import shutil
import os

class SpaceHASTENConfiguration:
      
    SPACEHASTEN_DIRECTORY = "/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])
    CONTROL_EXE = SPACEHASTEN_DIRECTORY+"/control.py"
    CHUNKPREDICT_EXE = SPACEHASTEN_DIRECTORY+"/chunkpredict.py"
    EXPORTPOSES_EXE = "$SCHRODINGER/run " + SPACEHASTEN_DIRECTORY + "/export_poses.py"
    SPACEHASTEN_VERSION=0.1
    MAX_CORES = 250
    EXE_SPACELIGHT_DEFAULT = "/data/programs/BiosolveIT/spacelight-1.3.0-Linux-x64/spacelight"
    EXE_FTREES_DEFAULT = "/data/programs/BiosolveIT/ftrees-6.11.0-Linux-x64/ftrees"
    SPACES_DIR_DEFAULT = "/data/programs/BiosolveIT/spaces"
    SPACES_FILE_DEFAULT = "/data/programs/BiosolveIT/spaces/REALSpace_48bn_2024-02.space"
    QUERIES_DEFAULT = 1000
    DOCKING_DEFAULT = 1000000
    DOCKING_CHUNK = 1000
    RDKIT_CHUNK_DEFAULT = 12345
    RDKIT_CPU_DEFAULT = 0
    PROP_MW_MIN_DEFAULT = 0.0
    PROP_MW_MAX_DEFAULT = 500.0
    PROP_SLOGP_MIN_DEFAULT = -10.0
    PROP_SLOGP_MAX_DEFAULT = 5.0
    PROP_HBA_MIN_DEFAULT = 0
    PROP_HBA_MAX_DEFAULT = 10
    PROP_HBD_MIN_DEFAULT = 0
    PROP_HBD_MAX_DEFAULT = 5
    PROP_ROTBONDS_MIN_DEFAULT = 0
    PROP_ROTBONDS_MAX_DEFAULT = 10
    PROP_TPSA_MIN_DEFAULT = 0.0
    PROP_TPSA_MAX_DEFAULT = 140.0
    CHEMPROP_CHUNK_DEFAULT = 12345
    CHEMPROP_CPU_DEFAULT = 250
    SCRATCH_DEFAULT = "/wrk"
    SIM_SPACELIGHT_DEFAULT = 0.5
    SIM_FTREES_DEFAULT = 0.9
    NNN_DEFAULT = 10000
    TRAIN_DOCKING_CUTOFF = 10.0
    FIELD_SMILES_DEFAULT = "SMILES"
    FIELD_TITLE_DEFAULT = "title"
    FIELD_SCORE_DEFAULT = "r_i_docking_score"
    SLURM_PARTITION = "jobs"
    SLURM_PREPARE_ANACONDA = "source /data/programs/oce/actoce"
    SLURM_ACTIVATE_CHEMPROP = "conda activate chemprop"
    SLURM_GPU_PARAMETER = "-w lfies-docki"
    SLURM_GPU_EXCLUSIVE = "1"
    SLURM_CPU_COUNT_SEARCH = "2"
    SLURM_CPU_COUNT_DOCK = "1"
    SLURM_CPU_COUNT_PREDICT = "1"
    SLURM_CPU_COUNT_CONTROL = "1"
    
    ENAMINEREAL_SEEDS = "/data/tuomo/PROJECTS/SPACEHASTEN/Enamine_Diverse_REAL_drug-like_48.2M_cxsmiles.cxsmiles.bz2"
    ENAMINEREAL_SEEDS_COUNT = 1000000
    ENAMINEREAL_SEEDS_CPU = 4
    
    def __init__(self):
        cparser = configparser.ConfigParser()
        cparser.read(self.SPACEHASTEN_DIRECTORY+"/spacehasten.ini")

        for setting in cparser["Paths"]:
            if setting == "exe_spacelight_default":
                self.EXE_SPACELIGHT_DEFAULT = cparser["Paths"][setting]
            elif setting == "exe_ftrees_default":
                self.EXE_FTREES_DEFAULT = cparser["Paths"][setting]
            elif setting == "scratch_default":
                self.SCRATCH_DEFAULT = cparser["Paths"][setting]
            elif setting == "spaces_dir_default":
                self.SPACES_DIR_DEFAULT = cparser["Paths"][setting]
            elif setting ==  "spaces_file_default":
                self.SPACES_FILE_DEFAULT = cparser["Paths"][setting]
            elif setting ==  "enaminereal_seeds":
                self.ENAMINEREAL_SEEDS = cparser["Paths"][setting]
            else:
                print("Error: Unknown setting in spacehasten.ini:",setting)
                sys.exit(1)

        for setting in cparser["Slurm"]:
            if setting == "slurm_partition":
                self.SLURM_PARTITION = cparser["Slurm"][setting]
            elif setting == "slurm_prepare_anaconda":
                self.SLURM_PREPARE_ANACONDA = cparser["Slurm"][setting]
            elif setting == "slurm_activate_chemprop":
                self.SLURM_ACTIVATE_CHEMPROP = cparser["Slurm"][setting]
            elif setting == "slurm_gpu_parameter":
                self.SLURM_GPU_PARAMETER = cparser["Slurm"][setting]
            elif setting == "slurm_gpu_exclusive":
                self.SLURM_GPU_EXCLUSIVE = cparser["Slurm"][setting]
            elif setting == "slurm_cpu_count_search":
                self.SLURM_CPU_COUNT_SEARCH = cparser["Slurm"][setting]
            elif setting == "slurm_cpu_count_dock":
                self.SLURM_CPU_COUNT_DOCK = cparser["Slurm"][setting]
            elif setting == "slurm_cpu_count_predict":
                self.SLURM_CPU_COUNT_PREDICT = cparser["Slurm"][setting]
            elif setting == "slurm_cpu_count_control":
                self.SLURM_CPU_COUNT_CONTROL = cparser["Slurm"][setting]
            else:
                print("Error: Unknown setting in spacehasten.ini:",setting)
                sys.exit(1)

        for setting in cparser["Properties"]:
            if setting ==  "mw_min":
                self.PROP_MW_MIN_DEFAULT = float(cparser["Properties"][setting])
            elif setting == "mw_max":
                self.PROP_MW_MAX_DEFAULT = float(cparser["Properties"][setting])
            elif setting ==  "slogp_min":
                self.PROP_SLOGP_MIN_DEFAULT = float(cparser["Properties"][setting])
            elif setting == "slogp_max":
                self.PROP_SLOGP_MAX_DEFAULT = float(cparser["Properties"][setting])
            elif setting == "hba_min":
                self.PROP_HBA_MIN_DEFAULT = int(cparser["Properties"][setting])
            elif setting == "hba_max":
                self.PROP_HBA_MAX_DEFAULT = int(cparser["Properties"][setting])
            elif setting == "hbd_min":
                self.PROP_HBD_MIN_DEFAULT = int(cparser["Properties"][setting])
            elif setting == "hbd_max":
                self.PROP_HBD_MAX_DEFAULT = int(cparser["Properties"][setting])
            elif setting == "rotbonds_min":
                self.PROP_ROTBONDS_MIN_DEFAULT = int(cparser["Properties"][setting])
            elif setting == "rotbonds_max":
                self.PROP_ROTBONDS_MAX_DEFAULT = int(cparser["Properties"][setting])
            elif setting == "tpsa_min":
                self.PROP_TPSA_MIN_DEFAULT = float(cparser["Properties"][setting])
            elif setting == "tpsa_max":
                self.PROP_TPSA_MAX_DEFAULT = float(cparser["Properties"][setting])
            else:
                print("Error: Unknown setting in spacehasten.ini:",setting)
                sys.exit(1)

        if shutil.which("chemprop_predict") is None:
            sys.exit("Error: chemprop_predict not available, please active chemprop environment before SpaceHASTEN!")
        if shutil.which("chemprop_train") is None:
            sys.exit("Error: chemprop_train not available, please active chemprop environment before SpaceHASTEN!")
        if shutil.which("bzcat") is None:
            sys.exit("Error: bzcat not available, please install bzip2 package!")
        if not os.path.exists(self.EXE_SPACELIGHT_DEFAULT):
            sys.exit("Error: "+self.EXE_SPACELIGHT_DEFAULT+" not available!")                      
        if not os.path.exists(self.EXE_FTREES_DEFAULT):
            sys.exit("Error: "+self.EXE_FTREES_DEFAULT+" not available!")                      
        if not os.path.exists(os.getenv("SCHRODINGER")+"/run"):
            sys.exit("Error: $SCHRODINGER/run not available.")
        if not os.path.exists(self.SPACEHASTEN_DIRECTORY+"/control.py"):
            sys.exit("Error: invalid SpaceHASTEN installation, control.py missing!")
        if not os.path.exists(self.SPACEHASTEN_DIRECTORY+"/chunkpredict.py"):
            sys.exit("Error: invalid SpaceHASTEN installation, chunkpredict.py missing!")
        if not os.path.exists(self.SPACEHASTEN_DIRECTORY+"/export_poses.py"):
            sys.exit("Error: invalid SpaceHASTEN installation, export_poses.py missing!")
        if not os.path.exists(self.SPACEHASTEN_DIRECTORY+"/spacehasten_logo.png"):
            sys.exit("Error: invalid SpaceHASTEN installation, spacehasten_logo missing!")
        
        
