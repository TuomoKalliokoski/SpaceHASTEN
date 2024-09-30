# SpaceHASTEN: docking functions
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
import sqlite3
import functions
import tqdm
import glob
import pandas as pd
import time
import slurm_functions
import multiprocessing as mp
import random

def process_docking_results(args,dock_iteration):
    dbname = args.name+".dbsh"
    dock_dir = os.getenv("HOME")+"/SPACEHASTEN/DOCKING_"+args.name+"_iter"+str(dock_iteration)
    print("\nDecompressing results to a local drive...")
    resdir=args.scratch+"/"+os.getenv("USER")+"/COLLECTdock_"+args.name+"_iter"+str(dock_iteration)
    os.system("rm -fr "+resdir)
    os.system("mkdir -p "+resdir)
    tars = []
    for result_tar in glob.glob(dock_dir+"/results-*.tar.gz"):
        tars.append("tar -xzf "+result_tar+" -C "+resdir)
    pool = mp.Pool(mp.cpu_count())
    list(tqdm.tqdm(pool.imap_unordered(os.system,tars),mininterval=1,total=len(tars)))
    print("Processing docking scores...")
    if len(glob.glob(resdir+"/glide_*.csv")) == 0:
        print("ERROR: No docking results found! Please check that your docking input file is valid.")
        os.system("rm -fr "+resdir)
        return
    docking_results = {}
    for resfile in glob.glob(resdir+"/glide_*.csv"):
        if "_skip.csv" in resfile:
            continue
        resdata = pd.read_csv(resfile)
        best_scores = resdata[["title","r_i_docking_score"]].groupby(["title"]).min().to_dict()["r_i_docking_score"]
        docking_results.update(best_scores)
    os.system("rm -fr "+resdir)        
    print("Storing docking scores to SQLite3...")
    to_db = []
    for spacehastenid in docking_results: to_db.append((docking_results[spacehastenid],dock_iteration,spacehastenid))
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    c.executemany("UPDATE data SET dock_score = ?,dock_iteration = ? WHERE spacehastenid = ?",to_db)
    conn.commit()
    conn.close()
    
def dock(args,importing_seeds=False,do_not_update_gui=False):
    """
    Do the confgen + docking

    :param args: the args
    :param importing_seeds: importing seeds or not
    :param do_not_update_gui: do not update GUI
    """
    os.system("date")
    print("Docking:")
    dbname = args.name+".dbsh"
    
    if args.top is None:
        raise SystemExit("Error: --top not defined.")
    if args.cpu is None:
        raise SystemExit("Error: number of cores not defined with --cpu.")
    if not os.path.exists(dbname):
        raise SystemExit("Error: SpaceHASTEN database ("+dbname+") missing, have you imported the seeds?")
    if args.sff>0:
        top_spacelight=round(float(args.top)*args.sff)
        top_ftrees=args.top-top_spacelight
        print("Spacelight picked compounds:",top_spacelight)
        print("Ftrees picked compounds:",top_ftrees)

    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    to_be_docked = []
    docked_ids = set()
    """
    if args.sff>0:
        # Obsolete code, but I keep it for reference
        for smiles,spacehastenid,pred_score in c.execute("SELECT smiles,spacehastenid,pred_score FROM data WHERE dock_score IS NULL AND spacelight IS NOT NULL ORDER BY pred_score LIMIT "+str(top_spacelight)).fetchall():
            to_be_docked.append(smiles.strip() + " " + str(spacehastenid) + "\n")
            docked_ids.add(spacehastenid)
        compounds_left = top_ftrees
        for smiles,spacehastenid,pred_score in c.execute("SELECT smiles,spacehastenid,pred_score FROM data WHERE dock_score IS NULL and ftrees IS NOT NULL ORDER BY pred_score"):
            if spacehastenid in docked_ids:
                continue
            compounds_left -= 1
            to_be_docked.append(smiles.strip() + " " + str(spacehastenid) + "\n")
            if compounds_left <= 0:
                break
    """
    
    if importing_seeds:
        print("Picking all imported compounds for docking...")
        sql_query = "SELECT smiles,spacehastenid FROM data WHERE dock_score IS NULL"
    else:
        print("Picking top predicted compounds for docking...")
        sql_query = "SELECT smiles,spacehastenid FROM data WHERE pred_score IS NOT NULL AND dock_score IS NULL ORDER BY pred_score LIMIT "+str(args.top)
    for smiles,spacehastenid in c.execute(sql_query).fetchall():
        to_be_docked.append(smiles.strip() + " " + str(spacehastenid) + "\n")
    conn.close()
    print("Shuffling compounds before docking...")
    random.shuffle(to_be_docked)

    #smiles_per_cpu = round(float(len(to_be_docked))/float(args.cpu))
    
    print("\nCompounds to be docked:")
    print("Total number:",len(to_be_docked))
    #print("SMILES per CPU:",smiles_per_cpu)
    print("SMILES per docking job:",args.c.DOCKING_CHUNK)

    if importing_seeds:
        dock_iteration = 0
    else:
        dock_iteration = functions.get_latest_iteration(args.name)+1

    dock_dir = os.getenv("HOME")+"/SPACEHASTEN/DOCKING_"+args.name+"_iter"+str(dock_iteration)
    print("Writing work files to",dock_dir)
    os.system("mkdir -p "+dock_dir)
    chunk_counter = 0
    smiles_per_chunk = 0
    new_file = True
    for line in tqdm.tqdm(to_be_docked,mininterval=1):
        if new_file:
            chunk_counter += 1
            smiles_per_chunk = 0
            w = open(dock_dir+"/dockinput_"+args.name+"_iter"+str(dock_iteration)+"_cpu"+str(chunk_counter)+".smi","wt")
            write_confgen_file(dock_dir+"/dockinput_"+args.name+"_iter"+str(dock_iteration)+"_cpu"+str(chunk_counter)+".inp")
            write_docking_file(dock_dir+"/dockinput_"+args.name+"_iter"+str(dock_iteration)+"_cpu"+str(chunk_counter)+".in",dbname,dock_dir)
            new_file = False
        w.write(line)
        smiles_per_chunk += 1
        if smiles_per_chunk >= args.c.DOCKING_CHUNK:
            w.close()
            new_file = True
    if not new_file:
        w.close()
    slurm_functions.write_dock_slurm(dock_dir,args,chunk_counter)
    curdir = os.getcwd()
    os.chdir(dock_dir)
    os.system("rm -f jobdone-"+args.name+"-CPU*")
    print("Running docking via slurm at "+dock_dir+" ...")
    os.system("sbatch submit_dockinput_"+args.name+"_iter"+str(dock_iteration)+".sh")
    os.chdir(curdir)
    #jobs_left = args.cpu - len(glob.glob(dock_dir+"/jobdone-"+args.name+"-CPU*"))
    jobs_left = chunk_counter - len(glob.glob(dock_dir+"/jobdone-"+args.name+"-CPU*"))
    while jobs_left>0:
        time.sleep(5)
        #jobs_left = args.cpu - len(glob.glob(dock_dir+"/jobdone-"+args.name+"-CPU*"))
        jobs_left = chunk_counter - len(glob.glob(dock_dir+"/jobdone-"+args.name+"-CPU*"))
    process_docking_results(args,dock_iteration)
    print("Docking complete!")
    os.system("date")
    if not do_not_update_gui:
        args.q.put("Percent:99.9")
        args.q.put("UpdateIteration:"+str(dock_iteration))
        args.q.put("DoneDocking")
    
def write_confgen_file(filename):
    """
    Write .inp for Phase/LigPrep

    :param filename: Filename
    """
    basename=filename.split("/")[-1]
    w = open(filename,"wt")
    w.write("[SET:ORIGINAL_LIGANDS]\n")
    w.write("    VARCLASS   Structures\n")
    w.write("    FILES   "+basename.replace(".inp",".smi")+",\n")
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
    w.write("    DATABASE  "+basename.replace(".inp",".phdb")+"\n")
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

def write_docking_file(filename,dbname,dock_dir):
    """
    Write .in and .zip for Glide

    :param filename: Filename
    :param dbase: Template is located in database
    :param dock_dir: File output directory
    """
    
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    template_blob = c.execute("SELECT dock_param FROM docking_param").fetchall()[0][0]
    grid_blob = c.execute("SELECT dock_grid FROM docking_grid").fetchall()[0][0]
    conn.close()
    with open("template.tmp","wb") as f:
        f.write(template_blob)
    with open(dock_dir+"/glide_grid.zip","wb") as f:
        f.write(grid_blob)
    template_lines = []
    for line in open("template.tmp","rt"):
        l = line.strip().split()
        if len(l)>0 and l[0] != "LIGANDFILE" and l[0] != "GRIDFILE":
            template_lines.append(line)
    os.system("rm template.tmp")
    basename=filename.split("/")[-1]
    w = open(dock_dir+"/glide_"+basename,"wt")
    w.write("LIGANDFILE   "+basename.replace(".in","_1.maegz\n"))
    w.write("GRIDFILE     glide_grid.zip\n")
    for template_line in template_lines:
        w.write(template_line)
    w.close()

