# SpaceHASTEN: functions to do the similarity searching
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
import functions
import sqlite3
import slurm_functions
import prediction_functions
import glob
import time
import pandas as pd
import tqdm
import gzip
import multiprocessing as mp

def remove_existing(filtered_mols,args):
    """
    Remove existing compounds from Spacelight/FTrees results

    :param filtered_mols: Searching results from Spacelight/FTrees
    :param args: the args
    :return: List of molecules
    """
    dbname = args.name + ".dbsh"
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    to_db = []
    for line in filtered_mols:
        l = line.split("§")
        reghash = l[0]
        smiles = l[1]
        title = l[2]
        existing_spacehastenid = c.execute("SELECT spacehastenid FROM data WHERE reghash = '"+reghash+"'").fetchall()
        if len(existing_spacehastenid)==0:
            to_db.append((reghash,smiles,title))
    return to_db

def process_sim_results(cycle_dir,args):
    """
    Process results from spacelight and ftrees

    :param cycle_dir: cycle dir
    :param args: the args
    """
    sim_methods = ["spacelight","ftrees"]
    raw_mols = {}
    sims = {}
    for sim_method in sim_methods:
        sims[sim_method] = {}
        print("Reading in "+sim_method+" results...")
        duplicates = 0
        # group by title here
        for resfile in tqdm.tqdm(glob.glob(cycle_dir+"/"+sim_method+"result_"+args.name+"_*_1.csv"),mininterval=1):
            resdata = pd.read_csv(resfile)
            for smiles,title,similarity in resdata[["#result-smiles","result-name","similarity"]].values.tolist():
                if title not in sims[sim_method] or sims[sim_method][title] <= similarity:
                    sims[sim_method][title] = similarity
                if title not in raw_mols:
                    raw_mols[title] = smiles + "§" + title
                else:
                    duplicates += 1
        print(len(raw_mols),"before property control added from",sim_method)
    # split data into smaller chunks for slurm
    smiles_per_cpu = round(float(len(raw_mols))/float(args.cpu))
    print("Splitting molecules, SMILES per CPU:",smiles_per_cpu)

    cpu_counter = 0
    
    os.system("mkdir -p "+cycle_dir+"/CONTROL")
    os.system("rm -f "+cycle_dir+"/CONTROL/control_*.smi.gz")
    model_version = functions.get_latest_model(args.name)
    conn = sqlite3.connect(args.name + ".dbsh")
    c = conn.cursor()
    model_blob = c.execute("SELECT model_tar FROM models WHERE model_version = "+str(model_version)).fetchall()[0][0]
    model_name = "model_"+args.name+"_ver"+str(model_version)
    with open(model_name+".tar.gz","wb") as f:
        f.write(model_blob)
    os.system("tar -xzf "+model_name+".tar.gz -C "+cycle_dir+"/CONTROL/")
    os.remove(model_name+".tar.gz")
    conn.close()

    new_file = True
    chunk_size = 0
    for title in tqdm.tqdm(raw_mols,mininterval=1):
        if new_file:
            cpu_counter += 1
            csv_filename = cycle_dir + "/CONTROL/control_"+args.name+"_cpu"+str(cpu_counter)+".smi.gz"
            w = gzip.open(csv_filename,"wt")
            new_file = False
        w.write(raw_mols[title]+"\n")
        chunk_size += 1
        if chunk_size >= smiles_per_cpu:
            chunk_size = 0
            w.close()
            new_file = True
    if not new_file:
        w.close()
    
    slurm_functions.write_control_slurm(cycle_dir+"/CONTROL",args)
    w = open(cycle_dir+"/CONTROL/control.param","wt")
    w.write(str(args.c.PROP_MW_MIN_DEFAULT)+"\n")
    w.write(str(args.c.PROP_MW_MAX_DEFAULT)+"\n")
    w.write(str(args.c.PROP_SLOGP_MIN_DEFAULT)+"\n")
    w.write(str(args.c.PROP_SLOGP_MAX_DEFAULT)+"\n")
    w.write(str(args.c.PROP_HBA_MIN_DEFAULT)+"\n")
    w.write(str(args.c.PROP_HBA_MAX_DEFAULT)+"\n")
    w.write(str(args.c.PROP_HBD_MIN_DEFAULT)+"\n")
    w.write(str(args.c.PROP_HBD_MAX_DEFAULT)+"\n")
    w.write(str(args.c.PROP_ROTBONDS_MIN_DEFAULT)+"\n")
    w.write(str(args.c.PROP_ROTBONDS_MAX_DEFAULT)+"\n")
    w.write(str(args.c.PROP_TPSA_MIN_DEFAULT)+"\n")
    w.write(str(args.c.PROP_TPSA_MAX_DEFAULT)+"\n")
    w.close()

    curdir = os.getcwd()
    os.chdir(cycle_dir+"/CONTROL")
    os.system("rm -f jobdone-"+args.name+"-CPU*")
    print("Controlling properties and predicting docking_score via slurm at "+cycle_dir+"/CONTROL ...")
    os.system("sbatch submit_ctrl_"+args.name+"_cycle"+str(functions.get_latest_cycle(args.name))+".sh")
    os.chdir(curdir)

    jobs_left = args.cpu - len(glob.glob(cycle_dir+"/CONTROL/jobdone-"+args.name+"-CPU*"))
    while jobs_left>0:
        time.sleep(5)
        jobs_left = args.cpu - len(glob.glob(cycle_dir+"/CONTROL/jobdone-"+args.name+"-CPU*"))

    print("Reading in predictions...")
    prop_inputs = []
    for prop_input in glob.glob(cycle_dir+"/CONTROL/predicted_propoutput_control_*.csv"):
        prop_inputs.append(prop_input)
    pool = mp.Pool(mp.cpu_count())
    filtered_mol_runs_with_predictions = list(filter(None,list(tqdm.tqdm(pool.imap_unordered(functions.get_rdkit_properties,prop_inputs),mininterval=1,total=len(prop_inputs)))))
    docking_scores = {}
    filtered_mols = []
    for cpu in range(len(filtered_mol_runs_with_predictions)):
        rawmol_list,docking_score_list = filtered_mol_runs_with_predictions[cpu]
        for comp_index in range(len(rawmol_list)):
            title = rawmol_list[comp_index].split("§")[2]
            docking_scores[title] = float(docking_score_list[comp_index])
            filtered_mols.append(rawmol_list[comp_index])

    filtered_sims = {}
    for sim_method in sim_methods:
        filtered_sims[sim_method] = {}
    for filtered_mol in filtered_mols:
        title = filtered_mol.split("§")[2]
        for sim_method in sim_methods:
            if title in sims[sim_method]:
                filtered_sims[sim_method][title] = sims[sim_method][title]

    return (filtered_mols,filtered_sims,docking_scores)

def simsearch(args,do_not_update_gui=False):
    """
    Do the similarity searching

    :param args: the args
    :param do_not_update_gui: do not update GUI
    """
    print("Similarity searching:")
    os.system("date")
    dbname = args.name + ".dbsh"
    
    if not os.path.exists(os.getenv("HOME")+"/SPACEHASTEN"):
        os.mkdir(os.getenv("HOME")+"/SPACEHASTEN")

    conn = sqlite3.connect(dbname)
    c = conn.cursor()

    cycle_number = functions.get_latest_cycle(args.name)+1
    cycle_dir = os.getenv("HOME")+"/SPACEHASTEN/SIMSEARCH_"+args.name+"_cycle"+str(cycle_number)
    os.system("mkdir -p "+cycle_dir)
    print("Starting searching cycle",cycle_number)
    print("Picking",args.top,"top docked/predicted compounds for similarity searching queries...")
    if not args.use_predicted:
        print("Using docking_score")
        to_search = c.execute("SELECT smiles,spacehastenid FROM data WHERE query IS NULL AND dock_score IS NOT NULL ORDER BY dock_score LIMIT ?",[args.top]).fetchall()
    else:
        print("Using predicted score")
        to_search = c.execute("SELECT smiles,spacehastenid FROM data WHERE query IS NULL AND pred_score IS NOT NULL AND dock_score IS NULL ORDER by pred_score LIMIT ?",[args.top]).fetchall()
    for smiles,spacehastenid in to_search: c.execute("UPDATE data SET query = "+str(cycle_number)+" WHERE spacehastenid = "+str(spacehastenid))
    conn.commit()
    conn.close()
    w = open(cycle_dir+"/queries_"+args.name+".smi","wt")
    for smiles,spacehastenid in to_search: w.write(smiles.strip() + " " + str(spacehastenid) + "\n")
    w.close()
    slurm_functions.write_search_slurm(cycle_dir,args)

    curdir = os.getcwd()
    os.chdir(cycle_dir)
    os.system("rm -f jobdone-"+args.name+"-CPU*")
    print("Running similarity searching via slurm...")
    os.system("sbatch submit_queries_"+args.name+".sh")
    os.chdir(curdir)
    jobs_left = args.top - len(glob.glob(cycle_dir+"/jobdone-"+args.name+"-CPU*"))
    while jobs_left>0:
        time.sleep(5)
        jobs_left = args.top - len(glob.glob(cycle_dir+"/jobdone-"+args.name+"-CPU*"))

    filtered_mols,filtered_sims,predicted_scores = process_sim_results(cycle_dir,args)
    print(len(filtered_mols),"property-matched mols")
    sim_methods = ["spacelight","ftrees"]
    for sim_method in sim_methods:
        print(sim_method,":",len(filtered_sims[sim_method]))
    only_in_spacelight = 0
    only_in_ftrees = 0
    for smilesid in filtered_sims["spacelight"]:
        if smilesid not in filtered_sims["ftrees"]:
            only_in_spacelight += 1
    for smilesid in filtered_sims["ftrees"]:
        if smilesid not in filtered_sims["spacelight"]:
            only_in_ftrees += 1
    print("Unique by IDs, SpaceLight",only_in_spacelight,", FTrees:",only_in_ftrees)

    print("Making sure that new compounds have unique RegistrationHashes...")
    df_filtered_mols = pd.DataFrame(filtered_mols,columns=["rawmol"])
    df_filtered_mols[["RegHash","SMILES","SMILES_ID"]] = df_filtered_mols["rawmol"].str.split("§",expand=True)
    reghash_unique_filtered_mols = list(df_filtered_mols.groupby("RegHash").first()["rawmol"])
    print(df_filtered_mols.shape[0],"before unique reghash check")
    print(len(reghash_unique_filtered_mols),"after unique reghash check")

    print("Checking which are already discovered...")
    new_filtered_mols = remove_existing(reghash_unique_filtered_mols,args)
    
    print("Adding new compounds to SQLite3...")
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    to_db = []
    for reghash,smiles,smilesid in new_filtered_mols:
        spacelight_similarity = None
        ftrees_similarity = None
        if smilesid in filtered_sims["spacelight"]:
            spacelight_similarity = filtered_sims["spacelight"][smilesid]
        if smilesid in filtered_sims["ftrees"]:
            ftrees_similarity = filtered_sims["ftrees"][smilesid]
        to_db.append((reghash,smiles.strip(),smilesid,spacelight_similarity,ftrees_similarity,predicted_scores[smilesid],cycle_number))
    c.executemany("INSERT INTO data(reghash,smiles,smilesid,spacelight,ftrees,pred_score,simsearch_cycle) VALUES (?,?,?,?,?,?,?)",to_db)
    conn.commit()
    conn.close()
    print("\nSimilarity seaching done!")
    os.system("date")
    if not do_not_update_gui:
        args.q.put("DoneSimsearch")
