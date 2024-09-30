# SpaceHASTEN: functions to do the predicting
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
import pandas as pd
import tqdm
import glob
import time
import multiprocessing as mp
import slurm_functions

def get_chemprop_predictions(param):
    """
    Read chemprop output

    :param param: the csv inputfile, args
    :return identifiers added with pred_score
    """
    csv_filename = param[0]
    args = param[1]
    outputfile = "output_"+csv_filename
    os.system("OMP_NUM_THREADS=1 chemprop_predict --no_cache_mol --no_cuda --num_workers 0 --test_path "+csv_filename+" --checkpoint_dir model_"+args.name+"_ver"+str(functions.get_latest_model(args.name))+" --preds_path "+outputfile +" >/dev/null 2>/dev/null")
    predictions = pd.read_csv(outputfile)
    predictions.drop(columns=["smiles"],inplace=True)
    os.remove(outputfile)
    return dict(zip(predictions["smilesid"],predictions["docking_score"]))

def predict_dock(mols,args):
    """
    Predict docking scores using chemprop

    :param mols: molecules to predict
    :param args: the args
    :return mols added with pred_score
    """
    chunk_size = 0
    cpu_counter = 0

    new_file = True
    pred_scores = {}

    os.system("rm -f pred_"+args.name+"_cycle"+str(functions.get_latest_cycle(args.name))+"_*.csv")
    for reghash,smiles,title in mols:
        if new_file:
            cpu_counter += 1
            csv_filename = "pred_"+args.name+"_cycle"+str(functions.get_latest_cycle(args.name))+"_"+str(cpu_counter)+".csv"
            w = open(csv_filename,"wt")
            w.write("smiles,smilesid\n")
            new_file = False
        w.write(smiles.strip()+","+title+"\n")
        chunk_size += 1
        if chunk_size >= args.chemprop_chunk:
            chunk_size = 0 
            w.close()
            new_file = True
    if not new_file:
        w.close()

    pred_inputs = []
    for pred_input in glob.glob("pred_"+args.name+"_cycle"+str(functions.get_latest_cycle(args.name))+"_*.csv"):
        pred_inputs.append((pred_input,args))
    if args.rdkit_cpu > 0:
        cores = args.rdkit_cpu
    else:
        cores = mp.cpu_count()
    pool = mp.Pool(cores)
    predictions = list(filter(None,list(tqdm.tqdm(pool.imap_unordered(get_chemprop_predictions,pred_inputs),mininterval=1,total=len(pred_inputs)))))

    pred_scores = {}
    for prediction in predictions:
        pred_scores.update(prediction)

    os.system("rm -f pred_"+args.name+"_cycle"+str(functions.get_latest_cycle(args.name))+"_*.csv")
    return pred_scores

def update_predicted_scores(args):
    """
    Update predicted scores for all non-docked compounds in the database
    
    :param args: the args
    """
    dbname = args.name + ".dbsh"

    print("Predicting pred_score for all non-docked compounds with model version "+str(functions.get_latest_model(args.name)))
    conn = sqlite3.connect(dbname)
    c = conn.cursor()

    cycle_dir = os.getenv("HOME")+"/SPACEHASTEN/SIMSEARCH_"+args.name+"_cycle"+str(functions.get_latest_cycle(args.name))

    to_update = c.execute("SELECT smiles,spacehastenid FROM data WHERE dock_score IS NULL").fetchall()
    if len(to_update) == 0:
        print("Note: No undocked compounds to update.")
        return
    # split data into smaller chunks for slurm
    smiles_per_cpu = round(float(len(to_update))/float(args.cpu))

    os.system("mkdir -p "+cycle_dir+"/PREDICT")
    print("Loading model...")
    model_version = functions.get_latest_model(args.name)
    model_blob = c.execute("SELECT model_tar FROM models WHERE model_version = "+str(model_version)).fetchall()[0][0]
    model_name = "model_"+args.name+"_ver"+str(model_version)
    with open(model_name+".tar.gz","wb") as f:
        f.write(model_blob)
    os.system("tar -xzf "+model_name+".tar.gz -C "+cycle_dir+"/PREDICT/")
    os.remove(model_name+".tar.gz")

    print("Splitting molecules, SMILES per CPU:",smiles_per_cpu)
    os.system("rm -f "+cycle_dir+"/PREDICT/predict_*.csv")

    cpu_counter = 0
    new_file = True
    chunk_size = 0

    for smiles,spacehastenid in tqdm.tqdm(to_update,mininterval=1):
        if new_file:
            cpu_counter += 1
            csv_filename = cycle_dir + "/PREDICT/predict_"+args.name+"_cpu"+str(cpu_counter)+".csv"
            w = open(csv_filename,"wt")
            w.write("smiles,smilesid\n")
            new_file = False
        w.write(smiles.strip()+","+str(spacehastenid)+"\n")
        chunk_size += 1
        if chunk_size >= smiles_per_cpu:
            chunk_size = 0
            w.close()
            new_file = True
    if not new_file:
        w.close()
    
    slurm_functions.write_predict_slurm(cycle_dir+"/PREDICT",args)

    curdir = os.getcwd()
    os.chdir(cycle_dir+"/PREDICT")
    os.system("rm -f jobdone-"+args.name+"-CPU* predicted_predict_*.csv")
    print("Predicting docking scores via slurm at "+cycle_dir+"/PREDICT ...")
    os.system("sbatch submit_predict_"+args.name+"_cycle"+str(functions.get_latest_cycle(args.name))+".sh")
    os.chdir(curdir)

    jobs_left = args.cpu - len(glob.glob(cycle_dir+"/PREDICT/jobdone-"+args.name+"-CPU*"))
    while jobs_left>0:
        time.sleep(5)
        jobs_left = args.cpu - len(glob.glob(cycle_dir+"/PREDICT/jobdone-"+args.name+"-CPU*"))
    print("Reading in predictions...")
    prop_inputs = []
    for prop_input in glob.glob(cycle_dir+"/PREDICT/predicted_predict_*.csv"):
        prop_inputs.append(prop_input)
    pool = mp.Pool(mp.cpu_count())
    filtered_mol_runs_with_predictions = list(filter(None,list(tqdm.tqdm(pool.imap_unordered(functions.get_rdkit_properties,prop_inputs),mininterval=1,total=len(prop_inputs)))))

    docking_scores = {}
    
    for cpu in range(len(filtered_mol_runs_with_predictions)):
        rawmol_list,docking_score_list = filtered_mol_runs_with_predictions[cpu]
        for comp_index in range(len(rawmol_list)):
            title = rawmol_list[comp_index]
            docking_scores[title] = float(docking_score_list[comp_index])

    to_db = []

    for title in docking_scores:
        to_db.append((docking_scores[title],model_version,title))

    print("Updating predicted scores...")
    c.executemany("UPDATE data SET pred_score =?,pred_version=? WHERE spacehastenid = ?",to_db)
    conn.commit()
    conn.close()
    print("Done!")
