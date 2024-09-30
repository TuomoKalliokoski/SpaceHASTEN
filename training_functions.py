# SpaceHASTEN: functions to do the training
#
import os
import functions
import sqlite3
import pandas as pd
import slurm_functions
import time
import glob

def train_new_model(args):
    """
    Train new model

    :param args: the args
    """
    print("Training new model:")
    dbname = args.name + ".dbsh"
    
    model_version = functions.get_latest_model(args.name)+1
    print("Model version",model_version)

    modeldir=os.getenv("HOME")+"/SPACEHASTEN/TRAIN_"+args.name+"_ver"+str(model_version)
    os.system("rm -fr "+modeldir)
    os.system("mkdir -p "+modeldir)

    print("Extracting data for training...")
    data_filename=modeldir+"/train_"+args.name+"_ver"+str(model_version)+".csv"
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    # exclude failed dockings with arbitrary score as they screw up training
    training_data = pd.read_sql_query("SELECT smiles,dock_score FROM data WHERE dock_score IS NOT NULL AND dock_score < "+str(args.c.TRAIN_DOCKING_CUTOFF),conn)
    training_data["docking_score"] = training_data["dock_score"]
    training_data.drop(["dock_score"],axis=1,inplace=True)
    training_data.to_csv(data_filename,index=False)
    slurm_functions.write_train_slurm(modeldir,args)

    curdir=os.getcwd()
    os.chdir(modeldir)
    print("Running training via slurm at "+modeldir+" ...")
    os.system("sbatch submit_train_"+args.name+"_ver"+str(model_version)+".sh")
    jobs_left = 1 - len(glob.glob(modeldir+"/jobdone-train_"+args.name+"*"))
    while jobs_left>0:
        time.sleep(5)
        jobs_left = 1 - len(glob.glob(modeldir+"/jobdone-train_"+args.name+"*"))
    #print("Running chemprop...")
    #os.system("chemprop_train --dataset_type regression --target_columns docking_score --data_path "+data_filename+" --save_dir model_"+args.name+"_ver"+str(model_version)+" --batch_size 250 --no_cache_mol")
    os.system("rm "+data_filename)
    tarfile = "model_"+args.name+"_ver"+str(model_version)+".tar.gz"
    os.system("tar -czf "+tarfile+" model_"+args.name+"_ver"+str(model_version))
    with open(tarfile,"rb") as f:
        blob = f.read()
    c.execute("INSERT INTO models VALUES("+str(model_version)+",?)",[memoryview(blob)])
    conn.commit()
    conn.close()
    os.system("rm -r model_"+args.name+"_ver"+str(model_version)+" model_"+args.name+"_ver"+str(model_version)+".tar.gz")
    os.chdir(curdir)
    print("Training complete!")
