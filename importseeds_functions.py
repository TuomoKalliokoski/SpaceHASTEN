# SpaceHASTEN: functions to import seeds
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
import pandas as pd
import multiprocessing as mp
import tqdm
import functions
import sqlite3
import training_functions
import docking_functions
import numpy as np

def add_param_to_dbsh(c,glideinfile):
    with open(glideinfile,"rb") as f:
        blob = f.read()
    c.execute("INSERT INTO docking_param VALUES(?)",[memoryview(blob)])

def add_grid_to_dbsh(c,glidegridfile):
    with open(glidegridfile,"rb") as f:
        blob = f.read()
    c.execute("INSERT INTO docking_grid VALUES(?)",[memoryview(blob)])

def import_seeds(args):
    """
    Import seed docking results / dock seeds for SpaceHASTEN

    :param args: the args
    
    """
    args.q.put("Percent:0")
    os.system("date")
    print("Importing seeds:")
    dbname = args.name+".dbsh"

    print("1. Creating new database...")
    args.q.put("Percent:10")
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    c.execute("CREATE TABLE data (spacehastenid INTEGER PRIMARY KEY,reghash TEXT,smiles TEXT,smilesid TEXT,dock_score REAL,pred_score REAL,spacelight REAL,ftrees REAL,query INTEGER,dock_iteration INTEGER,pred_version INTEGER,simsearch_cycle INTEGER)")
    c.execute("CREATE TABLE docking_param (dock_param BLOB)")
    c.execute("CREATE TABLE docking_grid (dock_grid BLOB)")
    args.q.put("Percent:10")
    print("2. Importing docking parameters and grid...")
    add_param_to_dbsh(c,args.dock_param)
    add_grid_to_dbsh(c,args.glidegridfile)
    c.execute("CREATE TABLE models (model_version INTEGER UNIQUE,model_tar BLOB)")
    
    args.q.put("Percent:20")
    print("3. Loading data...")
    if args.input.endswith(".csv"):
        print("From CSV")
        data = pd.read_csv(args.input)
    else:
        print("From SMILES")
        data = pd.read_csv(args.input,sep=" ",header=None,names=[args.c.FIELD_SMILES_DEFAULT,args.c.FIELD_TITLE_DEFAULT])
        data[args.c.FIELD_SCORE_DEFAULT] = "NULL"
                                    
    if args.c.FIELD_SCORE_DEFAULT not in data.columns:
        print("Error: --seeds_scorefield ("+args.c.FIELD_SCORE_DEFAULT+") not found in seeds file --input ("+args.input+")")
        args.q.put("Error")
        return False
    if args.c.FIELD_TITLE_DEFAULT not in data.columns:
        print("Error: --seeds_title ("+args.c.FIELD_TITLE_DEFAULT+") not found in seeds file --input ("+args.input+")")
        args.q.put("Error")
        return False
    if args.c.FIELD_SMILES_DEFAULT not in data.columns:
        print("Error: --seeds_smiles ("+args.c.FIELD_SMILES_DEFAULT+") not found in seeds file --input ("+args.input+")")
        args.q.put("Error")
        return False
    data[args.c.FIELD_TITLE_DEFAULT] = data[args.c.FIELD_TITLE_DEFAULT].astype(str)
    data[args.c.FIELD_SCORE_DEFAULT] = data[args.c.FIELD_SCORE_DEFAULT].astype(str)
    raw_data = data[[args.c.FIELD_SMILES_DEFAULT,args.c.FIELD_TITLE_DEFAULT,args.c.FIELD_SCORE_DEFAULT]].stack().groupby(level=0).agg("§".join).values.tolist()
    args.q.put("Percent:33")
    print("4. Processing molecules...")
    if args.rdkit_cpu > 0:
        cores = args.rdkit_cpu
    else:
        cores = mp.cpu_count()
    pool = mp.Pool(cores)
    mols = list(filter(None,list(tqdm.tqdm(pool.imap_unordered(functions.mol2hash,raw_data),mininterval=1,total=len(raw_data)))))

    to_db = []
    print("5. Importing molecules...")
    for mol in mols:
        if args.input.endswith(".csv"):
            l = mol.split("§")[0:4]
        else:
            l = mol.split("§")[0:3]
        to_db.append(l)
    if args.input.endswith(".csv"):
        c.executemany("INSERT INTO data(reghash,smiles,smilesid,dock_score,dock_iteration) VALUES (?,?,?,?,0)",to_db)
    else:
        c.executemany("INSERT INTO data(reghash,smiles,smilesid) VALUES (?,?,?)",to_db)
    conn.commit()
    print("6. Indexing SQLite3...")
    c.execute("CREATE INDEX idx_reghash ON data(reghash)")
    conn.commit()
    conn.close()
    
    args.q.put("Percent:50")
    if not args.input.endswith(".csv"):    
        docking_functions.dock(args,importing_seeds=True,do_not_update_gui=True)

    args.q.put("Percent:70")
    training_functions.train_new_model(args)
    
    print("Import done!")
    os.system("date")
    args.q.put("Percent:99.9")
    args.q.put("Done")
    return True
