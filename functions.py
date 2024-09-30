# SpaceHASTEN: functions used by various parts of the program
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
from rdkit import Chem
from rdkit.Chem import RegistrationHash
import argparse
import shutil
import os
import glob
import pandas as pd
import sqlite3
import cfg

import sys

def get_rdkit_properties(csv_filename):
    """
    Read rdkit output

    :return molecules that passed the properties criteria
    """
    passed_rows = pd.read_csv(csv_filename)
    return (passed_rows["smilesid"].tolist(),passed_rows["docking_score"].tolist())

def get_latest_model(name):
    """
    Get the version number of latest model
    
    :param name: args.name usually
    :return int of the latest version, 0 if no models exist
    """
    dbname = name + ".dbsh"
    if not os.path.exists(dbname):
        raise SystemExit("Internal Error: SpaceHASTEN database ("+dbname+") missing when called get_latest_model()!!!")
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    model_version = c.execute("SELECT COUNT(*) FROM models").fetchall()[0][0]
    conn.close()
    return model_version



def get_latest_cycle(name):
    """
    Get the latest searching cycle that has been done
    
    :param name: usually args.name
    :return int of the latest cycle, 0 if no exists
    """
    cycles = glob.glob(os.getenv("HOME")+"/SPACEHASTEN/SIMSEARCH_"+name+"_cycle*")
    if len(cycles) == 0:
        cycle = 0
    else:
        ranked_cycles = []
        for t in cycles:
            ranked_cycles.append((int(t.split("_cycle")[-1]),t))
        cycle = sorted(ranked_cycles)[-1][0]
    return cycle

def get_latest_iteration(name):
    """
    Get the latest docking iteration that has been done
    
    :param name: usually args.name
    :return int of the latest iteration, 0 if no exists
    """
    iterations = glob.glob(os.getenv("HOME")+"/SPACEHASTEN/DOCKING_"+name+"_iter*")
    if len(iterations) == 0:
        iteration = 0
    else:
        ranked_iterations = []
        for t in iterations:
            ranked_iterations.append((int(t.split("_iter")[-1]),t))
        iteration = sorted(ranked_iterations)[-1][0]
    return iteration

def mol2hash(line):
    """
    Take a § limited line and generate RegistrationHash

    :param line: The line to be parsed
    :return: RegistrationHash + input-line (§ limited)
    """
    l = line.split("§")
    smiles = l[0]
    mol_id = l[1]
    docking_score = l[2]

    mol = Chem.MolFromSmiles(smiles)
    if mol == None:
        return None

    return RegistrationHash.GetMolLayers(mol)[RegistrationHash.HashLayer.TAUTOMER_HASH] + "§" + line

def cxsmi2smi(cxsmiles):
    """
    Take a cxsmiles (no ID) and return smiles

    :param cxsmiles: The cxsmiles to be parsed
    :return: SMILES
    """
    
    mol = Chem.MolFromSmiles(cxsmiles)
    if mol == None:
        return None
    return Chem.MolToSmiles(mol)
