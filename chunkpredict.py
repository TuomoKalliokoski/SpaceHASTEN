# SpaceHASTEN: run chemprop prediction in chunks
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
import sys
import pandas as pd

def pred_chunk(filename,model):
    os.system("rm -f preds.csv")
    os.system("OMP_NUM_THREADS=1 chemprop_predict --no_cache_mol --no_cuda --num_workers 0 --test_path "+filename+" --checkpoint_dir "+model+" --preds_path preds.csv >/dev/null 2>/dev/null")
    predictions = pd.read_csv("preds.csv")
    os.system("rm -f preds.csv")
    predictions.drop(columns=["smiles"],inplace=True)
    return predictions

chunk = []
chunk_max = int(sys.argv[1])
chunk_size = 0
new_file = True
os.system("rm -f chunk_"+sys.argv[2].split(".")[0]+".csv")
all_preds = None
for line in open(sys.argv[2]).readlines()[1:]:
    l = line.strip().split(",")
    smiles = l[0]
    title = l[1]
    if new_file:
        csv_filename = "chunk_"+sys.argv[2].split(".")[0]+".csv"
        w = open(csv_filename,"wt")
        w.write("smiles,smilesid\n")
        new_file = False
    w.write(smiles.strip()+","+title+"\n")
    chunk_size += 1
    if chunk_size >= chunk_max:
        chunk_size = 0 
        w.close()
        new_file = True
        preds = pred_chunk("chunk_"+sys.argv[2].split(".")[0]+".csv",sys.argv[3])
        if all_preds is not None:
            all_preds = pd.concat([all_preds,preds])
        else:
            all_preds = preds
if not new_file:
    w.close()
    preds = pred_chunk("chunk_"+sys.argv[2].split(".")[0]+".csv",sys.argv[3])
    all_preds = pd.concat([all_preds,preds])

all_preds.to_csv("predicted_"+sys.argv[2].split(".")[0]+".csv",index=False)
