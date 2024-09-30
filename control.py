# SpaceHASTEN: control parameters of retrieved compounds
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
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import RegistrationHash
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors
import gzip

r = open(sys.argv[2],"rt")
prop_mw_min = float(r.readline())
prop_mw_max = float(r.readline())
prop_slogp_min = float(r.readline())
prop_slogp_max = float(r.readline())
prop_hba_min = int(r.readline())
prop_hba_max = int(r.readline())
prop_hbd_min = int(r.readline())
prop_hbd_max = int(r.readline())
prop_rotbonds_min = int(r.readline())
prop_rotbonds_max = int(r.readline())
prop_tpsa_min = float(r.readline())
prop_tpsa_max = float(r.readline())
r.close()

w = open("propoutput_"+sys.argv[1].replace(".smi",".csv"),"wt")
w.write("smiles,rawmol\n")
for line in gzip.open(sys.argv[1],"rt"):
    smiles = line.strip().split("§")[0]
    mol = Chem.MolFromSmiles(smiles)
    if mol == None:
        continue
    molwt = Descriptors.MolWt(mol)
    if molwt < prop_mw_min or molwt > prop_mw_max:
        continue
    logp = Crippen.MolLogP(mol)
    if logp < prop_slogp_min or logp > prop_slogp_max:
        continue
    hba = rdMolDescriptors.CalcNumHBA(mol)
    if hba < prop_hba_min or hba > prop_hba_max:
        continue
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    if hbd < prop_hbd_min or hbd > prop_hbd_max:
        continue
    rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rotbonds < prop_rotbonds_min or rotbonds > prop_rotbonds_max:
        continue
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    if tpsa < prop_tpsa_min or tpsa > prop_tpsa_max:
        continue
    w.write(smiles + "," + RegistrationHash.GetMolLayers(mol)[RegistrationHash.HashLayer.TAUTOMER_HASH] + "§" + line)
w.close()    
