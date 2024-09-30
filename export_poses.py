# SpaceHASTEN: export poses using Schrödinger Suite
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
import sqlite3

try:
    from schrodinger import structure
except:
    raise SystemExit("Schrödinger Suite Python environment not detected, start with $SCHRODINGER/run!")

name = sys.argv[1]
cutoff = float(sys.argv[2])

resfilename = "/".join(sys.argv[1].split("/")[0:-1])+"/spacehasten_virtual_hits_"+sys.argv[1].split("/")[-1].replace(".maegz",".mae")

conn = sqlite3.connect(sys.argv[3])
c = conn.cursor()

writer = structure.StructureWriter(resfilename)
hits = []
hits_found = set()
first = True
for st in structure.StructureReader(name):
    if first:
        first = False
        continue
    docking_score = float(st.property["r_i_docking_score"])
    hit_id = int(st.property["s_m_title"])
    if docking_score <= cutoff:
        compound_id = c.execute("SELECT smilesid FROM data WHERE spacehastenid = ?",[hit_id]).fetchone()[0]
        st.property["s_m_title"] = compound_id
        writer.append(st)
    else:
        # pv files are sorted by docking_score
        break
conn.close()
writer.close()
