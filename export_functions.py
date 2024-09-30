# SpaceHASTEN: functions to do the exporting
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
import tqdm
import multiprocessing as mp
import glob

def export_results(args):
    """
    Export results

    :param args: the args
    """
    dbname = args.name+".dbsh"
    
    print("Exporting results:")
    args.q.put("Percent:50.0")
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    w = open(args.resfilename,"wt")
    w.write("smiles,smilesid,dock_score,pred_score,spacelight,ftrees,dock_iteration\n")
    #for smiles,smilesid,dock_score,pred_score,spacehasten,ftrees,dock_iteration in c.execute("SELECT smiles,smilesid,dock_score,pred_score,spacelight,ftrees,dock_iteration FROM data WHERE dock_score <= "+str(args.cutoff)+" AND (spacelight IS NOT NULL OR ftrees is NOT NULL) ORDER BY dock_score"):
    for smiles,smilesid,dock_score,pred_score,spacehasten,ftrees,dock_iteration in c.execute("SELECT smiles,smilesid,dock_score,pred_score,spacelight,ftrees,dock_iteration FROM data WHERE dock_score <= "+str(args.cutoff)+" ORDER BY dock_score"):
        w.write(smiles.strip()+","+smilesid.strip()+","+str(dock_score)+","+str(pred_score)+","+str(spacehasten)+","+str(ftrees)+","+str(dock_iteration)+"\n")
    w.close()
    conn.close()
    print("CSV created:",args.resfilename)
    args.q.put("Percent:99.9")
    args.q.put("DoneExport")

def export_poses(args):
    """
    Export poses using Schrodinger toolkit

    :param args: the args
    """
    dock_dir = os.getenv("HOME")+"/SPACEHASTEN/DOCKING_"+args.name+"_iter"+str(args.iteration)
    
    print("Exporting results:")
    print("\nDecompressing results to a local drive...")
    resdir=args.scratch+"/"+os.getenv("USER")+"/COLLECT_"+args.name+"_iter"+str(args.iteration)
    os.system("rm -fr "+resdir)
    os.system("mkdir -p "+resdir)
    pool = mp.Pool(mp.cpu_count())
    result_tars = []
    for result_tar in glob.glob(dock_dir+"/results-*.tar.gz"):
        result_tars.append("tar xzf "+result_tar+" -C "+resdir)
    jobs = list(tqdm.tqdm(pool.imap_unordered(os.system,result_tars),mininterval=1,total=len(result_tars)))

    print("Extracting poses...")
    result_pvs = []
    for result_pv in glob.glob(resdir+"/*_pv.maegz"):
        result_pvs.append(args.c.EXPORTPOSES_EXE + " "+ result_pv + " " + str(args.cutoff)+" "+args.name+".dbsh")
    jobs = list(tqdm.tqdm(pool.imap_unordered(os.system,result_pvs),mininterval=1,total=len(result_pvs)))
    os.system("cat "+resdir+"/spacehasten_virtual_hits_*.mae >> "+args.resfilename)
    os.system("rm -fr "+resdir)
    print("Poses exported to",args.resfilename)
    

