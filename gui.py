# SpaceHASTEN: graphical user interface (GUI)
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
import functions
import training_functions
import importseeds_functions
import export_functions
import docking_functions
import prediction_functions
import simsearch_functions
import cfg

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import simpledialog
from tkinter.ttk import Radiobutton
from PIL import Image, ImageTk
import threading
from types import SimpleNamespace
import os
import glob
import pandas as pd
import tqdm
import multiprocessing as mp
# NOTE: it is important to avoid all tkinter calls in the worker threads
import queue

class SpaceHASTENGUI(tk.Tk):

    # class attributes
    logo = None 
    gui_dbname = None
    gui_search_cycle = None
    gui_docking_iteration = None
    gui_model_version = None
    gui_cpu = None
    gui_spacename = None
    gui_cycle_mode = None
    gui_export_mode = None
    gui_cutoff = None
    frame_main = None
    frame_working = None
    frame_task_menu = None
    frame_export_menu = None
    c = None
    progressbar_working = None
    progressbar_value = None
    q = None
    button_vs = None

    # font definitions
    font_buttons = ("Arial",12)
    font_text = ("Arial",10)

    # color definitions
    color_background = "#fff"
    color_button_bg = "#fff"
    color_button_fg = "#000"
    color_button_active_bg = "#0f0"
    color_button_active_fg = "#fff"

    def __init__(self,*args,**kwargs):
        super().__init__(*args, **kwargs)

        self.c = cfg.SpaceHASTENConfiguration()
        self.q = queue.Queue()

        print("Running SpaceHASTEN at "+self.c.SPACEHASTEN_DIRECTORY)

        self.title("SpaceHASTEN "+str(self.c.SPACEHASTEN_VERSION))
        self.geometry("541x400")
        self.resizable(0,0)
        self.config(background=self.color_background)
      
        self.logo = ImageTk.PhotoImage(Image.open(self.c.SPACEHASTEN_DIRECTORY + "/spacehasten_logo.png"))       
        self.gui_dbname = tk.StringVar(self,"<not defined>")
        self.gui_cutoff = tk.DoubleVar(self,-10.0)
        self.gui_search_cycle = tk.IntVar(self,0)
        self.gui_docking_iteration = tk.IntVar(self,0)
        self.gui_model_version = tk.IntVar(self,0)
        self.gui_cpu = tk.IntVar(self,self.c.MAX_CORES)
        self.gui_spacename = tk.StringVar(self,self.c.SPACES_FILE_DEFAULT)
        self.gui_cycle_mode = tk.StringVar(self,"1")
        self.gui_export_mode = tk.IntVar(self,1)
        self.progressbar_value = tk.DoubleVar(self,0)
                
        self.build_main_frame()
        self.build_working_frame()
        self.build_task_menu()
        self.build_export_menu()

        self.frame_main.grid()

        self.perioidic_call()

    def perioidic_call(self):
        self.after(200,self.perioidic_call)
        self.check_queue()

    def check_queue(self):
        while self.q.qsize():
            try:
                msg = self.q.get_nowait()
                if msg == "Error":
                    messagebox.showerror(title="Error",message="Error occured, see shell window for details.")
                    self.frame_working.grid_forget()
                    self.frame_main.grid()
                if msg == "Done":
                    messagebox.showinfo(title="Note",message="Task done.")
                    self.frame_working.grid_forget()
                    self.frame_main.grid()    
                if msg == "DoneTaskmenu":
                    messagebox.showinfo(title="Note",message="Task done.")
                    self.frame_working.grid_forget()
                    self.frame_task_menu.grid()    
                if msg == "DoneExport":
                    messagebox.showinfo(title="Note",message="Task done.")
                    self.frame_working.grid_forget()
                    self.frame_export_menu.grid()    
                if msg.startswith("Percent:"):
                    self.progressbar_value.set(float(msg.split(":")[1]))
                if msg.startswith("UpdateModel:"):
                    self.gui_model_version.set(int(msg.split(":")[1]))
                if msg.startswith("UpdateIteration:"):
                    self.gui_docking_iteration.set(int(msg.split(":")[1]))
                if msg.startswith("UpdateCycle:"):
                    self.gui_search_cycle.set(int(msg.split(":")[1]))
                if msg.startswith("UpdateVS"):
                    self.button_vs.config(text="Continue virtual screening")
            except queue.Empty:
                pass

    def build_main_frame(self):
        self.frame_main = tk.Frame(self)        
        tk.Label(self.frame_main,image=self.logo).grid(row=0)
        tk.Button(self.frame_main,text="Pick Enamine REALSpace seeds",font=self.font_buttons,command=self.gui_pickseeds,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=1)
        tk.Button(self.frame_main,text="New job",font=self.font_buttons,command=self.gui_new,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=2)
        tk.Button(self.frame_main,text="Load existing job",font=self.font_buttons,command=self.gui_load,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=3)
        tk.Button(self.frame_main,text="Quit",font=self.font_buttons,command=quit,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=4)
        tk.Label(self.frame_main,text="").grid(row=5)
        tk.Label(self.frame_main,text="Developed by Tuomo Kalliokoski, Orion Pharma",font=self.font_text).grid(row=6)
        tk.Label(self.frame_main,text="Powered by SpaceLight, FTrees, chemprop, Ligprep, Glide and RDKit",font=self.font_text).grid(row=7)

    def build_working_frame(self):
        self.frame_working = tk.Frame(self)
        tk.Label(self.frame_working,image=self.logo).grid(row=0)
        self.progressbar_working = ttk.Progressbar(self.frame_working,variable=self.progressbar_value,orient=tk.HORIZONTAL,length=500,maximum=100,mode="determinate")
        self.progressbar_working.grid(row=1)
        tk.Label(self.frame_working,text="Working, see shell window for progress...",font=self.font_text).grid(row=2)

    def build_task_menu(self):
        self.frame_task_menu = tk.Frame(self)
        tk.Label(self.frame_task_menu,image=self.logo).grid(row=0,columnspan=4)        
        tk.Label(self.frame_task_menu,text="Similarity search cycles:",font=self.font_text).grid(row=1,column=0,columnspan=2)
        tk.Label(self.frame_task_menu,textvariable=self.gui_search_cycle,font=self.font_text).grid(row=1,column=2,columnspan=2)
        tk.Label(self.frame_task_menu,text="Docking iterations:",font=self.font_text).grid(row=2,column=0,columnspan=2)
        tk.Label(self.frame_task_menu,textvariable=self.gui_docking_iteration,font=self.font_text).grid(row=2,column=2,columnspan=2)
        tk.Label(self.frame_task_menu,text="Model versions:",font=self.font_text).grid(row=3,column=0,columnspan=2)
        tk.Label(self.frame_task_menu,textvariable=self.gui_model_version,font=self.font_text).grid(row=3,column=2,columnspan=2)
        self.button_vs = tk.Button(self.frame_task_menu,text="Start virtual screening",command=self.gui_virtual_screening,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg)
        self.button_vs.grid(row=4,columnspan=4)
        tk.Button(self.frame_task_menu,text="Export results",command=self.gui_export_menu,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=5,columnspan=4)
        ttk.Separator(self.frame_task_menu,orient=tk.HORIZONTAL).grid(row=6,columnspan=4,sticky="ew")
        tk.Label(self.frame_task_menu,text="Run tasks manually:",font=self.font_text).grid(row=7,column=0)
        tk.Button(self.frame_task_menu,text="Search",command=self.gui_similarity,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=7,column=1)
        tk.Button(self.frame_task_menu,text="Dock",command=self.gui_docking,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=7,column=2)
        tk.Button(self.frame_task_menu,text="Train",command=self.gui_train,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=7,column=3)
        ttk.Separator(self.frame_task_menu,orient=tk.HORIZONTAL).grid(row=8,columnspan=4,sticky="ew")
        tk.Button(self.frame_task_menu,text="Go back",command=self.gui_goto_main_menu,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=9,columnspan=4)
        tk.Button(self.frame_task_menu,text="Quit",command=quit,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=10,columnspan=4)

    def validate_float_value(self,P):
        if P == "" or P == "-":
            return True
        try:
            float(P)
        except ValueError:
            return False
        return True
        
    def build_export_menu(self):
        self.frame_export_menu = tk.Frame(self)
        tk.Label(self.frame_export_menu,image=self.logo,font=self.font_text).grid(row=0)
        options = {"2D (.csv)":1,"3D (.maegz)":2}
        c = 0
        for (txt,val) in options.items():
            c+=1
            Radiobutton(self.frame_export_menu,text=txt,variable=self.gui_export_mode,value=val).grid(row=c)

        ttk.Separator(self.frame_export_menu,orient=tk.HORIZONTAL).grid(row=c+1,columnspan=2,sticky="ew")
        tk.Label(self.frame_export_menu,text="Virtual hit cutoff:",font=self.font_text).grid(row=c+2)
        tk.Entry(self.frame_export_menu,textvariable=self.gui_cutoff,width=5,font=self.font_text,validate="all",validatecommand=(self.register(self.validate_float_value),"%P")).grid(row=c+3)
        ttk.Separator(self.frame_export_menu,orient=tk.HORIZONTAL).grid(row=c+4,columnspan=2,sticky="ew")
        tk.Button(self.frame_export_menu,text="Export data",command=self.gui_export,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=c+5)
        tk.Button(self.frame_export_menu,text="Go back",command=self.gui_goto_task_menu,font=self.font_buttons,bg=self.color_button_bg,fg=self.color_button_fg,activebackground=self.color_button_active_bg,activeforeground=self.color_button_active_fg).grid(row=c+6)

    def gui_pickseeds(self):
        self.q.put("Percent:0")
        self.frame_main.grid_forget()
        self.frame_working.grid()
        number_of_seeds = simpledialog.askinteger("Number of seeds to pick","How many seeds to pick?",initialvalue=self.c.ENAMINEREAL_SEEDS_COUNT,parent=self)
        if number_of_seeds is None:
            messagebox.showwarning(title="Note",message="Picking seeds cancelled")
            self.frame_working.grid_forget()
            self.frame_main.grid()
            return
        
        number_of_local_cores = simpledialog.askinteger("Number of local cores to use","How many local cores to use?",initialvalue=self.c.ENAMINEREAL_SEEDS_CPU,parent=self)
        if number_of_local_cores is None:
            messagebox.showwarning(title="Note",message="Picking seeds cancelled")
            self.frame_working.grid_forget()
            self.frame_main.grid()
            return
        
        seedsname = filedialog.asksaveasfilename(filetypes=[("SMILES","*.smi")],defaultextension=".smi",title="Save seeds as...")
        if len(seedsname)==0:
            messagebox.showwarning(title="Note",message="Picking seeds cancelled")
            self.frame_working.grid_forget()
            self.frame_main.grid()
            return
        if os.path.exists(seedsname):
            self.frame_working.grid_forget()
            self.frame_main.grid()
            messagebox.showerror(title="Error!",message="Refusing to overwrite existing file!")
            return
        
        worker = threading.Thread(target=self.gui_thread_pickseeds,args=(seedsname,number_of_seeds,number_of_local_cores,))
        worker.start()

    def gui_thread_pickseeds(self,seedsname,num_seeds,num_cores):
        self.q.put("Percent:0")
        print("Loading seeds from cxsmiles...")
        os.system("date")
        leadlike = pd.read_csv(self.c.ENAMINEREAL_SEEDS,compression="bz2",sep="\t")
        self.q.put("Percent:10")
        leadlike.drop(columns=["Type"],inplace=True)
        print("Saving and picking seeds...")
        self.q.put("Percent:30")
        cx_seeds = leadlike.sample(n=num_seeds)["smiles"].to_list()
        # this consumes a lot of RAM, so not a good idea to use lots of CPUs
        self.q.put("Percent:50")
        pool = mp.Pool(num_cores)
        seeds = list(filter(None,list(tqdm.tqdm(pool.imap_unordered(functions.cxsmi2smi,cx_seeds),mininterval=1,total=len(cx_seeds)))))
        self.q.put("Percent:70")
        w = open(seedsname,"wt")
        seed_counter = 0
        for seed in seeds:
            seed_counter += 1
            w.write(seed + " S" + str(seed_counter)+"\n")
        w.close()
        self.q.put("Percent:99.9")
        self.q.put("Done")
        os.system("date")

    def gui_ask_num(self,title_message,question_message,inivalue,cancel_message):
        ask_n_queries = simpledialog.askinteger(title_message,question_message,initialvalue=inivalue,parent=self)
        if ask_n_queries is None:
            messagebox.showwarning(title="Note",message=cancel_message)
            self.frame_working.grid_forget()
            self.frame_task_menu.grid()
        return ask_n_queries

    def gui_ask_space(self,cancel_message):
        space_name = filedialog.askopenfilename(filetypes=[("Chemical space (.space)","*.space")],defaultextension=".in",title="Choose chemical space...",initialdir=self.c.SPACES_DIR_DEFAULT,initialfile=self.c.SPACES_FILE_DEFAULT)
        if len(space_name)==0:
            messagebox.showwarning(title="Note",message=cancel_message)
            self.frame_working.grid_forget()
            self.frame_task_menu.grid()
        return space_name
                
    def gui_virtual_screening(self):
        self.q.put("Percent:0")
        self.frame_task_menu.grid_forget()
        self.frame_working.grid()

        yes = messagebox.askyesno("Confirm","This might take days to run. Are you sure?")
        if yes:
            space_name = self.gui_ask_space("Virtual screening cancelled")
            if space_name is None or len(space_name)==0: return
            self.gui_spacename.set(space_name)
        
            ask_n_queries = self.gui_ask_num("Number of SimSearch queries","How many SimSearch queries to run?",self.c.QUERIES_DEFAULT,"Virtual screening cancelled")
            if ask_n_queries is None: return
            ask_n_docked = self.gui_ask_num("Number of docked molecules","How many molecules to dock?",self.c.DOCKING_DEFAULT,"Virtual screening cancelled")
            if ask_n_docked is None: return
            ask_cpu_count_simsearch = self.gui_ask_num("SimSearch tasks","How many parallel SimSearch tasks ("+self.c.SLURM_CPU_COUNT_SEARCH+" core(s) each)?",250,"Virtual screening cancelled")
            if ask_cpu_count_simsearch is None: return
            ask_cpu_count_docking = self.gui_ask_num("Docking tasks","How many parallel docking tasks ("+self.c.SLURM_CPU_COUNT_DOCK+" core(s) each)?",250,"Virtual screening cancelled")
            if ask_cpu_count_docking is None: return

            job_args = SimpleNamespace()
            job_args.c = self.c
            job_args.q = self.q
            job_args.name = self.gui_dbname.get().split("/")[-1].split(".")[0]
            job_args.space = self.gui_spacename.get()
            job_args.exe_spacelight = self.c.EXE_SPACELIGHT_DEFAULT
            job_args.exe_ftrees = self.c.EXE_FTREES_DEFAULT
            job_args.rdkit_chunk = self.c.RDKIT_CHUNK_DEFAULT
            job_args.rdkit_cpu = self.c.RDKIT_CPU_DEFAULT      
            job_args.chemprop_chunk = self.c.CHEMPROP_CHUNK_DEFAULT
            job_args.chemprop_cpu = self.c.CHEMPROP_CPU_DEFAULT
            job_args.scratch = self.c.SCRATCH_DEFAULT
            job_args.nnn = self.c.NNN_DEFAULT
            job_args.sim_spacelight = self.c.SIM_SPACELIGHT_DEFAULT
            job_args.sim_ftrees = self.c.SIM_FTREES_DEFAULT
            worker = threading.Thread(target=self.gui_thread_virtual_screening,args=(job_args,ask_cpu_count_simsearch,ask_cpu_count_docking,ask_n_queries,ask_n_docked))
            worker.start()
        else:
            messagebox.showwarning(title="Note",message="Virtual screening cancelled")
            self.frame_working.grid_forget()
            self.frame_task_menu.grid()
    
    def gui_thread_virtual_screening(self,job_args,simsearch_cpu,docking_cpu,n_queries,n_docked):
        if functions.get_latest_iteration(job_args.name) == 0:
            print("Starting virtual screening...")
        else:
            print("Continuing virtual screening...")
            print("Training new model as the first step...")
            training_functions.train_new_model(job_args)
            self.q.put("UpdateModel:"+str(functions.get_latest_model(job_args.name)))
      
        job_args.top = n_queries
        job_args.cpu = simsearch_cpu
        job_args.use_predicted = False
        simsearch_functions.simsearch(job_args,do_not_update_gui=True)
        self.q.put("Percent:25.0")
        job_args.use_predicted = True
        simsearch_functions.simsearch(job_args,do_not_update_gui=True)
        self.q.put("Percent:50.0")
        job_args.use_predicted = True
        simsearch_functions.simsearch(job_args,do_not_update_gui=True)
        self.q.put("Percent:75.0")

        job_args.top = n_docked
        job_args.sff = -1.0
        job_args.cpu = docking_cpu
        docking_functions.dock(job_args,do_not_update_gui=True)
        self.q.put("UpdateIteration:"+str(functions.get_latest_iteration(job_args.name)))
        self.q.put("UpdateCycle:"+str(functions.get_latest_cycle(job_args.name)))
        self.q.put("UpdateVS")
        self.q.put("Percent:99.9")
        self.q.put("DoneTaskmenu")
        
    def gui_train(self):
        self.q.put("Percent:0.0")
        self.frame_task_menu.grid_forget()
        self.frame_working.grid()
        job_args = SimpleNamespace()
        job_args.c = self.c
        job_args.name = self.gui_dbname.get().split("/")[-1].split(".")[0]
        
        job_args.scratch = self.c.SCRATCH_DEFAULT
        job_args.chemprop_chunk = self.c.CHEMPROP_CHUNK_DEFAULT
        yes = messagebox.askyesno("Confirm","Training can take hours. Are you sure?")
        if yes:
            ask_cpu_count_pred = self.gui_ask_num("Number of prediction tasks","How many parallel prediction tasks (1 core each)?",self.c.CHEMPROP_CPU_DEFAULT,"Training cancelled")
            if ask_cpu_count_pred is None: return
            job_args.cpu = ask_cpu_count_pred
            worker = threading.Thread(target=self.gui_thread_train,args=(job_args,))
            worker.start()
        else:
            messagebox.showwarning(title="Note",message="Training cancelled")
            self.frame_working.grid_forget()
            self.frame_task_menu.grid()
        
    def gui_thread_train(self,job_args):
        self.q.put("Percent:10.0")
        training_functions.train_new_model(job_args)
        self.q.put("Percent:70.0")
        prediction_functions.update_predicted_scores(job_args)
        self.q.put("UpdateModel:"+str(functions.get_latest_model(job_args.name)))
        self.q.put("Percent:99.9")
        self.q.put("DoneTaskmenu")
            
    def gui_new(self):
        self.frame_main.grid_forget()
        self.frame_working.grid()
        csv_or_smiles = messagebox.askyesnocancel("Docked already?","Have you already docked the seeds?")
        if csv_or_smiles is None:
            messagebox.showwarning(title="Note",message="New job cancelled")
            self.frame_working.grid_forget()
            self.frame_main.grid()
            return
        glideinfile = filedialog.askopenfilename(filetypes=[("Glide docking settings","*.in")],title="Select glide docking file...")
        if len(glideinfile)==0:
            self.frame_working.grid_forget()
            self.frame_main.grid()
            messagebox.showwarning(title="Note",message="New job cancelled")
            return
        glidegridfile = filedialog.askopenfilename(filetypes=[("Glide docking grid","*.zip")],title="Select glide grid file...")
        if len(glidegridfile)==0:
            self.frame_working.grid_forget()
            self.frame_main.grid()
            messagebox.showwarning(title="Note",message="New job cancelled")
            return
        if csv_or_smiles == False:
            datafile = filedialog.askopenfilename(filetypes=[("SMILES","*.smi")],title="Load seeds for docking...")
        else:
            datafile = filedialog.askopenfilename(filetypes=[("Glide docking results","*.csv")],title="Load docked seeds...")

        if len(datafile)==0:
            self.frame_working.grid_forget()
            self.frame_main.grid()
            messagebox.showwarning(title="Note",message="New job cancelled")
            return
        dbname = filedialog.asksaveasfilename(filetypes=[("SpaceHASTEN files","*.dbsh")],defaultextension=".dbsh",title="Save new job as...")
        if len(dbname)==0:
            messagebox.showwarning(title="Note",message="New job cancelled")
            self.frame_working.grid_forget()
            self.frame_main.grid()
            return
        if os.path.exists(dbname):
            self.frame_working.grid_forget()
            self.frame_main.grid()
            messagebox.showerror(title="Error!",message="Refusing to overwrite existing search!")
            return
        try:
            testfile = open(dbname,"wt")
            testfile.close()
            os.system("rm -f "+dbname)
        except:
            if not os.access(dbname, os.W_OK):
                self.frame_working.grid_forget()
                self.frame_main.grid()
                messagebox.showerror(title="Error!",message="Cannot write to "+dbname)
                return
        if os.path.exists(dbname):
            self.frame_working.grid_forget()
            self.frame_main.grid()
            messagebox.showerror(title="Error!",message="Refusing to overwrite existing search!")
            return
        job_args = SimpleNamespace()
        job_args.c = self.c
        job_args.name = dbname.split("/")[-1].split(".")[0]
        job_args.input = datafile
        job_args.seeds_scorefield = self.c.FIELD_SCORE_DEFAULT
        job_args.seeds_title = self.c.FIELD_TITLE_DEFAULT
        job_args.seeds_smiles = self.c.FIELD_SMILES_DEFAULT
        job_args.rdkit_cpu = 0
        job_args.dock_param = glideinfile
        job_args.q = self.q
        job_args.glidegridfile = glidegridfile
          
        if csv_or_smiles == False:
            ask_cpu_count = simpledialog.askinteger("CPUs","How many CPUs to use?",initialvalue=250,parent=self)
            if ask_cpu_count is None:
                self.frame_working.grid_forget()
                self.frame_main.grid()
                messagebox.showwarning(title="Note",message="New job cancelled")
                return
            job_args.cpu = ask_cpu_count
            # dock all seeds
            job_args.top = -1
            job_args.scratch = self.c.SCRATCH_DEFAULT
            job_args.sff = -1.0
        
        worker = threading.Thread(target=importseeds_functions.import_seeds,args=(job_args,))
        worker.start()
        
    def gui_load(self):
        dbname = filedialog.askopenfilename(filetypes=[("SpaceHASTEN files","*.dbsh")],defaultextension=".dbsh",title="Load existing job...")
        if len(dbname)==0:
            messagebox.showwarning(title="Note",message="Loading cancelled")
            return
        self.gui_dbname.set(dbname)
        name = dbname.split("/")[-1].split(".")[0]
        os.chdir("/".join(dbname.split("/")[0:-1]))
        self.gui_search_cycle.set(functions.get_latest_cycle(name))
        self.gui_docking_iteration.set(functions.get_latest_iteration(name))
        self.gui_model_version.set(functions.get_latest_model(name))
        if self.gui_docking_iteration.get() > 0:
            self.button_vs.config(text="Continue virtual screening")
        self.frame_main.grid_forget()
        self.frame_task_menu.grid()
        self.title("SpaceHASTEN "+str(self.c.SPACEHASTEN_VERSION)+ " -- "+dbname.split("/")[-1])

    def gui_export_menu(self):
        self.frame_task_menu.grid_forget()
        self.frame_export_menu.grid()

    def gui_goto_main_menu(self):
        self.frame_task_menu.grid_forget()
        self.frame_main.grid()

    def gui_goto_task_menu(self):
        self.frame_export_menu.grid_forget()
        self.frame_task_menu.grid()

    def gui_export(self):
        self.q.put("Percent:0.0")
        self.frame_export_menu.grid_forget()
        self.frame_working.grid()

        job_args = SimpleNamespace()
        job_args.c = self.c
        job_args.q = self.q
        job_args.name = self.gui_dbname.get().split("/")[-1].split(".")[0]
        job_args.cutoff = float(self.gui_cutoff.get())
        job_args.scratch = self.c.SCRATCH_DEFAULT

        if int(self.gui_export_mode.get()) == 1:
            export_filename = filedialog.asksaveasfilename(filetypes=[("CSV files","*.csv")],defaultextension=".csv",title="Export virtual hits to CSV...",initialfile="spacehasten_virtualhits_cutoff"+str(job_args.cutoff).replace(".","_")+"_"+job_args.name+"_iter"+str(functions.get_latest_iteration(job_args.name))+".csv")
            if len(export_filename)==0:
                messagebox.showwarning(title="Note",message="Exporting cancelled")
                self.frame_working.grid_forget()
                self.frame_export_menu.grid()
                return
            job_args.resfilename = export_filename
            worker = threading.Thread(target=export_functions.export_results,args=(job_args,))
        else:
            export_filename = filedialog.asksaveasfilename(filetypes=[("maegz files","*.maegz")],defaultextension=".csv",title="Export poses of virtual hits to maegz...",initialfile="spacehasten_virtualhits_cutoff"+str(job_args.cutoff).replace(".","_")+"_"+job_args.name+".maegz")
            if len(export_filename)==0:
                messagebox.showwarning(title="Note",message="Exporting cancelled")
                self.frame_working.grid_forget()
                self.frame_export_menu.grid()
                return
            job_args.resfilename = export_filename
            worker = threading.Thread(target=self.gui_thread_export_poses,args=(job_args,))
        worker.start()
                    
    def gui_thread_export_poses(self,job_args):
        os.system("rm -f "+job_args.resfilename)
        for dock_dir in glob.glob(os.getenv("HOME")+"/SPACEHASTEN/DOCKING_"+job_args.name+"_iter*"):
            job_args.iteration = int(dock_dir.split("_iter")[-1])
            print("Exporting iteration",job_args.iteration)
            export_functions.export_poses(job_args)
        self.q.put("Percent:50")
        os.system("mv "+job_args.resfilename+" tmp.mae")
        os.system("pigz -c tmp.mae > "+job_args.resfilename)
        os.system("rm -f tmp.mae")
        self.q.put("Percent:99.9")
        self.q.put("DoneExport")

    def gui_docking(self):
        self.q.put("Percent:0")
        self.frame_task_menu.grid_forget()
        self.frame_working.grid()

        ask_n_docked = self.gui_ask_num("Number of docked molecules","How many molecules to dock?",self.c.DOCKING_DEFAULT,"Docking cancelled")
        if ask_n_docked is None: return
        ask_cpu_count_docking = self.gui_ask_num("Docking tasks","How many parallel docking tasks ("+self.c.SLURM_CPU_COUNT_DOCK+" core(s) each)?",250,"Docking cancelled")
        if ask_cpu_count_docking is None: return

        job_args = SimpleNamespace()
        job_args.c = self.c
        job_args.q = self.q
        job_args.name = self.gui_dbname.get().split("/")[-1].split(".")[0]
        job_args.space = self.gui_spacename.get()
        job_args.exe_spacelight = self.c.EXE_SPACELIGHT_DEFAULT
        job_args.exe_ftrees = self.c.EXE_FTREES_DEFAULT
        job_args.rdkit_chunk = self.c.RDKIT_CHUNK_DEFAULT
        job_args.rdkit_cpu = self.c.RDKIT_CPU_DEFAULT      
        job_args.chemprop_chunk = self.c.CHEMPROP_CHUNK_DEFAULT
        job_args.chemprop_cpu = self.c.CHEMPROP_CPU_DEFAULT
        job_args.scratch = self.c.SCRATCH_DEFAULT
        job_args.nnn = self.c.NNN_DEFAULT
        job_args.top = ask_n_docked
        job_args.sff = -1.0
        job_args.cpu = ask_cpu_count_docking
        
        worker = threading.Thread(target=self.gui_thread_docking,args=(job_args,))
        worker.start()
    
    def gui_thread_docking(self,job_args):
        docking_functions.dock(job_args,do_not_update_gui=True)
        self.q.put("UpdateIteration:"+str(functions.get_latest_iteration(job_args.name)))
        self.q.put("Percent:99.9")
        self.q.put("DoneTaskmenu")

    def gui_thread_similarity(self,job_args):
        simsearch_functions.simsearch(job_args,do_not_update_gui=True)
        self.q.put("UpdateCycle:"+str(functions.get_latest_cycle(job_args.name)))
        self.q.put("Percent:99.9")
        self.q.put("DoneTaskmenu")
        
    def gui_similarity(self):
        self.q.put("Percent:0")
        self.frame_task_menu.grid_forget()
        self.frame_working.grid()

        space_name = self.gui_ask_space("Similarity search cancelled")
        if space_name is None or len(space_name)==0: return
        self.gui_spacename.set(space_name)

        # ask if to pick new queries from docked or predicted
        docked_or_predicted = messagebox.askyesnocancel("Used docked molecules for queries?","Do you wish to use docked molecules for queries? If you answer No, queries will be picked based on predicted docking scores.")
        if docked_or_predicted is None:
            messagebox.showwarning(title="Note",message="Similarity search cancelled")
            self.frame_working.grid_forget()
            self.frame_task_menu.grid()
            return
        
        ask_n_queries = self.gui_ask_num("Number of SimSearch queries","How many SimSearch queries to run?",self.c.QUERIES_DEFAULT,"Virtual screening cancelled")
        if ask_n_queries is None: return
        ask_cpu_count_simsearch = self.gui_ask_num("SimSearch tasks","How many parallel SimSearch tasks ("+self.c.SLURM_CPU_COUNT_SEARCH+" core(s) each)?",250,"Virtual screening cancelled")
        if ask_cpu_count_simsearch is None: return
        
        job_args = SimpleNamespace()
        job_args.c = self.c
        job_args.q = self.q
        job_args.name = self.gui_dbname.get().split("/")[-1].split(".")[0]
        job_args.space = self.gui_spacename.get()
        job_args.exe_spacelight = self.c.EXE_SPACELIGHT_DEFAULT
        job_args.exe_ftrees = self.c.EXE_FTREES_DEFAULT
        job_args.rdkit_chunk = self.c.RDKIT_CHUNK_DEFAULT
        job_args.rdkit_cpu = self.c.RDKIT_CPU_DEFAULT      
        job_args.chemprop_chunk = self.c.CHEMPROP_CHUNK_DEFAULT
        job_args.chemprop_cpu = self.c.CHEMPROP_CPU_DEFAULT
        job_args.scratch = self.c.SCRATCH_DEFAULT
        job_args.nnn = self.c.NNN_DEFAULT
        job_args.sim_spacelight = self.c.SIM_SPACELIGHT_DEFAULT
        job_args.sim_ftrees = self.c.SIM_FTREES_DEFAULT
        job_args.use_predicted = not docked_or_predicted
        job_args.top = ask_n_queries
        job_args.cpu = ask_cpu_count_simsearch
        
        worker = threading.Thread(target=self.gui_thread_similarity,args=(job_args,))
        worker.start()