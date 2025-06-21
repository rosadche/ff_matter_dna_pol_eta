#===============================================================================
#                               IMPORTS
#===============================================================================
import os
import shutil
import MDAnalysis as mda


#===============================================================================
#                               PARAMETERS
#===============================================================================
base_dir        = "/projectnb/cui-buchem/rosadche/dnap/qm_mm/"
charmm_exe      = "/projectnb/cui-buchem/rosadche/development/charmm_make_rmax/bin/charmm"
charmm_submission_script = "/projectnb/cui-buchem/rosadche/dnap/qm_mm/qcharmm_plumed2_1cr_150hr.sh"

charmm_required = """
module load intel/2021.1
module load openmpi/3.1.4_intel-2021
module use /projectnb/cui-buchem/rosadche/development/plumed_make/lib/plumed/
module load plumed_2.7.4_reilly
"""


#===============================================================================
#                               FUNCTIONS
#===============================================================================
def setup_minimization(pdb, prot, mg, state):
    """
    initial minimization after setup
    """
    global base_dir, charmm_exe, charmm_submission_script, charges
    template_dir = base_dir + "templates/"
    files = ["minimization.inp"]
    
    working_dir = base_dir + "setup/" + "{}_{}_{}mg_{}/".format(pdb, prot, mg, state)
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    
    [ shutil.copy(template_dir + file, working_dir + file) for file in files ]
    
    os.chdir(working_dir)
    for file in files:
        os.popen('sed -i "s/SEDPDB/{}/g" {}'.format(pdb, file) ).readlines()
        os.popen('sed -i "s/SEDPROT/{}/g" {}'.format(prot, file) ).readlines()
        os.popen('sed -i "s/SEDMG/{}/g" {}'.format(mg, file) ).readlines()
        os.popen('sed -i "s/SEDSTATE/{}/g" {}'.format(state, file) ).readlines()
        
        charge = charges["{}_{}mg".format(state, mg)]
        os.popen('sed -i "s/SEDCHARGE/{}/g" {}'.format(charge, file) ).readlines()
        
        to_do = "bash {} {}".format(charmm_submission_script, file)
        os.popen(to_do).readlines()
    os.chdir(base_dir)


def equilibration_runs(pdb, prot, mg, state):
    """
    5 independent equilibration runs
    """
    global base_dir, charmm_exe, charmm_submission_script, charges
    template_dir = base_dir + "templates/"
    files = ["equilibration.inp"]
    
    sccdftb_dir = base_dir + "setup/" + "{}_{}_{}mg_{}/".format(pdb, prot, mg, state)
    
    indi_runs = 5
    
    for run in range(1, indi_runs+1):
        working_dir = base_dir + "md/" + "{}_{}_{}mg_{}/".format(pdb, prot, mg, state) + "run{}/".format(run)
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        
        [ shutil.copy(template_dir + file, working_dir + file) for file in files ]
        shutil.copy(sccdftb_dir + "sccdftb.dat",  working_dir + "sccdftb.dat")
        
        os.chdir(working_dir)
        for file in files:
            os.popen('sed -i "s/SEDPDB/{}/g" {}'.format(pdb, file) ).readlines()
            os.popen('sed -i "s/SEDPROT/{}/g" {}'.format(prot, file) ).readlines()
            os.popen('sed -i "s/SEDMG/{}/g" {}'.format(mg, file) ).readlines()
            os.popen('sed -i "s/SEDSTATE/{}/g" {}'.format(state, file) ).readlines()
            
            charge = charges["{}_{}mg".format(state, mg)]
            os.popen('sed -i "s/SEDCHARGE/{}/g" {}'.format(charge, file) ).readlines()
        
            to_do = "bash {} {}".format(charmm_submission_script, file)
            os.popen(to_do).readlines()
        os.chdir(base_dir)


def production_runs(pdb, prot, mg, state):
    """
    5 independent production runs, daisy chain for X runs
    """
    global base_dir, charmm_exe, charmm_submission_script, charges
    template_dir = base_dir + "templates/"
    
    submission_file = "submit_production.sh"
    files = ["production.inp", submission_file]
    
    indi_runs = 5
    daisy_chain = 37
    
    for run in range(1, indi_runs+1):
        working_dir = base_dir + "md/" + "{}_{}_{}mg_{}/".format(pdb, prot, mg, state) + "run{}/".format(run)
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        
        [ shutil.copy(template_dir + file, working_dir + file) for file in files ]
        #shutil.copy(sccdftb_dir + "sccdftb.dat",  working_dir + "sccdftb.dat")
        
        jobbase = "prod_{}_{}_{}mg_{}".format(pdb, prot, mg, state)
        
        os.chdir(working_dir)
        for file in files:
            os.popen('sed -i "s/SEDPDB/{}/g" {}'.format(pdb, file) ).readlines()
            os.popen('sed -i "s/SEDPROT/{}/g" {}'.format(prot, file) ).readlines()
            os.popen('sed -i "s/SEDMG/{}/g" {}'.format(mg, file) ).readlines()
            os.popen('sed -i "s/SEDSTATE/{}/g" {}'.format(state, file) ).readlines()
            charge = charges["{}_{}mg".format(state, mg)]
            os.popen('sed -i "s/SEDCHARGE/{}/g" {}'.format(charge, file) ).readlines()
            
            
            os.popen('sed -i "s/SEDRUN/{}/g" {}'.format(run, file) ).readlines()
            
            os.popen('sed -i "s/SEDJOBNAME/{}/g" {}'.format(jobbase, file) ).readlines()
            os.popen('sed -i "s|SEDWRKDIR|{}|g" {}'.format(working_dir, file) ).readlines()
            os.popen('sed -i "s/SEDTEMPLATE/{}/g" {}'.format("production.inp", file) ).readlines()
            os.popen('sed -i "s/SEDTOTALROUNDS/{}/g" {}'.format(daisy_chain, file) ).readlines()
        
        to_do = "qsub {}".format(submission_file)
        os.popen(to_do).readlines()
        os.chdir(base_dir)

    

#===============================================================================
#                               MAIN
#===============================================================================
pdbs = ["5kfg", "5kfh"]
prots = ["wt", "s113a"]
mgs = [2, 3]

charges = dict()
# S113 vs. WT have the same charge
# 2 vs. 3 Mg^(2+) gives difference of 2 charge
charges["rs_2mg"] = +0.0
charges["rs_3mg"] = +2.0
# ps different from rs as the rs have the 3'OH proton, missing from ps
charges["ps_2mg"] = -1.0
charges["ps_3mg"] = +1.0

for pdb in pdbs:
    if pdb == "5kfg": state="rs"
    if pdb == "5kfh": state="ps"
    
    for prot in prots:
        for mg in mgs:
            
            #setup_minimization(pdb, prot, mg, state)
            #equilibration_runs(pdb, prot, mg, state)
            production_runs(pdb, prot, mg, state)
            
            
            


