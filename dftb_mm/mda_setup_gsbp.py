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
charmm_submission_script = "/projectnb/cui-buchem/rosadche/ap/dftb/ap/qcharmm_plumed2_1cr_72hr.sh"

charmm_required = """
module load intel/2021.1
module load openmpi/3.1.4_intel-2021
module use /projectnb/cui-buchem/rosadche/development/plumed_make/lib/plumed/
module load plumed_2.7.4_reilly
"""


#===============================================================================
#                               FUNCTIONS
#===============================================================================
def setup_steps123(pdb, prot, mg, state):
    """
    droplet setup & gsbp region and restraints
    """
    global base_dir, charmm_exe, charmm_submission_script
    template_dir = base_dir + "templates/"
    files = ["step1_setup_qmmm_droplet.inp", "step2_gsbp_region_partitioning.inp", "step3_gsbp_restraints.inp"]
    
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
        
        outfile = file.split(".")[0] + ".out"
        to_do = "{} -i {} > {}".format(charmm_exe, file, outfile)
        os.popen(to_do).readlines()
    os.chdir(base_dir)


def setup_step4(pdb, prot, mg, state):
    """
    gsbp Mij matrix setup
    """
    global base_dir, charmm_exe, charmm_submission_script
    template_dir = base_dir + "templates/"
    files = ["step4_gsbp_Mij.inp"]
    
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
        
        to_do = "bash {} {}".format(charmm_submission_script, file)
        os.popen(to_do).readlines()
    os.chdir(base_dir)


def setup_step5(pdb, prot, mg, state):
    """
    gsbp Phi(x) setup
    """
    global base_dir, charmm_exe, charmm_submission_script
    template_dir = base_dir + "templates/"
    files = ["step5_gsbp_phix.inp"]
    
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
        
        to_do = "bash {} {}".format(charmm_submission_script, file)
        os.popen(to_do).readlines()
    os.chdir(base_dir)

    

#===============================================================================
#                               MAIN
#===============================================================================
pdbs = ["5kfg", "5kfh"]
prots = ["wt", "s113a"]
mgs = [2, 3]

for pdb in pdbs:
    if pdb == "5kfg": state="rs"
    if pdb == "5kfh": state="ps"
    
    for prot in prots:
        for mg in mgs:
            
            #setup_steps123(pdb, prot, mg, state)
            #setup_step4(pdb, prot, mg, state)
            setup_step5(pdb, prot, mg, state)
            
            
            