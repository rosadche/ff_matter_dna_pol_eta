#==============================================================================
#                                   IMPORTS
#==============================================================================
import os
#----------------------------------
import numpy as np
#----------------------------------
import MDAnalysis as mda
from MDAnalysis.analysis import align, nuclinfo, distances


#==============================================================================
#                                   FUNCTIONS
#==============================================================================
def create_G09CMD(direc, details, charge, spin):
    text = ""
    text += "$method" + "\n"
    text += details["method"] + "\n" + "\n"
    
    text += "$basis" + "\n"
    text += details["basis"] + "\n" + "\n"
    
    text += "$title" + "\n"
    text += details["title"] + "\n" + "\n"
    
    text += "$charge" + "\n"
    text += str(charge) + "\n" + "\n"
    
    text += "$spin" + "\n"
    text += str(spin) + "\n" + "\n"
    
    text += "$pop" + "\n"
    text += details["pop"] + "\n" + "\n"
    
    text += "$extra" + "\n"
    text += details["extra"] + "\n" + "\n"
    
    text += "$nproc" + "\n"
    text += str(details["nproc"]) + "\n" + "\n"
    
    text += "$mem" + "\n"
    text += details["mem"] + "\n" + "\n"
    
    text += "$extend" + "\n"
    text += details["extend"] + "\n" + "\n"
    
    with open(direc + details["g09cmd"], "w") as outfile:
        outfile.write(text)
    

def create_G09PROFILE(direc, details):
    text = "\n"
    with open(direc + details["g09profile"], "w") as outfile:
        outfile.write(text)
    
    
    
#==============================================================================
#                                  QM FILES
#==============================================================================    
base_dir = "/projectnb/cui-buchem/rosadche/dnap/sugar_dft/"
template = base_dir + "template_g09_sp_of_traj_no_scr.inp"
toppar   = base_dir + "toppar.str"
submit_script = base_dir + "qcharmm_g09_28cr_15day.sh"


#------------- molecules dict setup --------------------------------------------
molecules = dict()

"""
molecules["THFO"] = dict()
molecules["THFO"]["dir"] = 'thfo_neutral'
molecules["THFO"]["psf"] = "ligandrm.psf"
molecules["THFO"]["cor"] = "ligandrm.crd"
molecules["THFO"]["charge"] = 0
molecules["THFO"]["multiplicity"] = 1
molecules["THFO"]["proton"] = "name H32' "

molecules["RIBE"] = dict()
molecules["RIBE"]["dir"] = "ribe_neutral"
molecules["RIBE"]["psf"] = "ligandrm.psf"
molecules["RIBE"]["cor"] = "ligandrm.crd"
molecules["RIBE"]["charge"] = 0
molecules["RIBE"]["multiplicity"] = 1
molecules["RIBE"]["proton"] = "name H4"


# set 2

molecules["THFT"] = dict()
molecules["THFT"]["dir"] = 'thft_neutral'
molecules["THFT"]["psf"] = "ligandrm.psf"
molecules["THFT"]["cor"] = "ligandrm.crd"
molecules["THFT"]["charge"] = 0
molecules["THFT"]["multiplicity"] = 1
molecules["THFT"]["proton"] = "name H8"

molecules["RIBET"] = dict()
molecules["RIBET"]["dir"] = "ribet_neutral"
molecules["RIBET"]["psf"] = "ligandrm.psf"
molecules["RIBET"]["cor"] = "ligandrm.crd"
molecules["RIBET"]["charge"] = 0
molecules["RIBET"]["multiplicity"] = 1
molecules["RIBET"]["proton"] = "name H7"


# set 3 -------------------------------
molecules["DTHYNS"] = dict()
molecules["DTHYNS"]["dir"] = 'dthyns_neutral'
molecules["DTHYNS"]["charge"] = 0
molecules["DTHYNS"]["multiplicity"] = 1

molecules["RTHYNS"] = dict()
molecules["RTHYNS"]["dir"] = 'rthyns_neutral'
molecules["RTHYNS"]["charge"] = 0
molecules["RTHYNS"]["multiplicity"] = 1
"""
# set 4 -------------------------------
"""
molecules["DTHYMP"] = dict()
molecules["DTHYMP"]["dir"] = 'dthymp_neutral'
molecules["DTHYMP"]["charge"] = 0
molecules["DTHYMP"]["multiplicity"] = 1
"""

molecules["RTHYMP"] = dict()
molecules["RTHYMP"]["dir"] = 'rthymp_neutral'
molecules["RTHYMP"]["charge"] = 0
molecules["RTHYMP"]["multiplicity"] = 1


todo_lst = ["neutral", "deprot"]

correctly_labeled = ["THFO", "DTHYNS", "RTHYNS"]
pucker_atom_template = ["ATOM_O4", "ATOM_C1", "ATOM_C2", "ATOM_C3", "ATOM_C4"]
pucker_atoms_1       = ["O4'", "C1'", "C2'", "C3'", "C4'"]
pucker_atoms_2       = ["O1", "C4", "C3", "C2", "C1"]



details = dict()
details["g09cmd"]     = "g09cmd_mp2.txt"
details["g09profile"] = "g09profile.txt"
details["g09inp"]     = "g09_mp2"

details["method"]   = "mp2"
details["basis"]    = "6-311++G(d,p)"
details["extra"]    = ""

details["title"]    = "Title: gaussian-mp2"
details["nproc"]    = 28
details["mem"]      = "28GB"

details["pop"]      = ""
details["extend"]   = ""


#==============================================================================
#                                   Setup files
#==============================================================================
for mol in molecules.keys():
    for case in todo_lst:
        
        spin_multiplicity = 1
        
        scan_dir  = base_dir + mol + "_" + case + "_scans" + "/"
        if not os.path.exists(scan_dir):
            os.makedir(scan_dir)
        
        scan_file_name = mol + "_" + case + "_mp2_sp_scan.inp"
        scan_file = scan_dir + scan_file_name
        
        if case == "neutral":
            charge = 0
            psf = "ligandrm.psf"
            cor = "ligandrm.crd"
        else:
            charge = -1
            psf = "ligandrm_deprot.psf"
            cor = "ligandrm_deprot.crd"
            
        u = mda.Universe(base_dir + molecules[mol]["dir"] + "/" + psf, base_dir + molecules[mol]["dir"] + "/" + cor, topology_format="PSF", format="CRD")
        
        os.popen('cp {} {}'.format(template, scan_file)).readlines()
        os.popen('cp {} {}/toppar.str'.format(toppar, scan_dir)).readlines()
        
        os.popen('sed -i "s/{}/{}/g" {}'.format("MOLECULE", mol, scan_file)).readlines()
        os.popen('sed -i "s/{}/{}/g" {}'.format("PROTONATION", case, scan_file)).readlines()
        os.popen('sed -i "s/{}/{}/g" {}'.format("WHICHCORFILE", cor, scan_file)).readlines()
        os.popen('sed -i "s/{}/{}/g" {}'.format("WHICHPSFFILE", psf, scan_file)).readlines()
        os.popen('sed -i "s/{}/{}/g" {}'.format("QMCHARGE", charge, scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("INPUTDIR", base_dir + molecules[mol]["dir"], scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("G09CMDFILE", details["g09cmd"], scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("G09CMDPROFILEVAL", details["g09profile"], scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("G09INPUTFILE", details["g09inp"], scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("WHICHMETHOD", details["method"], scan_file)).readlines()
        
        
        for i in range(len(pucker_atom_template)):
            old = pucker_atom_template[i]
            if mol in correctly_labeled:
                new = pucker_atoms_1[i]
            else:
                new = pucker_atoms_2[i]
            
            replace = "{} 1 {}".format(mol, new)
            os.popen('sed -i "s/{}/{}/g" {}'.format(old, replace, scan_file)).readlines()
        
        create_G09CMD(scan_dir, details, charge, spin_multiplicity)
        create_G09PROFILE(scan_dir, details)
           
        os.chdir(scan_dir)
        os.popen('bash {} {}'.format(submit_script, scan_file_name)).readlines()



    
