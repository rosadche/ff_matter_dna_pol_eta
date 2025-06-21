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
def get_dftb_details(u, scan_dir):
    #-------------
    # SETUP SCCDFTB.DAT FILE BASED ON QM REGION
    #-------------
    dftb_par_dir = "/projectnb/cui-buchem/rosadche/ap/dftb/ap/" + "dftb_par"
    dft_par_phos_dir = dftb_par_dir + "/l02" # Improved phosphate hydolysis parameters
    
    # Zeta parameter
    ZETA = 4.00
    
    #hubbard paramters
    #----------------------
    HUBBARD = dict()
    #----------------------
    HUBBARD["C"]  = "-0.1492"
    HUBBARD["H"]  = "-0.1857"
    HUBBARD["N"]  = "-0.1535"
    HUBBARD["O"]  = "-0.1575"
    HUBBARD["S"]  = "-0.1100"
    HUBBARD["P"]  = "-0.1400"
    #----------------------
    # HALOGENS
    #----------------------
    HUBBARD["F"]  = "-0.1623"
    HUBBARD["CL"] = "-0.0697"
    HUBBARD["BR"] = "-0.0573"
    HUBBARD["I"]  = "-0.0433"
    HUBBARD["CU"] = "-0.1400"
    #----------------------
    # METALS
    #----------------------
    HUBBARD["MG"] = "-0.02"
    HUBBARD["ZN"] = "-0.03"
    
    
    # Guess the atom types for the structure
    # always check the output before using this blindly
    u_element_guess = mda.topology.guessers.guess_types( u.atoms.names )
    # we know for a fact the QQH are guessed as "Q" so lets replace these
    u_element_guess = [x if x!="Q" else "H" for x in u_element_guess]
    #add the guessed elements to the traj
    u.add_TopologyAttr('elements', u_element_guess)
    
    all_atoms = u.select_atoms("all")
    unique_atom_types = sorted( list( set(all_atoms.elements) ) )
    print("Unique QM atom types", unique_atom_types)
    
    #start by writing spl file location
    text = ""
    for atom1 in unique_atom_types:
        for atom2 in unique_atom_types:
            directory = dftb_par_dir
            file = "{}/{}{}.spl".format(directory, atom1.lower(), atom2.lower())
            if os.path.exists(file) == False:
                file = "{}/{}{}.spl".format(directory, atom2.lower(), atom1.lower())
            else:
                if os.path.exists(file) == False:
                    print("Error, {}{}.spl doesn't appear to exist".format(atom2.lower(), atom1.lower()))
                    exit()
            
            text += "\'" + file + "\'\n"
    
    #then by writing hubbard paramters
    for atom1 in unique_atom_types:
        text += "\'" + atom1.upper() + "\'  " + "{}\n".format(HUBBARD[atom1.upper()])
      
    #then write zeta for HBON correction, which you probably should be using  
    text += "{}\n".format(ZETA)
    
    with open(scan_dir + "sccdftb.dat", "w") as file:
        file.write(text)
    
    i = 1
    with open(scan_dir + "sccdftb_wmain.str", "w") as file:
        file.write("* Fake title\n*\n")
        for atom1 in unique_atom_types:
            text = "scalar WMAIN set {:.1f} sele (qm_all) .and. type {}*  SHOW end\n"
            if atom1 == "H":
                file.write( text.format(i, "QQ") )
            file.write( text.format(i, atom1) )
            i += 1
    
#==============================================================================
#                                  QM FILES
#==============================================================================    
base_dir = "/projectnb/cui-buchem/rosadche/dnap/sugar_dft/"
template = base_dir + "template_dftb_scan.inp"
toppar   = base_dir + "toppar.str"
submit_script = base_dir + "qcharmm_plumed2_1cr_24hr.sh"


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
molecules["DTHYMP"] = dict()
molecules["DTHYMP"]["dir"] = 'dthymp_neutral'
molecules["DTHYMP"]["charge"] = 0
molecules["DTHYMP"]["multiplicity"] = 1

molecules["RTHYMP"] = dict()
molecules["RTHYMP"]["dir"] = 'rthymp_neutral'
molecules["RTHYMP"]["charge"] = 0
molecules["RTHYMP"]["multiplicity"] = 1


# set 3 -------------------------------
todo_lst = ["neutral", "deprot"]

correctly_labeled = ["THFO", "RTHYNS", "DTHYNS"]
reverse_labeled   = ["THFT", "RIBE", "RIBET", "DTHYMP", "RTHYMP"]

pucker_atom_template = ["ATOM_O4", "ATOM_C1", "ATOM_C2", "ATOM_C3", "ATOM_C4"]
pucker_atoms_correct      = ["O4'", "C1'", "C2'", "C3'", "C4'"]
pucker_atoms_reverse      = ["O1", "C4", "C3", "C2", "C1"]

#==============================================================================
#                                   Setup files
#==============================================================================
for mol in molecules.keys():
    for case in todo_lst:
        
        scan_dir  = base_dir + mol + "_" + case + "_scans" + "/"
        if not os.path.exists(scan_dir):
            os.makedirs(scan_dir)
        
        scan_file_name = mol + "_" + case + "_dftb_scan.inp"
        scan_file = scan_dir + scan_file_name
        
        if case == "neutral":
            charge = molecules[mol]["charge"] 
            psf = "ligandrm.psf"
            cor = "ligandrm.crd"
        else:
            charge = molecules[mol]["charge"] - 1
            psf = "ligandrm_deprot.psf"
            cor = "ligandrm_deprot.crd"
            
        u = mda.Universe(base_dir + molecules[mol]["dir"] + "/" + psf, base_dir + molecules[mol]["dir"] + "/" + cor, topology_format="PSF", format="CRD")
        get_dftb_details(u, scan_dir)
        
        os.popen('cp {} {}'.format(template, scan_file)).readlines()
        os.popen('cp {} {}/toppar.str'.format(toppar, scan_dir)).readlines()
        
        os.popen('sed -i "s/{}/{}/g" {}'.format("MOLECULE", mol, scan_file)).readlines()
        os.popen('sed -i "s/{}/{}/g" {}'.format("PROTONATION", case, scan_file)).readlines()
        os.popen('sed -i "s/{}/{}/g" {}'.format("WHICHCORFILE", cor, scan_file)).readlines()
        os.popen('sed -i "s/{}/{}/g" {}'.format("WHICHPSFFILE", psf, scan_file)).readlines()
        os.popen('sed -i "s/{}/{}/g" {}'.format("QMCHARGE", charge, scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("INPUTDIR", base_dir + molecules[mol]["dir"], scan_file)).readlines()
        
        alc_O3_constraint = ""
        ribose_constraint = ""
        handle_constraint = ""
        nucleo_constraint = ""
        
        
        if case == "neutral":
            alc_O3_constraint += r'define O3 sele type O* .and. .bonded. atom @C3 end \n'
            alc_O3_constraint += r'set O3 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            alc_O3_constraint += r'define HO3 sele type H* .and. .bonded. atom @O3 end \n'
            alc_O3_constraint += r'set HO3 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            alc_O3_constraint += r'define HC3 sele type H* .and. .bonded. atom @C3 end \n'
            alc_O3_constraint += r'set HC3 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            alc_O3_constraint += r'cons dihe @HC3 @C3 @O3 @HO3 force 50.0 min 60.0 width 15.0 period 0 \n'
                
        
        if mol in ["RIBE", "RIBET", "RTHYNS", "RTHYMP"]:
            ribose_constraint += r'define O2 sele type O* .and. .bonded. atom @C2 end \n'
            ribose_constraint += r'set O2 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            ribose_constraint += r'define OH2 sele type H* .and. .bonded. O2 end \n'
            ribose_constraint += r'set OH2 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            ribose_constraint += r'cons dihe @C3 @C2 @O2 @OH2 force 50.0 min 180.0 width 20.0 period 0 \n'
                
            
        if mol in ["RTHYNS", "DTHYNS"]:
            handle_constraint += r'define C5 sele (type C* .and. .bonded. type O*) .and. .bonded. atom @C4 .and. .not. atom @C3 show end \n'
            handle_constraint += r'set C5 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'define O5 sele type O* .and. .bonded. C5 end \n'
            handle_constraint += r'set O5 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'define OH5 sele type H* .and. .bonded. O5 end \n'
            handle_constraint += r'set OH5 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'cons dihe @C3 @C4 @C5 @O5 force 50.0 min 180.0 width 20.0 period 0 \n'
            handle_constraint += r'cons dihe @C4 @C5 @O5 @OH5 force 50.0 min 180.0 width 20.0 period 0 \n'

        
        if mol in ["RTHYMP", "DTHYMP"]:
            handle_constraint += r'define C5 sele (type C* .and. .bonded. type O*) .and. .bonded. atom @C4 .and. .not. atom @C3 show end \n'
            handle_constraint += r'set C5 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'define O5 sele type O* .and. .bonded. C5 end \n'
            handle_constraint += r'set O5 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'define P5 sele type P* .and. .bonded. O5 end \n'
            handle_constraint += r'set P5 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'cons dihe @C3 @C4 @C5 @O5 force 50.0 min 180.0 width 20.0 period 0 \n'
            handle_constraint += r'cons dihe @C4 @C5 @O5 @P5 force 50.0 min 180.0 width 20.0 period 0 \n'

            handle_constraint += r'define O6 sele type O* .and. .bonded. atom @P5 .and. .not. .bonded. atom @C5 show end \n'
            handle_constraint += r'set O6 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'define C6 sele type C* .and. .bonded. atom @O6 end \n'
            handle_constraint += r'set C6 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'cons dihe @C5 @O5 @P5 @O6 force 50.0 min 180.0 width 20.0 period 0 \n'
            handle_constraint += r'cons dihe @O5 @P5 @O6 @C6 force 50.0 min 60.0 width 20.0 period 0 \n'

            handle_constraint += r'define O7 sele type O* .and. .bonded. atom @P5 .and. .bonded. type H* end \n'
            handle_constraint += r'set O7 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'define HO7 sele type H* .and. .bonded. atom @O7 end \n'
            handle_constraint += r'set HO7 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            handle_constraint += r'cons dihe @O5 @P5 @O7 @HO7 force 100.0 min 180.0 width 10.0 period 0 \n'
            
            handle_constraint += r'NOE \n'
            handle_constraint += r'assign sele type O* .and. .bonded. atom @P5 end - \n'
            handle_constraint += r'sele type H* .and. .bonded. type N* end - \n'
            handle_constraint += r'MINDIST kmin 200.0 rmin 3.5 rmax 100.0 KMAX 0.0 FMAX 1000.0 \n'
            handle_constraint += r'END \n'
            
            

        if mol in ["RIBET", "THFT", "RTHYNS", "DTHYNS", "RTHYMP", "DTHYMP"]:
            nucleo_constraint += r'define HC1 sele type H* .and. .bonded. atom @C1 end \n'
            nucleo_constraint += r'set HC1 = ?SELSEGI ?SELRESI ?SELTYPE \n'
            nucleo_constraint += r'define NBASE sele type N* .and. .bonded. atom @C1 end \n'
            nucleo_constraint += r'set NBASE = ?SELSEGI ?SELRESI ?SELTYPE \n'
            nucleo_constraint += r'define CBASE sele type C* .and. .bonded. atom @NBASE .and. .not. .bonded. type O* end \n'
            nucleo_constraint += r'set CBASE = ?SELSEGI ?SELRESI ?SELTYPE \n'
            nucleo_constraint += r'cons dihe @HC1 @C1 @NBASE @CBASE force 50.0 min 270.0 width 60.0 period 0 \n'
        
        
        os.popen('sed -i "s|{}|{}|g" {}'.format("ALC_O3CONS", alc_O3_constraint, scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("RIBOSECONS", ribose_constraint, scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("HANDLECONS", handle_constraint, scan_file)).readlines()
        os.popen('sed -i "s|{}|{}|g" {}'.format("NUCLEOCONS", nucleo_constraint, scan_file)).readlines()


        for i in range(len(pucker_atom_template)):
            old = pucker_atom_template[i]
            if mol in correctly_labeled:
                new = pucker_atoms_correct[i]
            elif mol in reverse_labeled:
                new = pucker_atoms_reverse[i]
            elif mol in star_labeled:
                new = pucker_atoms_star[i]
            else:
                exit()
            
            replace = "{} 1 {}".format(mol, new)
            os.popen('sed -i "s/{}/{}/g" {}'.format(old, replace, scan_file)).readlines()
            
        os.chdir(scan_dir)
        os.popen('bash {} {}'.format(submit_script, scan_file_name)).readlines()



    
