#==============================================================================
#                                   IMPORTS
#==============================================================================
import os
import pickle
#----------------------------------
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis import transformations as trans
import MDAnalysis.analysis.pca as pca
from MDAnalysis.analysis import align, nuclinfo, distances
from MDAnalysis.lib.util import convert_aa_code
#----------------------------------

def mda_selection(pdb, prot, mg, state):
    #==============================================================================
    #                                  UNIVERSE
    #==============================================================================
    base_dir = "/projectnb/cui-buchem/rosadche/dnap/qm_mm/"
    dftb_file_dir = base_dir + "setup/" + "{}_{}_{}mg_{}/".format(pdb, prot, mg, state)
    
    simname = "{}_{}_{}mg_{}".format(pdb, prot, mg, state)
    term = "gsbp_setup_pt1"
    total = simname + "_{}".format(term)
    
    psf = base_dir + "psf/" + total + ".psf"
    cor = base_dir + "cor/" + total + ".cor"
    
    u = mda.Universe(psf, cor, topology_format="PSF", format="CRD")
    
    #-------------
    # Selections
    #-------------
    segname = "PROA"
    qm_residues = ["ASP_13", "ARG_55", "ARG_61", "SorA_113", "ASP_115", "GLU_116", "K_231"]
    selec_dict = dict()
    
    for res in qm_residues:
        resnum = res.split("_")[1]
        selec_dict[res] = u.select_atoms("segid {} and resid {} and not (name N HN CA HA C O)".format(segname, resnum))
    
    # backbone of M_14 to C_16
    backbone_14 = u.select_atoms("segid PROA and resid 14 and (name QQH* C O)")
    backbone_15 = u.select_atoms("segid PROA and resid 15 and (name N HN CA HA C O QQH*)")
    backbone_16 = u.select_atoms("segid PROA and resid 16 and (name QQH* N HN CA HA)")
    selec_dict["back_14"] = backbone_14
    selec_dict["back_15"] = backbone_15
    selec_dict["back_16"] = backbone_16

    # Magnesiums
    selec_dict["mgs"] = u.select_atoms("(segid MGA MGB MG3) and type MG")
    
    # ligand, segid LIG for all cases, 
    # just depends on whether that was origanlly DTP or DPO from MM
    selec_dict["lig"] = u.select_atoms("segid LIG")
    
    #DNA QM
    # for reactant & product, need DNAB 8 but not all
    selec_dict["DNAB_8"] = u.select_atoms("segid DNAB and resid 8 and not (name C5' O5' O1P O2P P H5' H5'')")
        
    if state == "ps":
        # for product we also need DNAB 9 (all)
        selec_dict["DNAB_9"] = u.select_atoms("segid DNAB and resid 9")
    
    # dist=wat/tot: 7.5=21/188, 7.65=30/197, 7.7=33/200, 7.8=39/206
    fire_center = u.select_atoms("segid MGA") 
    selec_dict["fire_water"] = u.select_atoms("byres (resname TIP3 and type O* and around 7.65 (group NAMEGROUP))", NAMEGROUP=fire_center)
    
    i = 0
    for key in selec_dict.keys():
        if i == 0:
            select_qm_all = selec_dict[key]
        else:
            select_qm_all += selec_dict[key]
        i += 1
     
    selec_dict["qm_all"] = select_qm_all
    
    print(simname)
    print(selec_dict["fire_water"].n_atoms)
    print(selec_dict["qm_all"].n_atoms)
    
    # if this is not QM then it will suddenly become a QM atom = bad
    # however, if it is QM it wil double count in the qm_all selection files
    # but it can't contribute twice to anything so its fine
    selec_dict["fire_center"] = fire_center
    
    #-------------
    # Create selection direc and write out
    #-------------
    sele_dir = base_dir +"select/"
    if not os.path.exists(sele_dir):
        os.makedirs(sele_dir)
    
    preamble = "\n* Fake title\n*\n"
    with mda.selections.charmm.SelectionWriter(sele_dir + simname + "_charmm_selections.str", mode='w', preamble=preamble) as ndx:
        for key in selec_dict.keys():
            ndx.write(selec_dict[key], name=key)
    
    key = "qm_all"
    preamble="source LOCALPATH/" + simname + "_vmd_selections.vmd set sel [atomselect top {}]".format(key)
    with mda.selections.vmd.SelectionWriter(sele_dir + simname + "_vmd_selections.vmd", mode='w', preamble=preamble) as ndx:
            ndx.write(selec_dict[key], name=key)
    
    preamble="source LOCALPATH/@" + simname + "_pymol_selections.pml"
    with mda.selections.pymol.SelectionWriter(sele_dir + simname + "_pymol_selections.pml", mode='w', preamble=preamble) as ndx:
            ndx.write(selec_dict[key], name=key)
    
    #-------------
    # SETUP SCCDFTB.DAT FILE BASED ON QM REGION
    #-------------
    dftb_par_dir = base_dir + "dftb_par"
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
    
    unique_atom_types = sorted( list( set(selec_dict["qm_all"].elements) ) )
    print("Unique QM atom types", unique_atom_types)
    print('-------')
        
    #start by writing spl file location
    text = ""
    for atom1 in unique_atom_types:
        for atom2 in unique_atom_types:
            if "P" in [atom1, atom2]:
                directory = dft_par_phos_dir
                file = "{}/{}{}.spl".format(directory, atom1.lower(), atom2.lower())
                if os.path.exists(file) == False:
                    directory = dftb_par_dir
                    file = "{}/{}{}.spl".format(directory, atom1.lower(), atom2.lower())
                    if os.path.exists(file) == False:
                        print("Error, {}{}.spl doesn't appear to exist".format(atom2.lower(), atom1.lower()))
                        exit()
            else:
                directory = dftb_par_dir
                file = "{}/{}{}.spl".format(directory, atom1.lower(), atom2.lower())
                if os.path.exists(file) == False:
                    print("Error, {}{}.spl doesn't appear to exist".format(atom2.lower(), atom1.lower()))
                    exit()
            
            text += "\'" + file + "\'\n"
    
    #then by writing hubbard paramters
    for atom1 in unique_atom_types:
        text += "\'" + atom1.upper() + "\'  " + "{}\n".format(HUBBARD[atom1.upper()])
      
    #then write zeta for HBON correction, which you probably should be using  
    text += "{}\n".format(ZETA)
    
    with open(dftb_file_dir + "sccdftb.dat", "w") as file:
        file.write(text)
    
    i = 1
    with open(sele_dir + simname + "_sccdftb_wmain.str", "w") as file:
        file.write("* Fake title\n*\n")
        for atom1 in unique_atom_types:
            text = "scalar WMAIN set {:.1f} sele (qm_all) .and. type {}*  SHOW end\n"
            if atom1 == "H":
                file.write( text.format(i, "QQ") )
            file.write( text.format(i, atom1) )
            i += 1
    
    
    crd = base_dir + "cor/" + simname + "_qm_region.crd"
    with mda.coordinates.CRD.CRDWriter(crd, extended=True) as CRDWRITER:
        for ts in u.trajectory:
            CRDWRITER.write( selec_dict["qm_all"] )



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
            mda_selection(pdb, prot, mg, state)
            
            
            


