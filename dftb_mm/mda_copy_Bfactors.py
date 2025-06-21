#==============================================================================
#                                   IMPORTS
#==============================================================================
import os
import shutil
#----------------------------------
import numpy as np
#----------------------------------
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis import transformations as trans
import MDAnalysis.analysis.pca as pca
from MDAnalysis.analysis import align, nuclinfo, distances
from MDAnalysis.lib.util import convert_aa_code


#==============================================================================
#                                  FUNCTIONS
#==============================================================================
def copy_bfactors(pdb, prot, mg, state):
    global base_dir
    
    this = "{}_{}_{}mg_{}_unsolvated".format(pdb, prot, mg, state)
    psf_file = base_dir + "psf/" + this + ".psf"
    cor_file = base_dir + "cor/" + this + ".cor"
    pdb_file = base_dir + "pdb/" + pdb + ".pdb"
    
    u = mda.Universe(psf_file, cor_file, topology_format="PSF", format="CRD")
    ref = mda.Universe(pdb_file)
    
    # if you select altloc A then it ignores things which only have one position
    # selecting (not altloc B) gets all altloc A and things with only one position
    ref_select = ref.select_atoms("not (altloc B or altloc C or (resname MN or resname CA or resname DPO or resname DTP or resname HOH or resname EDO or resname GOL))")
    u_select = u.select_atoms("segid PROA or segid DNAA or segid DNAB")
    cterm_resid = [432]
    
    resnames = ref_select.atoms.resnames
    resids = ref_select.atoms.resids
    names = ref_select.atoms.names
    
    u_index = []
    ref_index = []
    for i in range(len(resnames)):
        try:
            # need to be aware of a few naming differences
            # HIS protonation state
            if resnames[i] == "HIS":
                resname_u = "HSD"
            else:
                resname_u = resnames[i]
            
            #ILE issue
            if resnames[i] == "ILE" and names[i] == "CD1":
                u_name = "CD"
                
            elif prot == "s113a" and resnames[i] == "SER" and resids[i] == 113:
                resname_u = "ALA"
                if names[i] == "OG":
                    continue
                else:
                    u_name = names[i]
            
            elif resnames[i] in ["DT", "DA", "DG", "DC"]:
                index = ["DT", "DA", "DG", "DC"].index(resnames[i])
                resname_u = ["THY", "ADE", "GUA", "CYT"][index]
                
                if names[i] == "OP1":
                    u_name = "O1P"
                elif names[i] == "OP2":
                    u_name = "O2P"
                elif names[i] == "C7":
                    u_name = "C5M"
                else:
                    u_name = names[i]
                
            # C_term carbocylic acid issue
            elif resids[i] in cterm_resid:
                if names[i] == "OXT":
                    u_name = "OT1"
                elif names[i] == "O":
                    u_name = "OT2"
                else:
                    u_name = names[i]
            #Normal naming, matches between the 2
            else:
                u_name = names[i]
                
            select_u = u_select.select_atoms("resname {} and resid {} and name {}".format(resname_u, resids[i], u_name) )
            select_ref = ref_select.select_atoms("resname {} and resid {} and name {}".format(resnames[i], resids[i], names[i]) )
            
            """
            if prot == "s113a" and resnames[i] == "SER" and resids[i] == 113:
                #print(select_u, select_ref)
                #print("u  :", select_u.n_atoms, resname_u, resids[i], u_name)
                #print("ref:", select_ref.n_atoms, resnames[i], resids[i], names[i])
            """
                
            if select_u.n_atoms == 0 or select_ref.n_atoms == 0:
               print("u  :", select_u.n_atoms, resname_u, resids[i], u_name)
               print("ref:", select_ref.n_atoms, resnames[i], resids[i], names[i])
               print("-----------------")
            # Just to highlight mismatches if any occurs
            elif select_u.n_atoms != select_ref.n_atoms:
                print(select_u, select_ref)
                print(select_u.n_atoms, select_ref.n_atoms)
                print("----------")
                """
                test0 =  u_select.select_atoms("resname {} and resid {} and not type H*".format(resname_u, resids[i]) )
                test1 =  ref_select.select_atoms("resname {} and resid {}".format(resnames[i], resids[i]) )
                for k in range(len(test0)):
                    print(test0[k], test1[k])
                    """
            else:
                u_index.append(select_u.ix)
                ref_index.append(select_ref.ix)
    
        except Exception as e:
            raise e
    
    ref_index = np.concatenate( ref_index, axis=0 )        
    u_index = np.concatenate( u_index, axis=0 )  
    
    u_bfactors = np.zeros(u.atoms.n_atoms)
    
    for i in range(len(u_index)):
        u_ix = u_index[i]
        ref_ix = ref_index[i]
        #print(ref.select_atoms("index {}".format(ref_ix)), u.select_atoms("index {}".format(u_ix)))
        u_bfactors[u_ix] = ref.select_atoms("index {}".format(ref_ix)).atoms.tempfactors
    
    u.add_TopologyAttr('tempfactors')
    for ts in u.trajectory:
        u.atoms.tempfactors = u_bfactors
        
    # JUST CHECKS IT WORKED
    #the top part of this loop is copied from above. 
    #Still have to deal with naming mismatches
    for i in range(len(resnames)):
        try:
            # need to be aware of a few naming differences
            # HIS protonation state
            if resnames[i] == "HIS":
                resname_u = "HSD"
            else:
                resname_u = resnames[i]
            
            #ILE issue
            if resnames[i] == "ILE" and names[i] == "CD1":
                u_name = "CD"
                
            elif prot == "s113a" and resnames[i] == "SER" and resids[i] == 113:
                resname_u = "ALA"
                u_name = names[i]
            
            elif resnames[i] in ["DT", "DA", "DG", "DC"]:
                index = ["DT", "DA", "DG", "DC"].index(resnames[i])
                resname_u = ["THY", "ADE", "GUA", "CYT"][index]
                
                if names[i] == "OP1":
                    u_name = "O1P"
                elif names[i] == "OP2":
                    u_name = "O2P"
                elif names[i] == "C7":
                    u_name = "C5M"
                else:
                    u_name = names[i]
                
            # C_term carbocylic acid issue
            elif resids[i] in cterm_resid:
                if names[i] == "OXT":
                    u_name = "OT1"
                elif names[i] == "O":
                    u_name = "OT2"
                else:
                    u_name = names[i]
            #Normal naming, matches between the 2
            else:
                u_name = names[i]
                
            select_u = u_select.select_atoms("resname {} and resid {} and name {}".format(resname_u, resids[i], u_name) )
            select_ref = ref_select.select_atoms("resname {} and resid {} and name {}".format(resnames[i], resids[i], names[i]) )
            
            u_b = select_u.atoms.tempfactors
            ref_b = select_ref.atoms.tempfactors
            
            # print out any instances where B fac was copied incorrectly
            if u_b.size != ref_b.size:
                print("DID NOT WORK SIZE ERROR")
                print(select_u, select_ref)
                print(u_b, ref_b)
            else:
                if  u_b != ref_b:
                    print("DID NOT WORK MATCH ERROR")
                    print(select_u, select_ref)
                    print(u_b, ref_b)
                    exit()
    
        except Exception as e:
            print(e)
            exit()
    
    # MDA will NOT use the .cor extension as much as I want it (VMD can't read
    # CHARMM .crd files) so thats why I copy and delte things to rename the output
    crd = base_dir + "cor/" + this + "_with_bfactors.crd"
    cor = base_dir + "cor/" + this + "_with_bfactors.cor"
    with mda.coordinates.CRD.CRDWriter(crd, extended=True) as CRDWRITER:
        for ts in u.trajectory:
            u.atoms.tempfactors = u_bfactors
            CRDWRITER.write( u.select_atoms('all') )
    shutil.move(crd, cor)


#==============================================================================
#                                  UNIVERSE
#==============================================================================
base_dir = "/projectnb/cui-buchem/rosadche/dnap/qm_mm/"

pdbs = ["5kfg", "5kfh"]
prots = ["wt", "s113a"]
mgs = [2, 3]

for pdb in pdbs:
    if pdb == "5kfg": state="rs"
    if pdb == "5kfh": state="ps"
    
    for prot in prots:
        for mg in mgs:
            #print(pdb, prot, mg, state)
            copy_bfactors(pdb, prot, mg, state)


