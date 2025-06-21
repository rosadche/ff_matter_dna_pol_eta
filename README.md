# ff_matter_dna_pol_eta
Repository to accompany the publication "Force Fields Matter in DNA Pol $\eta$ Mechanistic Analysis" Osadchey &amp; Cui (2025). Full citation to be uploaded upon publication. 

## Contents
**classical_umbrella_sampling**: Files for both reactant state (RS) and product state (PS)
- Restart file used at the start of umbrella sampling (from the end of classical simulation)
- Example AMBER input file
- Collective variable (CV) definition files
- Reweighting script

**classical_unbiased_md**: Starting topologies and coordinates for all classical simulations for both RS and PS, sorted by forcefield
- CHARMM36 (original setup from CHARMM-GUI)
- AMBER GAFF2, 12-6 ions (original setup from CHARMM-GUI)
- AMBER GAFF2, 12-6-4 ions (original setup from CHARMM-GUI)
- AMBER OL15, 12-6-4 ions (contains files to edit/setup from AMBER GAFF2, 12-6-4 ions topologies)
- AMBER OL15, m12-6-4 ions (contains files to edit/setup from AMBER OL15, 12-6-4 ions topologies)

**dftb_mm**: Scripts needed to setup and run DFTB/MM simulations with Generalized Solvent Boundary Potential (GSBP) from classical CHARMM36 topology and coordinate files

**qm_sugar_puckering**: Files needed to run DFTB/MM scans of model sugar systems and higher-level single points

## Notes:
- This repository does not contain scripts obtained directly from third parties (such as CHARMM-GUI) and does not provide executables or source code for programs (e.g. CHARMM, FermiONs++, Gaussian)

## Citation:
- To be updated
