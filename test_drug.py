#!/usr/bin/python3

"""
PARCE: Protocol for Amino acid Refinement through Computational Evolution

From publication "PARCE: Protocol for Amino acid Refinement through Computational Evolution"
Computer Physics Communications 
Authors: Rodrigo Ochoa, Miguel A. Soler, Alessandro Laio, Pilar Cossio
Year: 2020

Third-party tools required:

- Scwrl4: http://dunbrack.fccc.edu/scwrl4/license/index.html
- Gromacs 5.1.4 (tested version): http://manual.gromacs.org/documentation/5.1.4/download.html
- BioPython: https://biopython.org/wiki/Download - Ubuntu package: python3-rdkit
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "Miguel A. Soler", "Alessandro Laio", "Pilar Cossio"]
__license__ = "MIT"
__version__ = "1.5"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################
import os
from src import scoring
from src import mutation
from src import general
# BioPython
from Bio.PDB import *

########################################################################################
# Test
########################################################################################

if __name__ == '__main__':
    
    bash = "pwd | cut -f 1"
    route = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
    
    # List of variables that should be defined to create the class
    chain="B"
    pdbID="drug_CQRTRFFKWYRC"
    score_list=["vina","smina","cyscore","nnscore","dligand","dsxscore"]
    threshold=3
    half_flag=0
    sim_time="5"
    num_mutations=5
    try_mutations=10
    target="drug"
    t_effective=0.5
    scoring_s="consensus"
    src_route=route
    # Additional variables for the main script
    residues_mod=[2,3,4,5,6,7,8,9,10,11]
    
    # Variables sensible to the mode: start or restart
    folder="test_drug"
    mode="start"
    peptide_reference="CQRTRFFKWYRC"
    md_route=route+"/design_input/peptide_drug"
    md_original="drug_example"
    
    iteration=0
    score_dictionary={}
    score_mode="all"
    
    # Variables for mutation test
    scwrl_path="/usr/local/bin/scwrl4/Scwrl4"
    gmxrc_path="/usr/local/gromacs/bin/GMXRC"
    os.system(". {}".format(gmxrc_path))
    peptide_mutated="CQRTRSFKWYRC"
    mutation_method="faspr"
    # File containing the results
    report=open("report_test_drug.txt","w")
    flag_step1=0
    flag_step2=0
    flag_step3=0
    
    # 1. Test Gromacs functionality
    try:
        print("1. Starting the test of Gromacs functionalities ...")
        protein_complex=general.complex(chain,pdbID,iteration,score_list,threshold,t_effective,num_mutations,scoring_s,target,score_mode,sim_time,mode,try_mutations)
        protein_complex.configure_folder(folder,src_route,md_route,md_original)
        protein_complex.setup(folder)
        print("####################################")
        print("The Gromacs test passed successfully")
        report.write("The Gromacs test passed successfully\n")
        print("####################################")
        flag_step1=1
    except:
        print("####################################")
        print("The Gromacs test failed. Please verify that Gromacs is installed, as well as the GromacsWrapper package in python")
        report.write("The Gromacs test failed. Please verify that Gromacs is installed, as well as the GromacsWrapper package in python\n")
        print("####################################")
    
    # 2. Test mutation functionality
    if mutation_method=="scwrl4":
        try:
            print("2. Starting the test of FASPR functionalities ...")
            parser = PDBParser()
            reference = parser.get_structure('REF',"design_output/{}/{}.pdb".format(folder,pdbID))
            local_mutation=mutation.mutate_peptide("design_output/{}".format(folder),peptide_reference,peptide_reference,6,chain,reference,"F","S",0,["A"],0,scwrl_path,target,src_route)
            local_mutation.replace_amino_acid()
            local_mutation.mutate_faspr()
            if os.path.isfile('design_output/{}/complex.pdb'.format(folder)):
                print("####################################")
                print("The FASPR test passed successfully")
                report.write("The FASPR test passed successfully\n")
                print("####################################")
                flag_step2=1
        except:
            print("####################################")
            print("The FASPR test failed. Please verify the permissions to execute FASPR are activated")
            report.write("The FASPR test failed. Please verify the permissions to execute FASPR are activated\n")
            print("####################################")
    if mutation_method=="scwrl4":
        try:
            print("2. Starting the test of Scwrl4 functionalities ...")
            parser = PDBParser()
            reference = parser.get_structure('REF',"design_output/{}/{}.pdb".format(folder,pdbID))
            local_mutation=mutation.mutate_peptide("design_output/{}".format(folder),peptide_reference,peptide_reference,6,chain,reference,"F","S",0,["A"],0,scwrl_path,target,src_route)
            local_mutation.replace_amino_acid()
            local_mutation.mutate_scwrl()
            if os.path.isfile('design_output/{}/complex.pdb'.format(folder)):
                print("####################################")
                print("The Scwrl4 test passed successfully")
                report.write("The Scwrl4 test passed successfully\n")
                print("####################################")
                flag_step2=1
        except:
            print("####################################")
            print("The Scwrl4 test failed. Please verify that Scwrl4 is correctly installed, and the permissions to execute are activated")
            report.write("The Scwrl4 test failed. Please verify that Scwrl4 is correctly installed, and the permissions to execute are activated\n")
            print("####################################")
    
    # 3. Test the scoring functions    
    #if flag_step3==0:
    try:
        print("3. Starting the test of scoring functions ...")
        # Call the scoring functions
        sc=scoring.score_protein_protein("complex","design_output/{}".format(folder),src_route,["A"],chain)
        # Calculate the designated scores
        for s in score_list:
            if s=="vina":
                sc.computeVina()    
                print(float(sc.vina_score))
                report.write("Vina score: {}\n".format(sc.vina_score))
            if s=="smina":
                sc.computeSmina()    
                print(float(sc.smina_score))
                report.write("Smina score: {}\n".format(sc.smina_score))
            if s=="cyscore":
                sc.computeCyscore()    
                print(float(sc.cyscore_score))
                report.write("Cyscore score: {}\n".format(sc.cyscore_score))
            if s=="nnscore":
                sc.computeNNscore()    
                print(float(sc.nnscore_score))
                report.write("NNscore score: {}\n".format(sc.nnscore_score))
            if s=="dligand":
                sc.computeDligand()    
                print(float(sc.dligand_score))
                report.write("Dligand score: {}\n".format(sc.dligand_score))
            if s=="dsxscore":
                sc.computeDSXscore()    
                print(float(sc.dsxscore_score))
                report.write("DSXscore score: {}\n".format(sc.dsxscore_score))
        print("####################################")
        print("The scoring functions test passed successfully")
        report.write("The scoring functions test passed successfully\n")
        print("####################################")
        flag_step3=1
    except:
        print("####################################")
        print("The scoring functions test failed. Please verify that the required scoring files are present in the src/scores folder")
        report.write("The scoring functions test failed. Please verify that the required scoring files are present in the src/scores folder\n")
        print("####################################")
    
    if flag_step1==1 and flag_step2==1 and flag_step3==1:
        print("####################################")
        print("Everything is ready. You can start the protocol :)")
        print("####################################")
        report.write("####################################\n")
        report.write("Everything is ready. You can start the protocol :)\n")
        report.write("####################################\n")
        
    # Close the report
    report.close()
