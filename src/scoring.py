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
import subprocess
import warnings

# BioPython
from Bio.PDB import *
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

# Class and functions
class score_protein_protein:
    #################################################################### 
    def __init__(self,pdbID,path,src_route,chain_target,chain_binder):
        """
        Initializer
        
        Arguments:
        pdbID -- code of the system that will be used to calculate the score
        path -- path where the structure files are located
        chain_target -- contain the chain or chains the target is composed
        chain_binder -- chain of the peptide that will be modified
        
        Output:
        path_scores -- local path where the score codes are stored
        """
        self.pdbID=pdbID
        self.path=path
        self.path_scores="{}/src/scores".format(src_route)
        self.chain_target=chain_target
        self.chain_binder=chain_binder
    
    ####################################################################
    def computePisa(self):
        """
        Function to calculate Pisa score alone based on calls to the system programs.
        
        Output:
        pisa_score -- Score predicted by Pisa
        """
        path_target=self.path+"/"+self.pdbID
        target=",".join(self.chain_target)
        bash = "{}/pisaEnergy_linux {}.pdb {} {} {}/pisa.params".format(self.path_scores,path_target,target,self.chain_binder,self.path_scores)
        output = subprocess.check_output(['bash','-c', bash])
        self.pisa_score="{:.3f}".format(float(output))
        
    ####################################################################
    def computeBach(self):
        """
        Function to calculate BACH score based on calls to the system programs.
        
        Output:
        bach_score -- Score predicted by BACH
        """
        path_target=self.path+"/"+self.pdbID
        os.system("echo {}.pdb > {}/list.txt".format(path_target,self.path))
        os.system("{}/bach -PDBLIST {}/list.txt -COMPUTE_ENE -FILE_PAR {}/BSS.par -FILE_PAR_AT {}/ATOMIC_PARAMETERS_BSS -STRICT_INTERFACE -o {}/output.bss".format(self.path_scores,self.path,self.path_scores,self.path_scores,self.path))
        bash="head -2 {}/output.bss | tail -1 | awk '{{print $4}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/list.txt {}/output.bss".format(self.path,self.path))
        self.bach_score="{:.3f}".format(float(output))
    
    ####################################################################
    def computeFiredock(self):
        """
        Function to calculate Firedock score based on calls to the system programs.
        
        Output:
        firedock_score -- Score predicted by Firedock
        """
        path_target=self.path+"/"+self.pdbID
        os.system("python3 {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))
            
        os.system("cp {}/FireDock_params.txt {}/FireDock_params.temp".format(self.path_scores,self.path))
        os.system("cp -r {}/lib_firedock {}/lib".format(self.path_scores,self.path))
        os.system("sed -i 's#chainA\.pdb#{}_target\.pdb#g' {}/FireDock_params.temp".format(path_target,self.path))
        os.system("sed -i 's#chainB\.pdb#{}_{}\.pdb#g' {}/FireDock_params.temp".format(path_target,self.chain_binder,self.path))
        os.system("sed -i 's#firedock\.out#{}/firedock\.out#g' {}/FireDock_params.temp".format(self.path,self.path))
        os.system("sed -i 's#\./lib#{}/lib#g' {}/FireDock_params.temp".format(self.path,self.path))
        os.system("{}/FireDock {}/FireDock_params.temp".format(self.path_scores,self.path))
        bash="tail -n1 {}/firedock.out.unref | awk '{{print $11}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm -r {}/FireDock_params.temp {}/lib {}/firedock.out.* log {}_*.pdb {}/chains.seq".format(self.path,self.path,self.path,path_target,self.path))
        self.firedock_score="{:.3f}".format(float(output))
        
    
    ####################################################################
    def computeGoap(self):
        """
        Function to calculate GOAP scores based on calls to the system programs.
        
        Output:
        totalgoap_score -- Score predicted by DFIRE and GOAP together
        dfire_score -- Score predicted by DFIRE
        goap_score -- Score predicted by GOAP alone
        """
        path_target=self.path+"/"+self.pdbID
        bash="locate -b goap-alone | head -n1"
        out=subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
 
        os.system("cp {}/tmp.inp MD.inp".float(self.path_scores))
        os.system("cp {}.pdb .".format(path_target))
        os.system("sed -i 's#route#{}#g' MD.inp".format(out.strip()))
        os.system("sed -i 's#protein#{}\.pdb#g' MD.inp".format(self.pdbID))
        bash="{}/goap < MD.inp".format(self.path_scores)
        output=subprocess.check_output(['bash','-c', bash])
        goapData=output.split()
        os.system("rm MD.inp {}.pdb".format(self.pdbID))
        self.totalgoap_score="{:.3f}".format(float(goapData[2]))
        self.dfire_score="{:.3f}".format(float(goapData[3]))
        self.goap_score="{:.3f}".format(float(goapData[4]))
    
    #################################################################### 
    def computeZrank(self):
        """
        Function to calculate ZRANK score based on calls to the system programs.
        
        Output:
        zrank_score -- Score predicted by ZRANK
        """
        path_target=self.path+"/"+self.pdbID
        
        # Fix hydrogens
        parser = PDBParser()
        structure = parser.get_structure('PEP',"{}.pdb".format(path_target))
        io = PDBIO()
        io.set_structure(structure)
        io.save('{}_fixed.pdb'.format(path_target))
    
        # Run zrank
        os.system("echo {}_fixed.pdb > {}/list.txt".format(path_target,self.path))
        os.system("{}/zrank {}/list.txt".format(self.path_scores,self.path))
        bash="awk '{{print $2}}' {}/list.txt.zr.out".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/list.txt {}/list.txt.zr.out {}_fixed.pdb".format(self.path,self.path,path_target))
        self.zrank_score="{:.3f}".format(float(output))
        
    ####################################################################
    def computeIrad(self):
        """
        Function to calculate IRAD score based on calls to the system programs.
        
        Output:
        irad_score -- Score predicted by IRAD
        """        
        path_target=self.path+"/"+self.pdbID
        
        # Fix hydrogens
        parser = PDBParser()
        structure = parser.get_structure('PEP',"{}.pdb".format(path_target))
        io = PDBIO()
        io.set_structure(structure)
        io.save('{}_fixed.pdb'.format(path_target))
    
        # Run irad
        os.system("echo {}_fixed.pdb > {}/list.txt".format(path_target,self.path))
        os.system("{}/irad {}/list.txt".format(self.path_scores,self.path))
        bash="awk '{{print $2}}' {}/list.txt.irad.out".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/list.txt {}/list.txt.irad.out {}_fixed.pdb".format(self.path,self.path,path_target))
        self.irad_score="{:.3f}".format(float(output))
    
    #################################################################### 
    def computeBMF(self):
        """
        Function to calculate BMF-BLUUES scores based on calls to the system programs.
        
        Output:
        totalbmfbluues_score -- Score predicted by BMF and BLUUES together
        bmf_score -- Score predicted by BMF
        bluues_score -- Score predicted by BLUUES
        """
        
        path_target=self.pdbID
        pdb=self.path+"/"+self.pdbID
        os.system("cp {}.pdb .".format(pdb))
        os.system("python3 {}/get_chains.py {}.pdb .".format(self.path_scores,path_target))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))
            
        os.system("bash {}/compute_BMF_BLUEES.sh {} {}".format(self.path_scores,path_target,self.path_scores))
        os.system("bash {}/compute_BMF_BLUEES.sh {}_target {}".format(self.path_scores,path_target,self.path_scores))
        os.system("bash {}/compute_BMF_BLUEES.sh {}_{} {}".format(self.path_scores,path_target,self.chain_binder,self.path_scores))
        os.system("bash {}/final_BMF_BLUUES.sh {} {}_target {}_{} {}".format(self.path_scores,path_target,path_target,path_target,self.chain_binder,self.path_scores))
        bash="awk '{{print $1}}' {}.score_final".format(path_target)
        outputTotal = subprocess.check_output(['bash','-c', bash])
        bash="awk '{{print $1}}' {}.score_bmf_final".format(path_target)
        outputBMF = subprocess.check_output(['bash','-c', bash])
        bash="awk '{{print $1}}' {}.score_bluues_final".format(path_target)
        outputBLUUES = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}_*.pdb {}.score_final {}.score_bmf_final {}.score_bluues_final chains.seq {}.pdb".format(path_target,path_target,path_target,path_target,self.pdbID))
        
        self.totalbmfbluues_score="{:.3f}".format(float(outputTotal.strip().decode("utf-8")))
        self.bmf_score="{:.3f}".format(float(outputBMF.strip().decode("utf-8")))
        self.bluues_score="{:.3f}".format(float(outputBLUUES.strip().decode("utf-8")))
    
    ####################################################################
    def computeVina(self):
        """
        Function to calculate Vina score based on calls to the system programs.
        
        Output:
        vina_score -- Score predicted by Vina
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python3 {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))
                
        #os.system("./{}/prepare_receptor4.py -r {}_target.pdb -o {}_target.pdbqt".format(self.path_scores,path_target,path_target))
        #os.system("./{}/prepare_ligand4.py -l {}_{}.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_{}.pdbqt".format(self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        #os.system("./{}/vina --receptor {}_target.pdbqt --ligand {}_{}.pdbqt --score_only --log {}/score.log".format(self.path_scores,path_target,path_target,self.chain_binder,self.path))
        
        os.system("./{}/prepare_receptor4.py -r {}_{}.pdb -o {}_{}.pdbqt".format(self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("./{}/prepare_ligand4.py -l {}_target.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_target.pdbqt".format(self.path_scores,path_target,path_target))
        os.system("./{}/vina --receptor {}_{}.pdbqt --ligand {}_target.pdbqt --score_only --log {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path))
        
        # Check the number of CPUs to filter the VINA file result
        bash="lscpu | grep -E '^CPU\(' | awk '{print $2}'"
        cores = subprocess.check_output(['bash','-c', bash])
        if int(cores) <= 8: 
            bash="head -n18 {}/score.log | tail -n1 | awk '{{print $2}}'".format(self.path)
        else:
            bash="head -n19 {}/score.log | tail -n1 | awk '{{print $2}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.pdbqt {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.vina_score="%0.3f" %float(output)
    
    ####################################################################
    def computeSmina(self):
        """
        Function to calculate Smina score based on calls to the system programs.
        
        Output:
        smina_score -- Score predicted by Smina
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python3 {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))
                
        # os.system("./{}/prepare_receptor4.py -r {}_target.pdb -o {}_target.pdbqt".format(self.path_scores,path_target,path_target))
        # os.system("./{}/prepare_ligand4.py -l {}_{}.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_{}.pdbqt".format(self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        # os.system("./{}/smina.static --receptor {}_target.pdbqt --ligand {}_{}.pdbqt --score_only --log {}/score.log".format(self.path_scores,path_target,path_target,self.chain_binder,self.path))
        os.system("./{}/prepare_receptor4.py -r {}_{}.pdb -o {}_{}.pdbqt".format(self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("./{}/prepare_ligand4.py -l {}_target.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_target.pdbqt".format(self.path_scores,path_target,path_target))
        os.system("./{}/smina.static --receptor {}_{}.pdbqt --ligand {}_target.pdbqt --score_only --log {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path))
        
        # Filter the SMINA file result
        bash="head -n22 {}/score.log | tail -n1 | awk '{{print $2}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.pdbqt {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.smina_score="%0.3f" %float(output)
    
    ####################################################################
    def computeCyscore(self):
        """
        Function to calculate Cyscore based on calls to the system programs.
        
        Output:
        cyscore_score -- Score predicted by Cyscore
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python3 {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))
                
        #os.system("babel -ipdb {}_{}.pdb -omol2 {}_{}.mol2".format(path_target,self.chain_binder,path_target,self.chain_binder))
        #os.system("./{}/cyscore_bin/Cyscore {}_target.pdb {}_{}.mol2 > {}/score.log".format(self.path_scores,path_target,path_target,self.chain_binder,self.path))
        os.system("babel -ipdb {}_target.pdb -omol2 {}_target.mol2".format(path_target,path_target))
        os.system("./{}/cyscore_bin/Cyscore {}_{}.pdb {}_target.mol2 > {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path))
        
        # Filter the Cyscore file result
        bash="grep 'Cyscore=' {}/score.log | awk '{{print $2}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.mol2 {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.cyscore_score="%0.3f" %float(output)
    
    ####################################################################
    def computeNNscore(self):
        """
        Function to calculate NNscore2.0 score based on calls to the system programs.
        
        Output:
        nnscore_score -- Score predicted by NNscore
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python3 {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))
                
        #os.system("./{}/prepare_receptor4.py -r {}_target.pdb -o {}_target.pdbqt".format(self.path_scores,path_target,path_target))
        #os.system("./{}/prepare_ligand4.py -l {}_{}.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_{}.pdbqt".format(self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("./{}/prepare_receptor4.py -r {}_{}.pdb -o {}_{}.pdbqt".format(self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("./{}/prepare_ligand4.py -l {}_target.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_target.pdbqt".format(self.path_scores,path_target,path_target))
        os.system("python3.5 {}/NNScore2.py -receptor {}_{}.pdbqt -ligand {}_target.pdbqt -vina_executable {}/vina > {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path_scores,self.path))
        
        # Filter the NNscore file result
        bash="grep 'Best Score:' {}/score.log | awk '{{print $4}}' | cut -f 1 -d ','".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        output_neg=-1.0*float(output)
        os.system("rm {}/score.log {}/*.pdbqt {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.nnscore_score="%0.3f" %output_neg
    
    ####################################################################
    def computeDligand(self):
        """
        Function to calculate Dligand score based on calls to the system programs.
        
        Output:
        dligand_score -- Score predicted by Dligand
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python3 {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))
                
        #os.system("babel -ipdb {}_{}.pdb -omol2 {}_{}.mol2".format(path_target,self.chain_binder,path_target,self.chain_binder))
        #os.system("./{}/cyscore_bin/Cyscore {}_target.pdb {}_{}.mol2 > {}/score.log".format(self.path_scores,path_target,path_target,self.chain_binder,self.path))
        os.system("babel -ipdb {}_target.pdb -omol2 {}_target.mol2".format(path_target,path_target))
        os.system("cp {}/dfire.2 .".format(self.path_scores))
        os.system("cp {}/amino.mol2 .".format(self.path_scores))
        os.system("./{}/dligand2.intel -P {}_{}.pdb -L {}_target.mol2 > {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path))
        
        # Filter the dligand file result
        bash="head -n1 {}/score.log | tail -n1".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.mol2 {}_*.pdb {}/chains.seq dfire.2 amino.mol2".format(self.path,self.path,path_target,self.path))
        self.dligand_score="%0.3f" %float(output)
    
    ####################################################################
    def computeDSXscore(self):
        """
        Function to calculate DSXcore based on calls to the system programs.
        
        Output:
        dsxscore_score -- Score predicted by DSXscore
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python3 {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))
                
        #os.system("babel -ipdb {}_{}.pdb -omol2 {}_{}.mol2".format(path_target,self.chain_binder,path_target,self.chain_binder))
        #os.system("./{}/cyscore_bin/Cyscore {}_target.pdb {}_{}.mol2 > {}/score.log".format(self.path_scores,path_target,path_target,self.chain_binder,self.path))
        os.system("babel -ipdb {}_target.pdb -omol2 {}_target.mol2".format(path_target,path_target))
        os.system("./{}/dsx_linux_64.lnx -P {}_{}.pdb -L {}_target.mol2 -D {}/pdb_pot_0511 -F {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path_scores,self.path))
        
        # Filter the DSXscore file result
        bash="sed '/^$/d' {}/score.log | tail -n1 | awk '{{print $7}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.mol2 {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.dsxscore_score="%0.3f" %float(output)
