class chain_not_found(Exception):
    def __init__(self, chains):
        self.chains = chains
    def __str__(self):
        if len(self.chains) == 1:
            return("The chain %s, which has been retrieved from the PDBs, is "
            "not found in the fasta. It seems that the sequence of this "
            "chain does not align with any of the sequences of the "
            "subunits that have been found in the fasta. Check if the "
            "fasta is the correct one.\nIf it is, you can lower the "
            "threshold of this alignment manually or you can obtain the "
            "fasta directly from the PDBs using the PrediProt_getfasta "
            "script." %(self.chains))
        else:
            return("The chains%s, which have been retrieved from the PDBs, are "
            "not found in the fasta. It seems that the sequences of "
            "these chains do not align with any of the sequences of the "
            "subunits that have been found in the fasta. Check if the "
            "fasta is the correct one.\nIf it is, you can lower the "
            "threshold of this alignment manually or you can obtain the "
            "fasta directly from the PDBs using the PrediProt_getfasta "
            "script." %(self.chains))


class output_path_not_exists(Exception):
    def __init__(self,path):
        self.path=path
    def __str__(self):
        return "There is no path %s to create the output directory" %(self.path)


class fasta_not_exists(Exception):
    def __init__(self, fasta):
        self.fasta = fasta
    def __str__(self):
        return "There is no file called %s" %(self.fasta)


class directory_not_exists(Exception):
    def __init__(self, dir):
        self.dir = dir
    def __str__(self):
        return "There is no directory called %s" %(self.dir)


class invalid_clash_dist(Exception):
    def __init__(self, dist):
        self.dist = dist
    def __str__(self):
        return "The value %s is not valid, it must be from 0.1 to 1.0" %(self.dist)


class two_chains_each(Exception):
    def __init__(self, pdb):
        self.pdb = pdb
    def __str__(self):
        return "The file %s does not have two chains" %(self.pdb)


class no_pdbs(Exception):
    def __str__(self):
        return "There are not pdbs in the selected directory"


class pdb_error(Exception):
    def __init__(self, pdb):
        self.pdb = pdb
    def __str__(self):
        return("The PDB file %s can not be read. Please check if the format is "
        "correct") %(self.pdb)


class not_interactions_found(Exception):
    def __str__(self):
        return "No interactions found. Check the directory of the PDBs"


class not_fasta_file(Exception):
    def __str__(self):
        return "Please select a file with a correct fasta extension"


class not_chains_added(Exception):
    def __str__(self):
        return("No chains have been added to the model, please check the "
        "selected stoichiometry ")
