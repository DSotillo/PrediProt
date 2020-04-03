from Bio import pairwise2
import Bio.PDB
import copy
import random
import re
import os
import sys
from prediprot_classes import *


def get_dict_from_fasta(args):
    '''
    The fasta is read and stored in a dictionary with the IDs as keys and
    the sequence as values. It also returns a list with all the chains id.
    '''

    fasta_file = open(args.inputfasta)
    list_ID = []
    fasta_dic = {}
    cur_id = ''
    cur_seq = []
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">") and cur_id == '':
            cur_id = line.split(' ')[0][1:]
            list_ID.append(cur_id)
        elif line.startswith(">") and cur_id != '':
            fasta_dic[cur_id] = ''.join(cur_seq)
            cur_id = line.split(' ')[0][1:]
            list_ID.append(cur_id)
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    fasta_dic[cur_id] = ''.join(cur_seq)
    return fasta_dic, list_ID


def get_subunits_from_fasta(fasta_dic, list_ID):
    '''
    It obtains the different unique subunits of the complex. Two chains are
    considered the same subunit if they share more than 95% of identity. Unique
    chains are renamed with a new id starting from letter A. It is returned:
    a dictionary linking unique chains with sequence, a dictionary that links
    unique chains with the times that chain appears in the complex and a list
    with the unique chain ids.
    '''

    unique = []
    repeated = []
    dic_repeated = {}
    threshold = 0.95
    n_chain_name = 0
    fasta_dic_unique = {}
    unique_new_id = []
    for ID1 in list_ID:
        if ID1 in repeated:
            continue
        for ID2 in list_ID:
            if ID1 != ID2:
                alignment = pairwise2.align.globalxx(fasta_dic[ID1], fasta_dic[ID2],
                                               one_alignment_only = True)
                length_sequence = min(len(fasta_dic[ID1]), len(fasta_dic[ID2]))
                identity = alignment[0][2] / len(alignment[0][0])
                if ID1 not in unique:
                    unique.append(ID1)
                    id = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456"[n_chain_name:n_chain_name+1]
                    fasta_dic_unique[id] = fasta_dic[ID1]
                    n_chain_name += 1
                    dic_repeated[id] = 1
                    unique_new_id.append(id)
                if identity > threshold:
                    repeated.append(ID2)
                    dic_repeated[id] += 1

    # The subunits are sorted ir order to store the most repeated ones at the
    # beginning in the dictionaries, e.g., A:8,B:6,C:4,D:1...
    sorted_subunits = sorted(dic_repeated.items(), key=lambda x: x[1],reverse=True)
    n_chain_name = 0
    fasta_dic_unique_sorted = {}
    dic_repeated_sorted = {}
    for id in sorted_subunits:
        new_id = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789"[n_chain_name:n_chain_name+1]
        dic_repeated_sorted[new_id] = id[1]
        fasta_dic_unique_sorted[new_id] = fasta_dic_unique[id[0]]
        n_chain_name += 1
    return fasta_dic_unique_sorted, dic_repeated_sorted, unique_new_id


def create_subunits_fasta(args, fasta_dic_unique_sorted, dic_repeated_sorted,
                          unique_new_id, name_fasta):
    '''
    If the user wants info about the stoichiometry, it is created a fasta
    with the unique subunits.
    '''

    fasta_subunits = open(name_fasta + '_subunits.fasta','w')
    stoic = ''
    for id in unique_new_id:
        stoic += id + ':' + str(dic_repeated_sorted[id]) + ','
        fasta_subunits.write('>' + name_fasta + '_' + id + '\n'
                             + fasta_dic_unique_sorted[id] + '\n')
    # Also a messeage with the stoichiometry obtained from the fasta appears
    print('A file named '+ name_fasta
    + '_subunits.fasta with the unique subunits has been created' + '\n'
    + 'From the given fasta, the stoichiometry of the complex is ' + stoic[:-1]
    + '\n' + 'Please choose the subunits you want to see in the output:')
    # Then the user can insert the stoichiometry that wants to see
    args.sto = input()
    return args.sto


def check_stoic (args):
    '''
    It checks if the format of the given stoichiometry is correct. If it is,
    the selected stoichiometry is stored in a dictionary.
    '''

    while args.sto:
        tag = 0
        sto_list = args.sto.split(',')
        for res in sto_list:
            result = re.search('^(.+):([0-9]+)$',res)
            if not result:
                tag = 1
        if tag == 1:
            print("The format of the given stoichiometry is not valid, please "
            "insert it correctly. For example: A:2,B:2,C:1")
            args.sto = input()
        else:
            sto_dic = {}
            for dic in sto_list:
                sto_dic[dic.split(':')[0]] = int(dic.split(':')[1])
            return sto_dic
            break


def read_and_store_pdbs(args, pdb_list, unique_new_id, fasta_dic_unique_sorted):
    '''
    The dictionary AA is used to get every sequence in the pdbs and then they
    are aligned with the fasta of the unique subunits. If there is an alignment,
    store the seq in the chain_sequences dictionary. If a sequence does not
    align with any subunit, that means that the fasta given is not adeccuate,
    so the threshold of the alignment could be changed in order to identify all
    subunits from the pdbs. Also there is an additional script to create the
    fasta from the sequences of the initial pdb file, doing that there is no
    problem in identifying all the subunits.
    '''

    pdbparser = Bio.PDB.PDBParser()
    AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
          'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
          'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
          'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
          'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
          ' DA':'A', ' DG':'G', ' DC':'C', ' DT':'T',
          '  A':'A', '  G':'G', '  C':'C', '  T':'T',
          '  U':'U', 'UNK':'X', 'TERM':''}

    threshold=0.95
    # Dictionary that stores the pdb names as keys and the structure as value.
    id_structure_dict = {}
    # Dictionary that stores the pdb names + chain as keys and the chains as value.
    ID_subunit_dict = {}
    # Dictionary that stores the pdb names + chain as keys and the aa sequence as value.
    chain_sequences  = {}
    # Dictionary that stores the id of the chain as key and the id of the
    # corresponding subunit as a value.
    recover_id_from_fasta = {}
    # List that stores the chains that didn't aligned.
    chains_not_found = []
    directory = args.inputpdb

    for pdb_file in pdb_list:
        # Read every pdb and store it in the dictionary.
        # If it can not be read, raise an error.
        try:
            id_structure_dict[pdb_file] = pdbparser.get_structure(pdb_file, directory + "/" + pdb_file + ".pdb")[0]
        except:
            print(pdb_error(pdb_file))
            sys.exit()

        # If there are not 2 chains in a PDB raise an exception indicating where is the error.
        try:
            if len(list(id_structure_dict[pdb_file].get_chains()))!=2:
                raise two_chains_each()
        except:
            print(two_chains_each(pdb_file))
            sys.exit()

        for chain in id_structure_dict[pdb_file].get_chains():
            sequence = ""
            for residue in chain.get_residues():
                if residue.get_resname() in AA:
                    sequence += AA[residue.get_resname()]
            # Store in dictionaries
            for id in unique_new_id:
                alignment = pairwise2.align.globalxx(sequence, fasta_dic_unique_sorted[id],
                                               one_alignment_only = True)
                # alignment[0][2] = score. alignment[0][0] = alignment including gaps.
                # Simply: 1 point per match, 0 points missmatch.
                # If 100 aminoadids perfectly aligned, score = 100.
                identity = alignment[0][2] / len(alignment[0][0])
                if identity > threshold:
                    chain_sequences[pdb_file + chain.get_id()] = sequence
                    ID_subunit_dict[pdb_file + chain.get_id()] = chain
                    recover_id_from_fasta[chain.get_id()] = id

    # List the chains that are not aligned with any of the subunits of the fasta.
        if chain.get_id() not in recover_id_from_fasta:
            if chain.get_id() not in chains_not_found:
                chains_not_found.append(chain.get_id())

    chains_not_found = sorted(chains_not_found)

    try:
        if len(chains_not_found) != 0:
            chains = ''
            if len(chains_not_found) == 1:
                chains = chains_not_found[0]
            else:
                for chain in chains_not_found[:-1]:
                    chains += ' ' + chain + ','
                chains = chains[:-1]
                chains += ' and ' + chains_not_found[-1]
            raise chain_not_found()
    except:
        print(chain_not_found(chains))
        sys.exit()

    return id_structure_dict, ID_subunit_dict, chain_sequences, recover_id_from_fasta


def find_interactions(chain_sequences):
    '''
    Pairwise alignments between the chains is performed so interaction
    between structures is found and stored in a dictionary that is returned.
    '''

    # Sort the sequences to iterate through them
    sorted_sequences = sorted(chain_sequences.items())
    # Index is used in order not to repeat the pairwise alignments already done,
    # to just compare the sequences that are following sequence1.
    index = 0
    threshold = 0.95
    # List storing interactions.
    # Each interaction is: (structure 1, chain 1, position of the aminoacids
    # that are aligned of seq1, structure 2, chain 2, position of...).
    interacting_structures = []

    for chain1, sequence1 in sorted_sequences:
        # We check only the following sequences to avoid repetition.
        for chain2, sequence2 in sorted_sequences[index:]:
            # If the structures they belong to are not the same...
            if chain1[:-1] != chain2[:-1]:
                # Perform a global alignment between the two chain sequences
                alignment = pairwise2.align.globalxx(chain_sequences[chain1], chain_sequences[chain2],
                                                   one_alignment_only = True)
                identity = alignment[0][2]/len(alignment[0][0])
                # Sequences indicating matching aminoacids and gaps.
                sequence1 = alignment[0][0]
                sequence2 = alignment[0][1]
                if identity > threshold:
                # In this code what we do is: store the position of the residues
                # that aligns with the other chain. This is done because the
                # superposition needs to be done between two sequences with the
                # same number of atoms. If one residue is missmatched and we
                # take it, it will return an error.
                    seq1pos = 0
                    seq2pos = 0
                    list_seq1pos = []
                    list_seq2pos = []
                    # For each position in the alignment (including gaps).
                    for respos in range(0, len(sequence1)):
                        # If neither of the sequences are gaps, they are
                        # matching: store the positions where they match.
                        if sequence1[respos] != "-" and sequence2[respos] != "-":
                            list_seq1pos.append(seq1pos)
                            list_seq2pos.append(seq2pos)
                        # If sequence1 is not a gap, advance one position (this
                        # way, third "real" aminoacid (not gap) will be third
                        # position, fourth will be fourth...)
                        if sequence1[respos] != "-":
                            seq1pos += 1
                        # The same as sequence1
                        if sequence2[respos] != "-":
                            seq2pos += 1
                    # Store the information as an interaction
                    interacting_structures.append([chain1[:-1], chain1, list_seq1pos,
                                                  chain2[:-1], chain2, list_seq2pos])
        index += 1

    # If no interactions found, raise an error
    try:
        if len(interacting_structures)==0:
            raise not_interactions_found()
    except:
        print(not_interactions_found())
        sys.exit()

    return interacting_structures


def counting_interactions(interacting_structures, id_structure_dict):
    '''
    Counting the interactions each structure has with each other in order to
    take the most interacting structure as the template for the superimposition.
    '''
    count_interactions = {}
    for structure in id_structure_dict.keys():
    	count_interactions[structure] = 0
    	for interaction in interacting_structures:
    	    if structure in interaction:
    	        count_interactions[structure] += 1

    # Sort them by number of interactions
    sorted_count_interactions = sorted(count_interactions.items(),
                                       key=lambda x: x[1], reverse=True)

    return sorted_count_interactions


def preparing_randomization(args, id_structure_dict, sorted_count_interactions = False):
    '''
    This function creates a list of different seeds (numbers from 1 to 100000)
    based on an input seed given by the user (optional). This list will be used
    to randomize in each model created the order of the interactions added and,
    if -r specified, also the initial structure that will be taken for starting
    each new model.
    '''

    # Creating a list of random numbers (acting as seeds) based on the user
    # input seed. Each number will be used in each loop as a seed.
    number_of_models = int(args.models)
    random.seed(args.seed)
    list_of_seeds = random.sample(range(1,100000), number_of_models)

    # Put every structure in a list
    every_structure = []
    for structure in id_structure_dict.keys():
        every_structure.append(structure)

    # Create a list of the structures that will be used as initial structure in
    # the superimposition.
    starting_structure_id_list = []
    if args.random_seed:
        for model in range(0, number_of_models):
            random.seed(list_of_seeds[model])
            starting_structure_id_list.append(every_structure[random.randrange(0, len(every_structure))])
    else:
        for model in range(0, number_of_models):
            starting_structure_id_list.append(sorted_count_interactions[0][0])

    return(starting_structure_id_list, list_of_seeds)


def complex_builder(modeller, path, log, args, new_starting_structure, index,
                    ID_subunit_dict, id_structure_dict, interacting_structures,
                    pdb_list, list_of_seeds, recover_id_from_fasta,
                    sto_dic=False):
    '''
    This function takes an initial structure (new_starting_structure) and keeps
    adding new chains based on the interactions found before, using the pairwise
    alignment. Identical chains are superimosed and the resulting structure
    is checked for clashes between atoms. Chain names are changed
    for them to be unique.
    '''

    print("Building model nº" + str(index + 1))
    # We get the starting structure from the id_structure_dict
    starting_structure_id = [new_starting_structure]
    starting_structure = id_structure_dict[new_starting_structure]

    # We make a copy of the list of interacting structures to reset it every
    # time the loop starts (as we will remove interactions).
    interacting_structures_loop = copy.deepcopy(interacting_structures)

    # A seed from the list of seeds is used for shuffling the interactions.
    random.seed(list_of_seeds[index])
    random.shuffle(interacting_structures_loop)

    number_residues = 0
    if args.sto:
        # Initialize the model with the specific stoichiometry
        sto_structure = Bio.PDB.Model.Model('sto_model')
        # If the user wants to insert a stoichiometry, in sto_dic the subunits
        # are stored as keys and the values are the maximum number of times
        # that the subunits will appear
        sto_dic_loop = copy.copy(sto_dic)

    # Initialize the model object
    structure = Bio.PDB.Model.Model('model')
    n_chain_name = 0
    n_chain_name_sto = 0
    added_chains = []
    # We put inside the model the two first chains: the ones of the first interacting structure
    for chain in starting_structure:
        # We create a new chain object with the residues of the selected chain
        # and we give it a new ID
        id_select = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"
        id = id_select[n_chain_name:n_chain_name+1]
        n_chain_name += 1
        added_chain = Bio.PDB.Chain.Chain(id)
        added_chain.child_list=list(chain.get_residues())
        ID_subunit = recover_id_from_fasta[chain.id]
        structure.add(added_chain)

        # If stoichiometry is selected, the subunits are added to the model the
        # maximum number of times the user selects.
        if args.sto:
            # If the subunit is selected by the user with -sto or -stoinfo
            if ID_subunit in sto_dic_loop:
                # If there are still subunits to add of the specific stoichiometry
                if sto_dic_loop[ID_subunit] > 0:
                    id_sto = id_select[n_chain_name_sto:n_chain_name_sto+1]
                    n_chain_name_sto += 1
                    added_chain_sto = Bio.PDB.Chain.Chain(id_sto)
                    added_chain_sto.child_list = list(chain.get_residues())
                    sto_structure.add(added_chain_sto)
                    print("Subunit {} added".format(ID_subunit))
                    sto_dic_loop[ID_subunit] -= 1
                    added_chains.append(ID_subunit)
        else:
            print("Subunit {} added".format(ID_subunit))
            added_chains.append(ID_subunit)

    # Everytime an interaction is done and a chain is added (or not accomplished)
    # this loops, starting again looking in all remaining interactions.
    for interaction_found in range(0, len(pdb_list)):
        if interaction_found != 0 and not found:
            # Exit the main loop if no more interactions are found.
            break

        # Loop through every interaction stored
        for interaction in interacting_structures_loop:
            # Atoms that will be fixed in the superimposition
            fixed_atoms = []
            # Atoms that will move and rotate in the superimposition
            moving_atoms = []
            found = 0
            # If the first member of the interaction == the starting structure...
            if interaction[0] in starting_structure_id:
                fixed_chain = ID_subunit_dict[interaction[1]]                          # Chain that remains fixed
                fixed_residues_position = interaction[2]                               # Which residues of the chain are the ones that aligned
                moving_chain = ID_subunit_dict[interaction[4]]                         # Chain that is rotating
                moving_residues_position = interaction[5]                              # Which residues of the chain are the ones that aligned
                moving_structure = copy.deepcopy(id_structure_dict[interaction[3]])    # Copy of the moving structure (we make a copy to eliminate the chain that is superimposed (duplicated))
                moving_structure_id = interaction[3]                                   # ID of the moving structure (to add it to the starting_structure_id list)
                found = 1                                                              # An interaction was found
            # Else, if the second member of the interaction == the starting_structure structure...
            elif interaction[3] in starting_structure_id:
                fixed_chain = ID_subunit_dict[interaction[4]]
                fixed_residues_position = interaction[5]
                moving_chain = ID_subunit_dict[interaction[1]]
                moving_residues_position = interaction[2]
                moving_structure = copy.deepcopy(id_structure_dict[interaction[0]])
                moving_structure_id = interaction[0]
                found = 1

            if found:
                # We are going to store the atoms that aligned to make the
                # superposition.
                for pos, fixed_residue in enumerate(fixed_chain.get_residues()):
                    # If the aminoacid in that position did align, then store it
                    if pos in fixed_residues_position:
                        for atom in fixed_residue.get_atom():
                            # We store the alpha-carbon for aminoacids or
                            # phosphorus for DNA/RNA.
                            if atom.get_id() == "CA" or atom.get_id() == "P":
                                fixed_atoms.append(atom)
                for pos, moving_residue in enumerate(moving_chain.get_residues()):
                    if pos in moving_residues_position:
                        for atom in moving_residue.get_atom():
                            if atom.get_id() == "CA" or atom.get_id() == "P":
                                moving_atoms.append(atom)

                # Delete the chain that is being superimposed, as it is
                # duplicated (atoms are already stored).
                moving_structure.detach_child(moving_chain.id)

                super_imposer = Bio.PDB.Superimposer()

                # Make a rotation matrix of moving atoms over the fixed atoms
                # to minimize RMSE.
                super_imposer.set_atoms(fixed_atoms, moving_atoms)
                # Apply this rotation to the moving structure (there is only one
                # chain left, the not duplicated one).
                super_imposer.apply(moving_structure)

                # Now we will see if the superimposed structure will clash with
                # the starting structure.
                # We take the position of the starting structure atoms.
                neighbor = Bio.PDB.NeighborSearch(list(structure.get_atoms()))
                clashes = 0
                # For each atom in the moving structure we check if there are
                # close atoms to the starting structure.
                for atom in moving_structure.get_atoms():
                    close_atoms = neighbor.search(atom.get_coord(), float(args.clash_dist))
                    # If there are atoms within that distance, consider a clash.
                    if len(close_atoms) > 0:
                        clashes += 1
                # If we get more than 5 clashes, we consider its a clash and
                # continue to the next loop, aborting the superimposition of
                # the structure.
                if clashes > 5:
                    #print("Clash!")
                    found = 0
                    continue

                # If its not a clash, we get the remaining chain in the moving
                # structure and store it.
                for chain in moving_structure.get_chains():
                    # We loop through the alphabet looking for a not used letter
                    # for the new chain.
                    id = id_select[n_chain_name:n_chain_name + 1]
                    n_chain_name += 1
                    # A new chain object is created and the residues of the
                    # moving chain are copied to it.
                    added_chain = Bio.PDB.Chain.Chain(id)
                    added_chain.child_list = list(chain.get_residues())
                    # Here we get the ID of the subunit that has been added.
                    ID_subunit = recover_id_from_fasta[chain.id]
                # Then, we store the ID of the added structure to look for
                # the next interactions to superimpose.
                starting_structure_id.append(moving_structure_id)

                # If stoichiometry is selected, the subunits are added to the
                # model the maximum number of times the user selects.
                if args.sto:
                    # If the subunit is selected by the user with -sto or -stoinfo.
                    if ID_subunit in sto_dic_loop:
                        # If there are still subunits to add of the specific stoichiometry.
                        if sto_dic_loop[ID_subunit] > 0:
                            id_sto = id_select[n_chain_name_sto:n_chain_name_sto + 1]
                            n_chain_name_sto += 1
                            added_chain_sto = Bio.PDB.Chain.Chain(id_sto)
                            added_chain_sto.child_list = list(chain.get_residues())
                            added_chains.append(ID_subunit)
                            sto_structure.add(added_chain_sto)
                            # The dict containg the stoichiometry changes.
                            sto_dic_loop[ID_subunit] -= 1
                            print("Subunit {} added".format(ID_subunit))
                else:
                    # The info of the added subunit is stored.
                    added_chains.append(ID_subunit)
                    print("Subunit {} added".format(ID_subunit))

                # We add the new chain to the existing structure.
                structure.add(added_chain)
                interacting_structures_loop.remove(interaction)

                # We remove every interaction which subunits has already been added.
                for interaction_to_remove in interacting_structures_loop:
                    if interaction_to_remove[0] in starting_structure_id and \
                    interaction_to_remove[3] in starting_structure_id:
                        interacting_structures_loop.remove(interaction_to_remove)

                # We exit the loop, to start over again looking from the first
                # interaction now that the structure is bigger.
                break

    # If a specific stoichiometry is selected, store the specific model.
    if args.sto:
        store_results(modeller, path, log, args, structure, index, added_chains,
                      sto_structure)
    else:
        store_results(modeller, path, log, args, structure, index, added_chains)


def store_results(modeller, path, log, args, structure, index, added_chains,
                  sto_structure = False):
    '''
    This functions saves a PDB with the model created and, if specified, it
    optimizes the model with Modeller and saves a PDB with the optimized version.
    '''
    # Store the result.
    io = Bio.PDB.PDBIO()
    if args.sto:
        try:
            # Check if it is empty or not.
            if sto_structure:
                io.set_structure(sto_structure)
            # If it is empty an error appears.
            else:
                raise not_chains_added()
        except:
            print(not_chains_added())
            sys.exit()
    else:
        io.set_structure(structure)

    output_path = path + "/" + args.output + "_" + str(index + 1) + ".pdb"
    io.save(output_path)

    # Here we get the chains that have been added to the model to store later
    # the info.
    repeated = []
    added = ''
    added_chains = sorted(added_chains)
    for chain in added_chains:
        if chain not in repeated:
            added += chain + ':' + str(added_chains.count(chain)) + ','
            repeated.append(chain)

    # Optimization with modeller.
    # The info of the experiment is stored in the log file.
    output_filename = output_path.split('/')[-1]
    if modeller == True:
        from prediprot_optimize import Optimizemodel
        if args.optimize:
            init_energy, opt_energy = Optimizemodel(output_path, True)
            print("{:<35}\t{:>20.3f}\t{:<}".format(output_filename,float(init_energy),str(added[:-1])),
                  file = log)
            print("{:<35}\t{:>20.3f}\t{:<}".format(output_filename[:-4]+'_optimized.pdb',float(opt_energy),str(added[:-1])),
                  file = log)
        else:
            init_energy = Optimizemodel(output_path, False)
            print("{:<35}\t{:>20.3f}\t{:<}".format(output_filename,float(init_energy),str(added[:-1])),
                  file = log)
    else:
        print("{:<35}\t{:<}".format(output_filename,str(added[:-1])),
              file = log)
    print("Model nº" + str(index+1) + " completed")
