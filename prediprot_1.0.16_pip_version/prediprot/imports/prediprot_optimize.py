from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions
import sys
import os

def Optimizemodel(pdb_file, optimize):
    print("Calculating the energy of the model...")
    output_path_optimized = "/".join(pdb_file.split('/')[:-1])
    sys.stdout= open(os.devnull, 'w')
    env = environ()
    env.io.atom_files_directory = ['../atom_files']
    env.edat.dynamic_sphere = True
    env.libs.topology.read(file = '$(LIB)/top_heav.lib')
    env.libs.parameters.read(file = '$(LIB)/par.lib')
    filename = pdb_file.split("/")[-1]
    code = filename.split('.')[0]
    mdl = complete_pdb(env, pdb_file)
    mdl.write(file = output_path_optimized + "/" + code + '.ini')
	# Select all atoms:
    atmsel = selection(mdl)
    mpdf2 = atmsel.energy()
    sys.stdout=sys.__stdout__
    if optimize:
        print('Optimizing the model...')
        sys.stdout= open(os.devnull, 'w')
        # Generate the restraints:
        mdl.restraints.make(atmsel, restraint_type = 'stereo', spline_on_site = False)
        mdl.restraints.write(file = output_path_optimized + "/" + code + '.rsr')

		# Create optimizer objects and set defaults for all further optimizations
        cg = conjugate_gradients(output = 'REPORT')
        md = molecular_dynamics(output = 'REPORT')

		# Run CG on the all-atom selection
        cg.optimize(atmsel, max_iterations = 20)
		# Run MD
        md.optimize(atmsel, temperature = 300, max_iterations = 50)
		# Finish off with some more CG
        cg.optimize(atmsel, max_iterations = 20)

        mpdf = atmsel.energy()

		# Create .pdb
        mdl.write(file = output_path_optimized + "/" + code + '_optimized.pdb')

		# Remove interdmediate files
        os.remove(output_path_optimized + "/" + code + ".ini")
        os.remove(output_path_optimized + "/" + code + ".rsr")
        sys.stdout=sys.__stdout__
        return mpdf2[0], mpdf[0]

    else:
        os.remove(output_path_optimized + "/" + code + '.ini')
        return mpdf2[0]
