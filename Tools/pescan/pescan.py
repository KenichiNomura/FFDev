import sys, os, time
import glob, math
import subprocess
from collections import OrderedDict


comment='''
$comment
	potential surface scanner
$end
'''

rem='''
$rem
   JOBTYPE  PES_SCAN
   UNRESTRICTED TRUE
   SYM_IGNORE  TRUE
   GEOM_OPT_MAX_CYCLES 400
   SCF_MAX_CYCLES 400
   METHOD   B3LYP
   BASIS   6-31G**
   PDB_PRINT 2
   SCF_CONVERGENCE 6
   THRESH 14
   SCF_PRINT 1
$end
'''

rem2='''
$rem
   JOBTYPE  PES_SCAN
   UNRESTRICTED TRUE
   SYM_IGNORE  TRUE
   GEOM_OPT_MAX_CYCLES 400
   SCF_MAX_CYCLES 400
   METHOD   B3LYP
   BASIS   6-31G**
   PDB_PRINT 2
   SCF_CONVERGENCE 6
   THRESH 10
   VARTHRESH 14
   INCDFT FALSE
   SCF_PRINT 1
$end
'''

def get_molecule_section(xyz_data, electron_state='singlet'):

	if electron_state == 'singlet':
		estate = '0 1'
	elif electron_state == 'triplet':
		estate = '0 3'
	else:
		print('\nWARNING : unsupported electron state {}. Use singlet instead.\n'.format(electron_state))
		estate = '0 1'

	molecule_section = estate + '\n'

	for n in xyz_data:
		atom_data = xyz_data[n]
		elem = atom_data[0]
		x = atom_data[1]
		y = atom_data[2]
		z = atom_data[3]
		
		#print(elem,x,y,z)
		line = "%s %f %f %f\n"%(elem, x, y, z)
		molecule_section += line 

	return '$molecule\n' + molecule_section + '$end\n' 


def load_xyz(xyz_file):
	
	with open(xyz_file,'r') as fin:
		# first line is the total number of atoms
		natoms = int(fin.readline())
	
		# second line is key words, "bond_scan 1 2" "angle_scan 1 2 3" "torsion_scan 1 2 3 4"
		scan_atom_index = {}
		data = fin.readline().split()
		while len(data) > 0:
			d = data.pop(0)
			if d == "bond":
				i = int(data.pop(0))
				j = int(data.pop(0))
				scan_atom_index[d] = (i,j)
			elif d == "angle":
				i = int(data.pop(0))
				j = int(data.pop(0))
				k = int(data.pop(0))
				scan_atom_index[d] = (i,j,k)
			elif d == "torsion":
				i = int(data.pop(0))
				j = int(data.pop(0))
				k = int(data.pop(0))
				l = int(data.pop(0))
				scan_atom_index[d] = (i,j,k,l)
	
		xyz_data={}
		for n in range(natoms):
			data = fin.readline().split()
			elem = data[0]
			x = float(data[1])
			y = float(data[2])
			z = float(data[3])
			xyz_data[n]=[elem,x,y,z]

	return xyz_data, scan_atom_index


def get_scan_initial_values(xyz_data, scan_atom_index):
	print(scan_atom_index)

	initial_values = {}
	for key in scan_atom_index:
		if key == 'bond':
			i = scan_atom_index[key][0] - 1  # zero-indexed in python
			j = scan_atom_index[key][1] - 1  # zero-indexed in python
			iatom = xyz_data[i]
			jatom = xyz_data[j]

			r2 = 0.0
			for ia in range(1,4):
				dr = iatom[ia]-jatom[ia]
				r2 += dr*dr
			distance = math.sqrt(r2)	

			distance = 0.1*int(distance*10) # up to first decimal point, e.g. 0.9, 1.2, 1.8 etc
			#print(key, distance, i,j,iatom,jatom)
			initial_values[key] = distance
		else:
			initial_values[key] = None
			print('scan keyword ', key, 'is not supported in get_scan_initial_values')

	return initial_values

def get_bondscan_range(scan_atom_index, initial, final):

	increment = DEFAULT_INCREMENT

	if initial < final:
		incr = increment
	else:
		incr = -increment

	num_datapoints = int(abs(final - initial)/increment)
	print('scan_atom_index, initial, final, increment, num_datapoints : ', \
			scan_atom_index, initial, final, increment, num_datapoints)

	# get atom index pair for bond-scan
	i = scan_atom_index[0]
	j = scan_atom_index[1]

	bondscan_strings = {}
	for idx in range(num_datapoints):
		v1 = initial + idx*incr
		v2 = initial + idx*incr+1e-3
		start = min(v1,v2)
		end = max(v1,v2)

		#scanid = int(v1*1000) # unique ID number 
		scanid = '{:03.1f}'.format(v1)

		#print(scanid,idx,start,end,initial)

		bondscan_strings[scanid] = "stre %d %d %f %f %f"%(i,j,start,end,abs(incr))

	return bondscan_strings

def save_xyz(filename, xyz_data):

	with open(filename, 'w') as fout:
		fout.write('{}\n{}\n'.format(len(xyz_data), PESCAN_GENERATED))
		for key, p in xyz_data.items():
			fout.write('{} {} {} {}\n'.format(p[0],p[1],p[2],p[3]))

def get_scan_result(qchem_outfile):

	key1 = '------- Summary of potential scan: ------'

	with open(qchem_outfile, 'r') as fin:

		print('In get_scan_result with {}'.format(qchem_outfile))

		for line in fin:

			if key1 in line:
				data = fin.readline().split()
				coord = float(data[0])
				energy = float(data[1])
			
				return coord, energy

	print('WARNING: keyword {} was not found'.format(key1))

def get_optimized_sturcture(qchem_outfile):

	key1 = '**  OPTIMIZATION CONVERGED  **'	
	key2 = 'Coordinates (Angstroms)' 
	key3 = 'ATOM                X               Y               Z' 

	xyz_data = {}

	with open(qchem_outfile, 'r') as fin:

		print('In get_optimized_sturcture with {}'.format(qchem_outfile))

		for line in fin:
			if key1 in line:

				# after key1 is found, process header lines
				print("keyword found : {}".format(key1))

				header = fin.readline()
				header = fin.readline()

				header = fin.readline()
				assert key2 in header, "ERROR : keyword {} was found but not {} in {}".format(key1, key2, header)

				header = fin.readline()
				assert key3 in header, "ERROR : keyword {} was found but not {} in {}".format(key1, key3, header)

				while True:
					line = fin.readline()
					data = line.split()

					if len(data) > 0:
						index = int(data[0])
						elem = data[1]
						x = float(data[2])
						y = float(data[3])
						z = float(data[4])
						xyz_data[index]=[elem, x, y, z]
					else:
						# finish loading the optimized structure
						#for key in xyz_data:
						#	print(xyz_data[key])
						return xyz_data

	if len(xyz_data) == 0:
		print('WARNING : optimized structure was not found')


	return xyz_data

def run_PSS_bond(pss_data, scan_atom_index, base_xyz, electron_state='singlet', initial=1.0, final=1.0, update_coordinate=True):

	start = time.time()

	molecule = get_molecule_section(base_xyz, electron_state)

	bondscan_strings = get_bondscan_range(scan_atom_index, initial, final)

	coord = initial
	xyz_opt_save = {}

	for key in bondscan_strings:

		scan = '$scan\n' + bondscan_strings[key] + '\n$end\n'
		#print(scan)

		scan_basestr = 'bond_' + electron_state + '_' + str(key)

		qchem_input_filename = os.path.join(input_dir, scan_basestr + ".in")
		qchem_output_filename = os.path.join(output_dir, scan_basestr + ".out")

		print('\n--- {} {} ---\n'.format(qchem_input_filename, qchem_output_filename) )

		with open(qchem_input_filename,'w') as infile:
			infile.write(comment)
			infile.write(rem)
			infile.write(molecule)
			infile.write(scan)

		command = ['qchem', '-nt', '12', qchem_input_filename, qchem_output_filename]
		result = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

		xyz_opt = get_optimized_sturcture(qchem_output_filename)
		if len(xyz_opt) > 0:
			save_xyz(os.path.join(xyz_dir, scan_basestr + ".xyz"), xyz_opt)

			coord, energy = get_scan_result(qchem_output_filename)
			pss_data[coord] = energy
			xyz_opt_save = xyz_opt
		else:
			break
		
		if update_coordinate:
			molecule = get_molecule_section(xyz_opt, electron_state)

	end = time.time()

	return pss_data, xyz_opt_save, coord, end - start

def check_if_pescangenerated_xyz(filename):
	with open(filename,'r') as fin:
		for line in fin:
			if PESCAN_GENERATED in line:
				return True

	return False


input_dir = 'inputs'
output_dir = 'outputs'
xyz_dir = 'xyz'

DEFAULT_INCREMENT = 0.1

PESCAN_GENERATED = 'generated by pescan.py'

xyz_file=sys.argv[1]

if check_if_pescangenerated_xyz(xyz_file):
	print('WARNING: {} is generated by pescan.py')	
	sys.exit(0)

base_dir = xyz_file+'.data'
input_dir = os.path.join(base_dir,input_dir)
output_dir = os.path.join(base_dir,output_dir)
xyz_dir = os.path.join(base_dir,xyz_dir)

if not os.path.exists(input_dir):
    os.makedirs(input_dir)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(xyz_dir):
    os.makedirs(xyz_dir)

xyz_data, scan_atom_index = load_xyz(xyz_file)
print('--------------')
for key in scan_atom_index:
	print('scan_atom_index: ', key, scan_atom_index[key])

print('--------------')
initial_values = get_scan_initial_values(xyz_data,scan_atom_index)

print('--------------')
for key in initial_values:
	print('initial_values: ', key, initial_values[key])

print('--------------')

bondscan_data = {}

electron_state = 'singlet'
bondlength_opt = initial_values["bond"]
min_cutoff = bondlength_opt*0.5
max_cutoff = max(4.0,bondlength_opt*3.0)
bondlength_min = bondlength_opt*0.4
bondlength_max = bondlength_opt*3.1

pss_index = scan_atom_index["bond"]


# singlet inner scan
print('\n=== 1. inner scan === \n')

print('\n--- start incremental inner scan ---\n')
bondscan_data, xyz_opt, coord_opt, t_elapsed = \
	run_PSS_bond(bondscan_data, pss_index, xyz_data, electron_state, bondlength_opt, bondlength_min, True)

if coord_opt > min_cutoff:
	print('\n--- final coord {}. start descrete inner scan ---\n'.format(coord_opt))
	bondscan_data, xyz_opt, coord_opt, t_elapsed = \
		run_PSS_bond(bondscan_data, pss_index, xyz_data, electron_state, bondlength_opt, bondlength_min, False)


# singlet outer scan
print('\n=== 2. outer scan === n')

print('\n--- start incremental outer scan ---\n')
bondscan_data, xyz_opt, coord_opt, t_elapsed = \
	run_PSS_bond(bondscan_data, pss_index, xyz_data, electron_state, bondlength_opt, bondlength_max, True)

if coord_opt < max_cutoff:
	print('\n--- start descrete outer scan ---\n')
	bondscan_data, xyz_opt, coord_opt, t_elapsed = \
		run_PSS_bond(bondscan_data, pss_index, xyz_data, electron_state, coord_opt, bondlength_max, False)

with open(os.path.join(base_dir, 'singlet.dat'),'w') as fout:
	for key in sorted(bondscan_data):
		fout.write('{} {}\n'.format(key, bondscan_data[key]))

# triplet scan
print('\n=== 3. triplet scan start === n')

electron_state = 'triplet'

xyz_files = glob.glob(xyz_dir+'/*.xyz')
xyz_files.sort()
xyz_opt, _ = load_xyz(xyz_files[len(xyz_files)-1])

bondscan_data = {}
bondscan_data, xyz_opt, coord_opt, t_elapsed = \
	run_PSS_bond(bondscan_data, pss_index, xyz_opt, electron_state, max_cutoff, bondlength_opt, True)

with open(os.path.join(base_dir, 'triplet.dat'),'w') as fout:
	for key in sorted(bondscan_data):
		fout.write('{} {}\n'.format(key, bondscan_data[key]))

print('\n=== pescan completed ===n')
