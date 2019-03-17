from __future__ import print_function
import numpy as np # numpy for managing arrays and basic math
import re
import collections # For OrderedDicts
import sys
class generator:

    def __init__(self):
        # Defining variables etc
        self.optxyz = []
        self.charge = []
        self.scan_energies = collections.OrderedDict()
        self.bond_stretch_xyz = collections.OrderedDict()
        self.energy_diff = collections.OrderedDict()
        self.other_param = collections.OrderedDict()
        self.weights = collections.OrderedDict()
        
    # This function extracts information from the file between two lines. We 
    # will use this function many times during the process to extract information
    # decided to write a function for it instead of coding this over and over
    @staticmethod
    def extractor(lines, marker_start, marker_end):
        counter = False
        a = 0
        b = [] #list that will contain starting indices of data in master lines list
        c = [] #list that contins endf point of indices 
        for line in lines:
            if(marker_start in line):
                counter = True
                b.append(a)
            elif(marker_end in line):
                counter = False
                c.append(a)
            a += 1
        d = [(i,j) for (i,j) in zip(b,c)]
        return d
    
    
    def read_files(self, opt_file, stretch_file):
        # Read those files
        # Return opt structure xyz coords as a list (float elements) - optxyz - [[1, 'atom_type1', X, Y, Z], [2, 'atom_type2', X, Y, Z]...]
        # Return opt mulliken charge as a list (float elements) - charge - [[1, 'atomtype1', charge], [2, 'atomtype2', charge]...]
        # Return length scan energies as a dict - scan_energies - {length1:energy1, length2:energy2, ...}
        # Return bond stretch xyz files as a dict - bond_stretch_xyz -  {length1:[[atom1, X, Y, Z], [atom2, X, Y, Z]...], length2:[[atom1, X, Y, Z], [atom2, X, Y, Z]...]...}
        # Return difference between scan energies and opt energy as a dict - energy_diff - {length1:e_diff, length2:e_diff, ...}
        # Return other parameters such as number of atoms, type of atoms, min and max stretch length etc as a dict - other_param - {} 
        
        # This cleans up the string i/p from bond stretch output file 
        with open(stretch_file, 'r') as fin_stretch:
            raw = fin_stretch.readlines()
            lines_stretch = [re.sub(r'\s\s+', ' ', element.strip()) for element in raw] # This is better than strip, this takes care of excessive spacing between first and last char in a line of a file
            lines_stretch = filter(None, lines_stretch) 
        
        # This cleans up the string i/p from bond optimization output file 
        with open(opt_file, 'r') as fin_opt:
            raw = fin_opt.readlines()
            lines_opt = [re.sub(r'\s\s+', ' ', element.strip()) for element in raw] # This is better than strip, this takes care of excessive spacing between first and last char in a line of a file
            lines_opt = filter(None, lines_opt) 
        opt_index = self.extractor(lines_opt, '** OPT', 'Z-matrix')
        stretch_index = self.extractor(lines_stretch, '** OPT', 'Z-matrix')
        # This operation turns this parsed info from the file into manipulatble lists
        # Operation done on opt structure
        unsplit = []
        for i in opt_index:
            unsplit.append(lines_opt[i[0] + 4:i[1]])
        for i in unsplit[0]:
            self.optxyz.append(i.split())
        for i in self.optxyz:
            i[0] = int(i[0])
            i[2] = float(i[2])
            i[3] = float(i[3])
            i[4] = float(i[4])

        # Same operation for bond stretch xyz
        unsplit = []
        for i in stretch_index:
            unsplit.append(lines_stretch[i[0] + 4:i[1]])
        test1 = []
        for i in unsplit:
            test = []
            for j in i:
                j = j.split()
                j[0] = int(j[0])
                j[2] = float(j[2])
                j[3] = float(j[3])
                j[4] = float(j[4])
                test.append(j)
            test1.append(test)     
         
        # Extracting bond stretch energy values as a dict
        energy_scan_index = self.extractor(lines_stretch, '--- Summary', 'Archival')
        unsplit = []
        energy_list = []
        for i in energy_scan_index:
            unsplit.append(lines_stretch[i[0] + 1:i[1] - 1])
        for i in unsplit[0]:
            energy_list.append(float(i.split()[0]))
            self.scan_energies.update({float(i.split()[0]):float(i.split()[1])})
        
        # Updating bond stretch energies and the correspoding xyz as a dict
        counter = 0
        for i in energy_list:
            self.bond_stretch_xyz.update({i:test1[counter]})
            counter += 1
        
        # Extracting the energy of the optimized structure and storing in other_params
        for line in lines_opt:
            if('Final energy is' in line):
                self.other_param.update({'E_opt':float(line.split()[-1])})
        
        # Constructing the differences in energy between stretch energies and optimized energy as a dict converted to kcal/mol 
        for key in self.scan_energies:
            self.energy_diff.update({key: ((self.scan_energies[key] - self.other_param['E_opt'])*627.51)})
        
        # Extracting the mulliken charge population.
        # Keep in mind that this outputs more than 1 tuple. The last one is the optimized one; we take that 
        charge_index = self.extractor(lines_opt, 'Ground-State Mulliken', 'Sum of atomic charges')
        for i in lines_opt[charge_index[-1][0] + 3 :charge_index[-1][1] - 1]:
            self.charge.append(i.split())
        
        # Extracting other parameters necessary for writing output files
        # The two atoms involved in bond stretching, stretch start and stretch end; among others
        for line in lines_stretch:
            if('stre' in line):
                #print line.split()
                self.other_param.update({'Atomnumber1': int(line.split()[1])})
                self.other_param.update({'Atomnumber2': int(line.split()[2])})
                self.other_param.update({'Stretch_Start': float(line.split()[3])})
                self.other_param.update({'Stretch_End': float(line.split()[4])})
                self.other_param.update({'Steps': float(line.split()[5])})
                break
            
        for i in self.optxyz:
            if(self.other_param['Atomnumber1'] == i[0]):
                self.other_param.update({'Atomtype1': i[1]})
            elif(self.other_param['Atomnumber2'] == i[0]):
                self.other_param.update({'Atomtype2': i[1]})


    def bond_l_calc(self):
        # Returns the bond length between the specific two atoms - bond_l
        # Returns the bonf angle between the two specific atoms - bond_theta
        for i in self.optxyz:
            if(self.other_param['Atomnumber1'] == i[0]):
                a1 = i[2:5]
            elif(self.other_param['Atomnumber2'] == i[0]):
                a2 = i[2:5]

        bond_length = np.linalg.norm([m - n for m,n in zip(a1,a2)])
        return bond_length
    
    def weight_calc(self):
        # Generates the weights for energies
        # Weights for charge and geometry hard coded
        
        old_range = max(self.energy_diff.values()) - min(self.energy_diff.values())
        new_range = 10.0
        new_min = 0.1
        for i in self.energy_diff:
            weight = 1.0*(((self.energy_diff[i]-min(self.energy_diff.values()))*new_range/old_range) + new_min)
            self.weights.update({weight:[i, self.energy_diff[i]]})
        return self.weights
    
    
    
    def generatefiles(self):
        # General string manipulations.
        # Nothing big to explain
        #550 format ('HETATM',1x,i5,1x,a2,3x,1x,3x,1x,1x,1x,5x,3f10.5,1x,$a5,i3,i2,1x,f8.5)
        f = open('geo', 'w')
        f.write('BIOGRF 200 \n')
        f.write('DESCRP ' + self.other_param['Atomtype2'] + '_' + self.other_param['Atomtype1'] + '_opt' + '\n')
        f.write('REMARK File generated by script\n')
        for i in self.optxyz:
            print('HETATM' + ' ' + str('%5d' % i[0]) + ' ' + i[1].ljust(2) + ' '*15 + str('%10.5f' % i[2]) +  str('%10.5f' % i[3]) + str('%10.5f' % i[4]) + ' '+ i[1].ljust(5) + '  0 0 0.00000',end = '\n', file = f )   
        f.write('END\n')
        for i in self.bond_stretch_xyz:
            f.write('\n')
            f.write('BIOGRF 200 \n')
            f.write('DESCRP ' + self.other_param['Atomtype2'] + '_' + self.other_param['Atomtype1'] + str(int(i*100)) + '\n')
            f.write('REMARK File generated by script\n')
            #BOND RESTRAINT    1  2  2.4500 7500.00  5.0000  0.0000000       0       0
            f.write('BOND RESTRAINT    ' + str(self.other_param['Atomnumber1']) + '  '+ str(self.other_param['Atomnumber2']) + '  ' + str((i)) + ' 7500.00  5.0000  0.0000000       0       0\n')
            #f.write()
            for j in self.bond_stretch_xyz[i]:
                print('HETATM' + ' ' + str('%5d' % j[0]) + ' ' + j[1].ljust(2) + ' '*15 + str('%10.5f' % j[2]) +  str('%10.5f' % j[3]) + str('%10.5f' % j[4]) + ' '+ j[1].ljust(5) + '  0 0 0.00000',end = '\n', file = f )   
            f.write('END\n')
        f.close()
        
        #Writing Geometry. 0.2 is weight; hardcoded
        f = open('trainset.in', 'w')
        f.write('GEOMETRY\n')
        f.write(self.other_param['Atomtype2']+ '_'+ self.other_param['Atomtype1'] + '_opt'+' 0.2 ' + str(self.other_param['Atomnumber1']) + ' ' + str(self.other_param['Atomnumber2']) + ' ' + str('%.3f' % self.bond_l_calc()))
        f.write('\nENDGEOMETRY\n')
        
        # Writing Charge, 0.01 is weight; hardcoded
        f.write('\nCHARGE\n')
        for i in self.charge:
            f.write(self.other_param['Atomtype2'] + '_' + self.other_param['Atomtype1'] + '_opt' + ' 0.01 ' + str(i[0]) + ' ' + str(i[2]) + '\n')
        f.write('ENDCHARGE\n')
        
        # Writing Energies with weights generated from the weight function
        f.write('\nENERGY\n')
        for i in self.weight_calc():
            print (str('%.4f' % i), '+', self.other_param['Atomtype2']+'_'+ self.other_param['Atomtype1'] + str(int(self.weight_calc()[i][0]*100)), '/1', '-', self.other_param['Atomtype2'] + '_' + self.other_param['Atomtype1'] + '_opt', '/1', str('%.4f' % self.weight_calc()[i][1]), sep = ' ', end = '\n', file = f)
        f.write('ENDENERGY')
        f.close()

        # Generates two files geo, trainset.in 

g = generator()
g.read_files(sys.argv[1], sys.argv[2])
g.generatefiles()
