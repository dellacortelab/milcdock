import MDAnalysis as mda
import sys

if len(sys.argv) == 1:
    print("Script to convert atom names that aren't unique to unique names that make Vina and Autodock output atoms comparable with other program outputs. ")
    print("Arg 1: path to mol2 file to change atom names of.")
    print("Will output a new file named {original_file_name}_renamed.mol2")
    sys.exit()

path = sys.argv[1]
uni = mda.Universe(path)

print('Fixing atom names')
atom_types=dict()
for i,atom in enumerate(uni.select_atoms('not type H')):
    old_name=atom.name
    if atom.name not in atom_types.keys():
        atom_types[atom.name]=1
    if not atom.name.endswith('*') and not atom.name == 'OXT' and not atom.name[-1].isdigit():
        atom.name = atom.name + str(atom_types[atom.name])
    atom_types[old_name]+=1

new_name = path[:path.rfind('.mol2')] + '_renamed.mol2'
print(f'Writing fixed atom names to {new_name}')
uni.atoms.write(new_name)