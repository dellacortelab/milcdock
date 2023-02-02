
#################################################################################################
# Download a tranche from the ZINC database and separate the combined SDF into individual files
#################################################################################################

import os
import argparse
import shutil
import os

parser = argparse.ArgumentParser()
parser.add_argument('--data_dir', type=str, default='./data/0_ligand_receptor/ligand')
parser.add_argument('--zinc_download_script_file', type=str, default='ZINC-downloader-3D-mol2.gz.wget')
args = parser.parse_args()

def split_mol2_file(multi_mol2_file, output_dir):
    with open(multi_mol2_file, 'r') as in_mol2file:
        line = in_mol2file.readline()
        while not in_mol2file.tell() == os.fstat(in_mol2file.fileno()).st_size:
            if line.startswith("@<TRIPOS>MOLECULE"):
                mol2cont = []
                mol2cont.append(line)
                line = in_mol2file.readline()
                molecule_id = line.strip()
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol2cont.append(line)
                    line = in_mol2file.readline()
                    if in_mol2file.tell() == os.fstat(in_mol2file.fileno()).st_size:
                        mol2cont.append(line)
                        break
                mol2cont[-1] = mol2cont[-1].rstrip() # removes blank line at file end
                with open(os.path.join(output_dir, molecule_id + ".mol2"), 'w') as out_mol2file:
                    out_mol2file.write("".join(mol2cont))

def gz_to_mol2(tranches_dir, ligands_dir):
    subdirs = [f.path for f in os.scandir(tranches_dir) if f.is_dir()]
    for subdir in subdirs:
        subsubdirs = [f.path for f in os.scandir(subdir) if f.is_dir()]
        for subsubdir in subsubdirs:
            for zipfile in os.listdir(subsubdir):
                if '.gz' not in zipfile:
                    continue
                zipfile = os.path.join(subsubdir, zipfile)
                mol2_file = ''.join(os.path.splitext(zipfile)[:-1])
                if not os.path.exists(mol2_file):
                    os.system(f'gunzip {zipfile}')
                # Save molecules to individual files
                split_mol2_file(mol2_file, ligands_dir)


if __name__ == "__main__":
    data_dir = os.path.abspath(args.data_dir)
    tranches_dir = os.path.join(data_dir, 'zinc_tranches')
    ligands_dir = os.path.join(data_dir, 'zinc_ligands')
    os.makedirs(tranches_dir, exist_ok=True)
    os.makedirs(ligands_dir, exist_ok=True)
    os.system(f'cd {tranches_dir} && bash {os.path.join(data_dir, args.zinc_download_script_file)}')
    gz_to_mol2(tranches_dir=tranches_dir, ligands_dir=ligands_dir)
    shutil.rmtree(f'{tranches_dir}')