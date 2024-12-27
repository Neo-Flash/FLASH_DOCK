import os
import pandas as pd
import torch
import numpy as np
from rdkit import Chem
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm  # 进度条库
from PLANET_model import PLANET
from chemutils import ProteinPocket, mol_batch_to_graph
import rdkit.Chem as Chem
import sys
from rdkit import RDLogger
import argparse

class PlanetEstimator:
    def __init__(self, device):
        self.model = PLANET(300, 8, 300, 300, 3, 10, 1, device=device)  # trained PLANET
        self.model.load_parameters()
        for param in self.model.parameters():
            param.requires_grad = False
        self.model.eval()

    def set_pocket_from_ligand(self, protein_pdb, ligand_sdf):
        try:
            self.pocket = ProteinPocket(protein_pdb=protein_pdb, ligand_sdf=ligand_sdf)
        except:
            raise RuntimeError('The protein PDB file needs to be fixed')
        self.res_features = self.model.cal_res_features_helper(self.pocket.res_features, self.pocket.alpha_coordinates)

    def set_pocket_from_coordinate(self, protein_pdb, centeriod_x, centeriod_y, centeriod_z):
        try:
            self.pocket = ProteinPocket(protein_pdb, centeriod_x, centeriod_y, centeriod_z)
        except:
            raise RuntimeError('The protein PDB file needs to be fixed')
        self.res_features = self.model.cal_res_features_helper(self.pocket.res_features, self.pocket.alpha_coordinates)

    def pre_cal_res_features(self):
        self.res_features = self.model.cal_res_features_helper(self.pocket.res_features, self.pocket.alpha_coordinates)

class SingleMoleculeDataset(Dataset):
    def __init__(self, sdf_file):
        self.sdf_supp = Chem.SDMolSupplier(sdf_file, removeHs=False, sanitize=True)
        self.mol = None
        for mol in self.sdf_supp:
            if mol is not None:
                self.mol = mol
                break
        if self.mol is None:
            raise ValueError("No valid molecule found in the provided SDF file.")

    def __len__(self):
        return 1

    def __getitem__(self, idx):
        return self.tensorize()

    def tensorize(self):
        try:
            mol = Chem.AddHs(self.mol)
            mol_feature = mol_batch_to_graph([mol])
            mol_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
            mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else 'Unnamed'
            return (mol_feature, mol_smiles, mol_name)
        except:
            return (None, None, None)

def workflow(protein_pdb, ligand_sdf, mol_file, device):
    estimator = PlanetEstimator(device)
    estimator.model.to(device)
    estimator.set_pocket_from_ligand(protein_pdb, ligand_sdf)

    dataset = SingleMoleculeDataset(mol_file)
    dataloader = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=0, drop_last=False, collate_fn=lambda x: x[0])

    predicted_affinities, mol_names, smis = [], [], []
    with torch.no_grad():
        for (mol_feature, smi, mol_name) in dataloader:
            if mol_feature is None:
                print(f"Skipping molecule {mol_name} due to processing error.")
                continue
            batch_size = 1
            fresidues_batch, res_scope = estimator.model.cal_res_features(estimator.res_features, batch_size)
            predicted_affinity = estimator.model.screening(fresidues_batch, res_scope, mol_feature)
            predicted_affinities.append(predicted_affinity.view([-1]).cpu().numpy())
            smis.append(smi)
            mol_names.append(mol_name)

    predicted_affinities = np.concatenate(predicted_affinities) if predicted_affinities else np.array([])
    return predicted_affinities, mol_names, smis

def result_to_csv(protein_pdb, predicted_affinities, mol_names, smis, csv_file):
    protein_name = os.path.basename(protein_pdb)  # 仅获取文件名
    csv_frame = pd.DataFrame([
        {
            'Protein_Name': protein_name,
            'Ligand_Smiles': smi,
            'Binding_Affinity': aff,
        }
        for aff, smi in zip(predicted_affinities, smis)
    ])
    csv_frame.to_csv(csv_file, mode='a', header=not os.path.exists(csv_file), index=False)

def single_process(protein_pdb, ligand_sdf, mol_file, output_csv):
    if torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
    predicted_affinities, mol_names, smis = workflow(protein_pdb, ligand_sdf, mol_file, device)
    if len(predicted_affinities) > 0:
        print(f"Saving results for {protein_pdb} and {mol_file}")
        result_to_csv(protein_pdb, predicted_affinities, mol_names, smis, output_csv)
        print(f"Finished processing {protein_pdb} and {mol_file}")
    else:
        print("No affinities were predicted.")

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')  # 禁用RDKit日志
    parser = argparse.ArgumentParser(description="Compute binding affinity between a protein and a ligand.")
    parser.add_argument('-p', '--protein', required=True, help="Protein PDB file")
    parser.add_argument('-l', '--ligand', required=True, help="Ligand SDF file (used as both ligand and molecule)")
    parser.add_argument('-m', '--mol_file', required=True, help="Molecule SDF file (same as ligand)")
    parser.add_argument('-o', '--output_csv', required=True, help="Output CSV file")

    args = parser.parse_args()

    # 验证ligand和mol_file是否相同
    if os.path.abspath(args.ligand) != os.path.abspath(args.mol_file):
        print("Warning: The ligand SDF file and molecule SDF file are different. They should be the same as per requirements.")
    
    single_process(args.protein, args.ligand, args.mol_file, args.output_csv)
    print("Processing completed!")
