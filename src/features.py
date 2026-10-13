# -*- coding: utf-8 -*-
"""
分子指纹提取与基础特征计算模块
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski

NBITS = 2048
RADIUS = 2

def smiles_to_ecfp4(smiles, n_bits=NBITS, radius=RADIUS):
    """
    计算分子的 ECFP4 (Morgan 半径=2, 2048位) 二进制指纹，忠实于文献与基准设定
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    bv = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros(n_bits, dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(bv, arr)
    return arr

def batch_smiles_to_ecfp4(smiles_list, n_bits=NBITS, radius=RADIUS):
    """
    批量计算 SMILES 列表的 ECFP4 指纹矩阵
    """
    X = np.zeros((len(smiles_list), n_bits), dtype=np.uint8)
    valid_mask = np.ones(len(smiles_list), dtype=bool)
    for i, s in enumerate(smiles_list):
        fp = smiles_to_ecfp4(s, n_bits=n_bits, radius=radius)
        if fp is not None:
            X[i] = fp
        else:
            valid_mask[i] = False
    return X, valid_mask

def calc_physicochemical_properties(smiles):
    """
    计算基础理化属性：MW, LogP, TPSA, RotB, HBD, HBA
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "rotb": Lipinski.NumRotatableBonds(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
    }
