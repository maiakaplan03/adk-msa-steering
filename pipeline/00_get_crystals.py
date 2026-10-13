#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_get_crystals.py — "OBTENEMOS LOS CRISTALES EXPERIMENTALES".
Downloads 1AKE (closed, Ap5A-bound) and 4AKE (open) from the PDB and keeps chain A,
standard amino acids only (waters, ions and the Ap5A ligand are dropped).
Usage: python pipeline/00_get_crystals.py [--out-dir pdb_chainA_protein_only]
"""
import argparse
import os
from Bio.PDB import PDBList, PDBParser, PDBIO, Select, is_aa

PDB_IDS = ["1AKE", "4AKE"]
OUT_DIR = "pdb_chainA_protein_only"
CHAIN_ID = "A"

class ChainAProteinOnly(Select):
    def __init__(self, chain_id="A"):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return 1 if chain.id == self.chain_id else 0

    def accept_residue(self, residue):
        # Keep only standard amino acids (protein). This excludes waters/ligands/ions.
        return 1 if is_aa(residue, standard=True) else 0

    def accept_atom(self, atom):
        # Keep all atoms of kept residues (including hydrogens if present)
        return 1

def download_and_clean(pdb_id: str, out_dir: str, chain_id: str = "A"):
    os.makedirs(out_dir, exist_ok=True)

    pdbl = PDBList()
    downloaded_path = pdbl.retrieve_pdb_file(pdb_id, pdir=out_dir, file_format="pdb")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, downloaded_path)

    # Optional: sanity check chain exists in first model
    model = structure[0]
    chains = {c.id for c in model.get_chains()}
    if chain_id not in chains:
        raise ValueError(f"{pdb_id}: chain '{chain_id}' not found. Available: {sorted(chains)}")

    out_path = os.path.join(out_dir, f"{pdb_id}_chain{chain_id}_protein_only.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, select=ChainAProteinOnly(chain_id))

    return downloaded_path, out_path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=OUT_DIR)
    out_dir = ap.parse_args().out_dir
    for pdb_id in PDB_IDS:
        dl, outp = download_and_clean(pdb_id, out_dir, CHAIN_ID)
        print(f"{pdb_id}: downloaded -> {dl}")
        print(f"{pdb_id}: cleaned (chain {CHAIN_ID}, protein only) -> {outp}")

if __name__ == "__main__":
    main()
