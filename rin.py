#!/usr/bin/env python3

from argparse import ArgumentParser
from os import path

from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.NeighborSearch import NeighborSearch


def parse_chains(structure_file, model_id=0, chain_ids=None):
    if structure_file.endswith(".cif"):
        parser = FastMMCIFParser(QUIET=True)
        structure_id = path.basename(structure_file).rstrip(".cif")
    model = parser.get_structure(structure_id, structure_file)[model_id]
    if chain_ids is not None:
        chains = [model[chain_id] for chain_id in chain_ids]
    else:
        chains = model.get_list()
    return chains


def extract_node_atoms(chains, node_atom_selection="CA", include_hetatms=False):
    node_atoms = []
    for chain in chains:
        for residue in chain:
            if node_atom_selection == "CA":
                if is_aa(residue) and residue.has_id(node_atom_selection):
                    node_atoms.append(residue[node_atom_selection])
            elif node_atom_selection == "CB":
                if is_aa(residue) and (residue.get_resname() != "GLY") and residue.has_id(node_atom_selection):
                    node_atoms.append(residue[node_atom_selection])
                elif is_aa(residue) and (residue.get_resname() == "GLY") and residue.has_id("CA"):
                    node_atoms.append(residue["CA"])
            elif node_atom_selection == "nonhydrogen":
                if not include_hetatms and residue.get_id()[0].strip():
                    continue
                for atom in residue:
                    if atom.element != "H":
                        node_atoms.append(atom)
    return node_atoms


def find_atom_contacts(node_atoms, cutoff=8.0):
    ns = NeighborSearch(node_atoms)
    atom_contacts = ns.search_all(cutoff)
    return atom_contacts


def _main(structure_file, model_id, chain_ids, node_atom_selection, include_hetatms, cutoff):
    chains = parse_chains(structure_file, model_id, chain_ids)
    node_atoms = extract_node_atoms(chains, node_atom_selection, include_hetatms)
    atom_contacts = find_atom_contacts(node_atoms, cutoff)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("structure_file")
    parser.add_argument("--model_id", default=0, type=int)
    parser.add_argument("--chain_ids", nargs="+")
    parser.add_argument("--node_atom_selection", default="CA")
    parser.add_argument("--include_hetatms", action="store_true")
    parser.add_argument("--cutoff", default=8.0, type=float)
    args = vars(parser.parse_args())
    _main(**args)
