#!/usr/bin/env python3

from argparse import ArgumentParser
from os import path

from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.NeighborSearch import NeighborSearch
import networkx as nx

_NONPOLAR_ATOMS = {"C"}


def _remove_hydrogens(chains):
    for chain in chains:
        for residue in chain:
            hydrogen_atom_ids = [atom.get_id() for atom in residue if atom.element == "H"]
            for hydrogen_atom_id in hydrogen_atom_ids:
                residue.detach_child(hydrogen_atom_id)
        assert [atom for atom in chain.get_atoms() if atom.element == "H"] == []


def parse_chains(structure_file, model_id=0, chain_ids=None):
    if structure_file.endswith(".cif"):
        parser = FastMMCIFParser(QUIET=True)
    structure_id = path.basename(structure_file)[:-4]
    model = parser.get_structure(structure_id, structure_file)[model_id]
    if chain_ids is not None:
        chains = [model[chain_id] for chain_id in chain_ids]
    else:
        chains = model.get_list()
    _remove_hydrogens(chains)
    return chains


def extract_node_atoms(chains, node_atom_selection="CA", include_hetatms=False):
    residue_ids = []
    node_atoms = []
    for chain in chains:
        for residue in chain:
            if node_atom_selection == "CA":
                if is_aa(residue) and residue.has_id(node_atom_selection):
                    residue_ids.append(residue.get_full_id())
                    node_atoms.append(residue[node_atom_selection])
            elif node_atom_selection == "CB":
                if is_aa(residue) and (residue.get_resname() != "GLY") and residue.has_id(node_atom_selection):
                    residue_ids.append(residue.get_full_id())
                    node_atoms.append(residue[node_atom_selection])
                elif is_aa(residue) and (residue.get_resname() == "GLY") and residue.has_id("CA"):
                    residue_ids.append(residue.get_full_id())
                    node_atoms.append(residue["CA"])
            elif node_atom_selection == "nonhydrogen":
                if not include_hetatms and residue.get_id()[0].strip():
                    continue
                residue_ids.append(residue.get_full_id())
                for atom in residue:
                    if atom.element != "H":
                        node_atoms.append(atom)
    return residue_ids, node_atoms


def find_atom_contacts(node_atoms, cutoff=8.0):
    ns = NeighborSearch(node_atoms)
    atom_contacts = ns.search_all(cutoff)
    return atom_contacts


def _get_interaction_type(atom_i, atom_j):
    if (atom_i.element in _NONPOLAR_ATOMS) and (atom_j.element in _NONPOLAR_ATOMS):
        interaction_type = "nonpolar"
    elif (atom_i.element not in _NONPOLAR_ATOMS) and (atom_j.element not in _NONPOLAR_ATOMS):
        interaction_type = "polar"
    else:
        interaction_type = "mixed"
    return interaction_type


def generate_residue_interaction_graph(residue_ids, atom_contacts):
    residue_interaction_graph = nx.Graph()
    residue_interaction_graph.add_nodes_from(residue_ids)
    assert len(residue_ids) == residue_interaction_graph.number_of_nodes()
    for atom_i, atom_j in atom_contacts:
        residue_i_id = atom_i.get_parent().get_full_id()
        residue_j_id = atom_j.get_parent().get_full_id()
        if residue_i_id != residue_j_id:
            assert residue_interaction_graph.has_node(residue_i_id)
            assert residue_interaction_graph.has_node(residue_j_id)
            interaction_type = _get_interaction_type(atom_i, atom_j)
            atom_distance = atom_i - atom_j
            if residue_interaction_graph.has_edge(residue_i_id, residue_j_id):
                residue_interaction_graph.edges[residue_i_id, residue_j_id]["interactions"][interaction_type].append(atom_distance)
            else:
                residue_interaction_graph.add_edge(residue_i_id, residue_j_id, interactions={"mixed": [], "nonpolar": [], "polar": []})
                residue_interaction_graph.edges[residue_i_id, residue_j_id]["interactions"][interaction_type].append(atom_distance)
    return residue_interaction_graph


def _main(structure_file, model_id, chain_ids, node_atom_selection, include_hetatms, cutoff):
    chains = parse_chains(structure_file, model_id, chain_ids)
    residue_ids, node_atoms = extract_node_atoms(chains, node_atom_selection, include_hetatms)
    atom_contacts = find_atom_contacts(node_atoms, cutoff)
    residue_interaction_graph = generate_residue_interaction_graph(residue_ids, atom_contacts)


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
