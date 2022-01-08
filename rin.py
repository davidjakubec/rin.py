#!/usr/bin/env python3

from argparse import ArgumentParser
from os import path

from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.NeighborSearch import NeighborSearch
import networkx as nx
from Bio.PDB.SASA import ShrakeRupley

_NONPOLAR_ATOMS = {"C"}


def parse_model(structure_file, model_id=0):
    if structure_file.endswith(".cif"):
        parser = FastMMCIFParser(QUIET=True)
    structure_id = path.basename(structure_file)[:-4]
    model = parser.get_structure(structure_id, structure_file)[model_id]
    return model


def extract_node_atoms(model, node_atom_selection="CA", include_hetatms=False):
    residue_ids = []
    node_atoms = []
    for chain in model:
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


def _remove_hydrogens(model):
    for chain in model:
        for residue in chain:
            hydrogen_atom_ids = [atom.get_id() for atom in residue if atom.element == "H"]
            for hydrogen_atom_id in hydrogen_atom_ids:
                residue.detach_child(hydrogen_atom_id)
    assert [atom for atom in model.get_atoms() if atom.element == "H"] == []


def _remove_hetatms(model):
    for chain in model:
        hetatm_residue_ids = [residue.get_id() for residue in chain if residue.get_id()[0].strip()]
        for hetatm_residue_id in hetatm_residue_ids:
            chain.detach_child(hetatm_residue_id)
    assert [residue for residue in model.get_residues() if residue.get_id()[0].strip()] == []


def _get_atom_type(atom):
    if atom.element in _NONPOLAR_ATOMS:
        atom_type = "nonpolar"
    else:
        atom_type = "polar"
    return atom_type


def calculate_residue_sasas(model, include_hetatms=False):
    _remove_hydrogens(model)
    if not include_hetatms:
        _remove_hetatms(model)
    sr = ShrakeRupley()
    sr.compute(model)
    residue_sasas = {}
    for chain in model:
        residues = chain.get_list()
        for i, residue in enumerate(residues):
            residue_id = residue.get_full_id()
            residue_sasas[residue_id] = {"free": {"nonpolar": 0.0, "polar": 0.0}, "screened": {"nonpolar": 0.0, "polar": 0.0}}
            for atom in residue:
                residue_sasas[residue_id]["screened"][_get_atom_type(atom)] += atom.sasa
            residue_copy = residue.copy()
            if is_aa(residue):
                if i != 0:
                    try:
                        CX = residues[i - 1]["C"].copy()
                        CX.id = "CX_"
                        residue_copy.add(CX)
                    except KeyError:
                        pass
                if i != (len(residues) - 1):
                    try:
                        NX = residues[i + 1]["N"].copy()
                        NX.id = "NX_"
                        residue_copy.add(NX)
                    except KeyError:
                        pass
            sr.compute(residue_copy)
            atom_ids = [atom.get_id() for atom in residue]
            for atom_id in atom_ids:
                atom = residue_copy[atom_id]
                residue_sasas[residue_id]["free"][_get_atom_type(atom)] += atom.sasa
    return residue_sasas


def _main(structure_file, model_id, node_atom_selection, include_hetatms, cutoff):
    model = parse_model(structure_file, model_id)
    residue_ids, node_atoms = extract_node_atoms(model, node_atom_selection, include_hetatms)
    atom_contacts = find_atom_contacts(node_atoms, cutoff)
    residue_interaction_graph = generate_residue_interaction_graph(residue_ids, atom_contacts)
    residue_sasas = calculate_residue_sasas(model.copy(), include_hetatms)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("structure_file")
    parser.add_argument("--model_id", default=0, type=int)
    parser.add_argument("--node_atom_selection", default="CA")
    parser.add_argument("--include_hetatms", action="store_true")
    parser.add_argument("--cutoff", default=8.0, type=float)
    args = vars(parser.parse_args())
    _main(**args)
