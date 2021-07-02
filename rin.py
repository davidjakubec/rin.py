#!/usr/bin/env python3

from argparse import ArgumentParser
from os import path

from Bio.PDB.MMCIFParser import FastMMCIFParser


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


def _main(structure_file, model_id, chain_ids):
    chains = parse_chains(structure_file, model_id, chain_ids)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("structure_file")
    parser.add_argument("--model_id", default=0, type=int)
    parser.add_argument("--chain_ids", nargs="+")
    args = vars(parser.parse_args())
    _main(**args)
