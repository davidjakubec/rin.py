#!/usr/bin/env python3

from argparse import ArgumentParser
# from os import path


def _main(structure_file):
    pass


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("structure_file")
    args = vars(parser.parse_args())
    _main(**args)
