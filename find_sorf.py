import itertools
import sys, os
import subprocess

import numpy as np
import pandas as pd
import argparse
import sorf_lib as sf

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Find and Annotate sORFs')
    parser.add_argument('--genome', required=True,
                        help='name of the genome')
    parser.add_argument('--minimum_length', default=10,
                        help='minimum length per ORF')

    parser.add_argument('--out', default='.',
                        help='output path')
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    accession = sf.get_accession(genome)
    df = find_and_annotate_sorf(genome, accession,
                                minimum_length=args.minimum_length, out_dpath= args.out)
    df.to_feather(os.path.join(args.out, f'{genome}.sorf.feather.gz'))

