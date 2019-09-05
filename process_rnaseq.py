import sep_lib as sep
import pysam
import sys, os
import argparse
import logging

import numpy as np
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--ann', "--annotation_csv", required=True, metavar='ANNFILE',
                        help="name of annotation csv file")
    parser.add_argument('--aln', "--alignment_bam", required=True, metavar='ALNFILE',
                        help="alignment file. Should be sorted and indexed BAM file")
    parser.add_argument('--o', "--output_dir", metavar='OUTDIR', default='.',
                        help="output directory")
    parser.add_argument('--d', "--debug", action='store_true',
                        help="debug mode")
    parser.add_argument('--head', action='store_true',
                        help="run only on annotation head")


    args = parser.parse_args()
    loglevel = logging.INFO
    if args.d:
        loglevel = logging.DEBUG
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=loglevel,
        datefmt='%Y-%m-%d %H:%M:%S')

    ann_fpath = args.ann
    aln_fpath = args.aln
    out_dpath = args.o


    if not os.path.isfile(ann_fpath):
        logging.error(f'annotation file {ann_fpath} not found')
        exit(1)
    if not os.path.isfile(aln_fpath):
        logging.error(f'alignment file {aln_fpath} not found')
        exit(1)

    os.makedirs(out_dpath, exist_ok=True)
    read_fpath = os.path.basename(aln_fpath)
    read_name = read_fpath.split('.')[0]

    counts_csv_fpath = os.path.join(out_dpath, f'{read_name}_counts.csv.gz')
    cover_csv_fpath = os.path.join(out_dpath, f'{read_name}_cover.csv.gz')

    logging.info(f"opening annotation file {ann_fpath}")
    ann_df = pd.read_csv(ann_fpath, sep='\t')
    ann_df = sep.add_useful_columns_to_annotation_df(ann_df)

    if args.head:
        ann_df = ann_df.head()
    logging.info(f"opening alignment file {aln_fpath}")
    samfile = pysam.AlignmentFile(aln_fpath, "rb")
    contig = samfile.references[0]

    logging.info(f"computing read counts")
    count_df = sep.add_read_counts_to_annotation_df(ann_df, samfile, contig)
    logging.info(f"done computing read counts")

    logging.info(f"saving to f{counts_csv_fpath}")
    count_df.to_csv(counts_csv_fpath)
    logging.info(f"done saving to f{counts_csv_fpath}")

    logging.info(f"computing read cover")
    cover_df = sep.create_cover_df(ann_df, samfile, contig)
    logging.info(f"done computing read cover")

    logging.info(f"saving to f{cover_csv_fpath}")
    cover_df.to_csv(cover_csv_fpath)
    logging.info(f"done saving to f{cover_csv_fpath}")

