#!/usr/bin/env python
# coding: utf-8

import itertools
import sys, os
import subprocess

import numpy as np
import pandas as pd


# 
# ![](http://oregonstate.edu/instruct/bb450/fall14/stryer7/2/table_02_02.jpg)

# https://www.uniprot.org/uniprot/Q7V735
# 
# http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=TIGR03798
# 
# https://www.ebi.ac.uk/training/online/course/interpro-functional-and-structural-analysis-protei/sequence-searching/searching-interpro-batc
# 
# http://www.ebi.ac.uk/interpro/sequencesearch/iprscan5-S20190707-131508-0462-76111813-p1m
# 
# https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/iprscan5-S20190707-131508-0462-76111813-p1m/json


def run_orffinder(accession, minimum_length='1',
                  dna_code='11', atg_only=True):
    """
    run orffinder utility to find all ORFs. Return the name of the file
    Only works in linux (orffinder is linux only program)
    """
    
    cwd = os.getcwd()
    orffinder_exe = os.path.join(cwd, 'bin','ORFfinder')
    out_dpath = os.path.join(cwd,'orffinder_tmp')
    out_fpath = os.path.join(out_dpath, f'{accession}.orffinder.fa')
    
    
    params = [orffinder_exe, 
              '-id', accession, 
              '-g', dna_code,
              '-ml', minimum_length,
              '-s', '0', # ATG only start codon
              '-out', out_fpath,              
    ]
    
    try:
        os.makedirs(out_dpath, exist_ok=True)
        print(' '.join(params))
        retcode = subprocess.check_call(' '.join(params), shell=True)
        if retcode < 0:
            print("Child was terminated by signal", -retcode)
    except OSError as e:
        print("Execution failed:", e)
    except subprocess.CalledProcessError as e:
        print("Execution failed:", e)
    return out_fpath
    

def _parse_orf_name(name):
    """ parse line in orffinder output fasta (header of entry) """
    n=name.split()[0]
    n=n.split('|')[1]
    geneid, start, stop = n.split(':')
    contig = geneid.split('_', 2)[1]
    start = int(start)
    stop = int(stop)
    l = start
    r= stop
    strand = '+'
    if l >= r:
        strand = '-'
        l = stop
        r = start
    return {
        'orfid' : n,
        'contig' : contig,
        'left' : l,
        'right' : r,
        'start' : start,
        'stop' : stop,
        'strand' : strand,
    }



def parse_orffinder_fasta(orf_fpath):
    """ parse ORF finder fasta. return pd DataFrame with the ORFs """
    result = list()
    cur_item = None
    cur_aaseq = list()
    with open(orf_fpath, 'r') as fh:
        for l in fh:
            l = l.strip()
            if l.startswith('>'):
                # new item
                #>lcl|ORF1_BX548174.1:2043:4382 unnamed protein product
                if cur_item is not None:
                    cur_item['aaseq'] = ''.join(cur_aaseq)
                    result.append(cur_item)
                cur_aaseq = list()
                cur_item = _parse_orf_name(l)
            else:
                cur_aaseq.append(l)
    if cur_item is not None:
        cur_item['aaseq'] = ''.join(cur_aaseq)
        result.append(cur_item)

    cols = ['orfid', 'contig', 'left', 'right', 
            'start', 'stop', 'strand', 'aaseq', ]
    df = pd.DataFrame(result, columns=cols)
    for c in ['stop', 'left','right','start']:
        df[c] = pd.to_numeric(df[c])
    df['len_nn'] = (df.right + 1 - df.left) 
    df['len_aa'] = df['aaseq'].str.len()
    df['len_aa3'] = (df['len_aa']+1)*3
    df['rast_left'] = df['left'] + 1
    df['rast_right'] = df['right'] + 1
    return df
    
def list_orfs(accession, minimum_length='1000'):
    """ given accession and minimum ORF length, run orf finder and return a data frame with list of orfs"""
    orf_fpath = run_orffinder(accession, minimum_length=minimum_length)
    return parse_orffinder_fasta(orf_fpath)

def load_genome(genome):
    """ load the genome of a given genome name. return dataframe """
    genomes_dpath = os.path.join('.', 'data', 
                                 'detailed_Prochlorococcus_genome_annotations')
    fname = f'{genome}.txt'
    df = pd.read_csv(os.path.join(genomes_dpath, fname), sep='\t')
    df['genome'] = os.path.basename(os.path.splitext(fname)[0])
    df['left'] = gdf['start']
    df['right'] = gdf['stop']
    df.loc[gdf.strand == '-', 'left']  = df.loc[gdf.strand == '-', 'stop']
    df.loc[gdf.strand == '-', 'right'] = df.loc[gdf.strand == '-', 'start']

    return df


def add_gene_annotation_to_sorf_df(sorf_df, genome_df):
    """ add gene annotation to the SORF DF (only for genes that have the same ORF
    TODO: handle contigs
     """
    # TODO handle contigs
    df_merge = sorf_df.merge(right=genome_df, how='left',
         left_on=['rast_left', 'rast_right'],
         right_on=['left', 'right'],
         suffixes=['','_rast']
        )
    return df_merge


if __name__ == '__main__':
    # tests
    _parse_orf_name('>lcl|ORF521_BX548174.1:17683:16679 unnamed protein product')
    _parse_orf_name('>lcl|ORF1_BX548174.1:2043:4382 unnamed protein product')

    df = list_orfs('BX548174')
    df.loc[df.len_aa3 != df.len_nn]
    gdf = load_genome('MED4')
    gdf.loc[gdf.start < 200]['aa_sequence'].unique() == df.loc[df.left < 200]['aaseq'].unique()
