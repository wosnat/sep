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


module_fpath = os.path.dirname(__file__)

def run_orffinder(accession, minimum_length='1',
                  dna_code='11', atg_only=True, out_dpath=None):
    """
    run orffinder utility to find all ORFs. Return the name of the file
    Only works in linux (orffinder is linux only program)
    """
    
    cwd = os.getcwd()
    orffinder_exe = os.path.join(module_fpath, 'bin','ORFfinder')
    if out_dpath is None:
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
    
def list_orfs(accession, minimum_length='1000', out_dpath= None):
    """ given accession and minimum ORF length, run orf finder and return a data frame with list of orfs"""
    orf_fpath = run_orffinder(accession, minimum_length=minimum_length, out_dpath=out_dpath)
    return parse_orffinder_fasta(orf_fpath)

def load_genome(genome):
    """ load the genome of a given genome name. return dataframe """
    genomes_dpath = os.path.join(module_fpath, 'data',
                                 'detailed_Prochlorococcus_genome_annotations')
    fname = f'{genome}.txt'
    df = pd.read_csv(os.path.join(genomes_dpath, fname), sep='\t')
    df['genome'] = os.path.basename(os.path.splitext(fname)[0])
    df['left'] = df['start']
    df['right'] = df['stop']
    df.loc[df.strand == '-', 'left']  = df.loc[df.strand == '-', 'stop']
    df.loc[df.strand == '-', 'right'] = df.loc[df.strand == '-', 'start']

    return df


def add_gene_annotation_to_sorf_df(sorf_df, genome_df):
    """ add gene annotation to the SORF DF (only for genes that have the same ORF
    TODO: handle contigs
     """
    # TODO handle contigs
    merge_cols = ['contig_id', 'gene_id', 'location', 'type', 
                  'start', 'stop', 'strand',
                  'function', 'figfam',
                  'nucleotide_sequence', 'aa_sequence', 'genome', 'left', 'right']

    df_merge = sorf_df.merge(right=genome_df.loc[:, merge_cols], how='left',
                        left_on=['rast_left', 'rast_right', 'strand'],
                        right_on=['left', 'right', 'strand'],
                        suffixes=['', '_r']
                        )
    return df_merge


def _gene_to_loclist(x):
    data = {
        'location' :x.location,
        'strand' :x.strand,
    }
    d = pd.DataFrame(index=range(x.left, x.right + 1), data=data)
    return d

def _analyze_overlap(x, overlap_genes):
    """ analyze the overlap between sorf and known genes
    """
    
    if (overlap_genes.shape[0] == 0):
        result =  { 
            'overlap_type' : 'standalone'
        }
    elif not pd.isna(x['location']):
        result =  { 
            'overlap_type' : 'known'
        }
    elif (overlap_genes.shape[0] > 1):
        result =  { 
            'overlap_type' : 'overlap_many'
        }
    else:
        # exactly one overlap and not the same as the orf
        y = overlap_genes.squeeze()
        # y should be a series - one found
        is_same_strand = (x.strand == y.strand)
        is_out_of_frame = ((x.start+1 - y.start) %3) != 0
        is_inside = (x.rast_left >= y.left) and (x.rast_right <= y.right)
        is_upstream = (is_same_strand and 
                       ((x.strand == '+' and (x.start+1 < y.start)) or
                        (x.strand == '-' and (x.start+1 > y.start))))
        is_downstream = (is_same_strand and 
                       ((x.strand == '+' and (x.stop+1 > y.stop)) or
                        (x.strand == '-' and (x.stop+1 < y.stop))))
        
        out_of_frame_str = 'out_frame' if is_out_of_frame else None
        same_strand_str = 'as' if not is_same_strand else None
        inside_str = 'internal' if is_inside else None
        upstream_str = 'upstream' if is_upstream else None
        downstream_str = 'downstream' if is_downstream else None
        overlap_type = '_'.join(filter(None,[same_strand_str, 
                                             inside_str, downstream_str, upstream_str, 
                                             out_of_frame_str ]))
                        
        result = {
            'is_same_strand': is_same_strand,
            'is_out_of_frame': is_out_of_frame,
            'is_inside': is_inside,
            'is_upstream': is_upstream,
            'is_downstream': is_downstream,
            'overlap_type' : overlap_type,
        }
    return result

def _to_str(df, col):
    return ','.join(map(str, df[col].unique()))

def _find_overlaps(x, list_df, genome_df):
    subset_df = list_df.loc[(list_df.index >= x.rast_left) & (list_df.index <= x.rast_right)]
    overlap_ids = subset_df['location'].unique()
    overlap_genes = genome_df.loc[genome_df['location'].isin(overlap_ids)]
    result = {
        'overlap_location': _to_str(overlap_genes, 'location'), 
        'overlap_strand':   _to_str(overlap_genes, 'strand'), 
        'overlap_gene_type':     _to_str(overlap_genes, 'type'), 
        'overlap_count':    overlap_genes.shape[0]
    }
    result.update(_analyze_overlap(x, overlap_genes))
    return result
    

def add_overlapping_genes(sorf_df, genome_df):
    """ add cverlapping genes to df"""
    l = [_gene_to_loclist(x) for x in genome_df[['left', 'right', 'location', 'strand']].itertuples()]
    list_df = pd.concat(l)
    overlap_df = sorf_df.apply(
        lambda x : _find_overlaps(x, list_df, genome_df),
        axis=1, result_type='expand')
    return sorf_df.join(overlap_df, rsuffix='_')

def find_and_annotate_sorf(genome, accession, minimum_length='1000', out_dpath= None):
    """ find all ORFs, and add gene overlap info from the genome
    :return df with the data
    """
    df = list_orfs(accession, minimum_length, out_dpath)
    gdf = load_genome(genome)
    df = add_gene_annotation_to_sorf_df(df, gdf)
    df = add_overlapping_genes(df, gdf)
    return df


def get_accession(genome):
    metadf1 = pd.read_excel(os.path.join(module_fpath, 'data', 'metadata_pro_biller.xlsx'), sheet_name='Sheet1')
    metadf2 = pd.read_excel(os.path.join(module_fpath, 'data', 'metadata_pro_biller.xlsx'), sheet_name='Sheet2')
    accession_list = metadf1.loc[metadf1.Name == genome, 'NCBI accession'].unique()
    if len(accession_list) == 0:
        accession_list = metadf2.loc[metadf2.Strain == genome, 'NCBI accession'].unique()
    assert len(accession_list) == 1
    return accession_list[0]

if __name__ == '__main__':
    # tests
    _parse_orf_name('>lcl|ORF521_BX548174.1:17683:16679 unnamed protein product')
    _parse_orf_name('>lcl|ORF1_BX548174.1:2043:4382 unnamed protein product')

    df = list_orfs('BX548174')
    df.loc[df.len_aa3 != df.len_nn]
    gdf = load_genome('MED4')
    gdf.loc[gdf.start < 200]['aa_sequence'].unique() == df.loc[df.left < 200]['aaseq'].unique()
