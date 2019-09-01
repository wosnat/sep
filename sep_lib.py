import itertools
import sys, os

from scipy.special import comb
from scipy import stats
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from sklearn.preprocessing import StandardScaler
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from sklearn.metrics import classification_report, accuracy_score
#import sklearn.metrics as metrics
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import auc

from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from scipy import stats

import pysam

def _read_calculate_overlap(read,
                            gene_start, gene_end,
                            gene_is_reversed):
    start = max(0, gene_start - 1)
    stop = gene_end
    # utr_start / utr_stop
    utr_start = start
    utr_stop = stop
    if gene_is_reversed:
        utr_stop = utr_stop + 100
    else:
        utr_start = max(0, utr_start - 100)
    overlap = read.get_overlap(start, stop)
    read_overflow = (read.reference_end - read.reference_start) - overlap
    utr_overlap = read.get_overlap(utr_start, utr_stop)
    read_utr_overflow = (read.reference_end - read.reference_start) - utr_overlap

    return read_overflow, read_utr_overflow


def _read_add_to_dict(read_dict, read,
                      gene_start, gene_end,
                      gene_is_reversed):
    read_start = read.reference_start + 1
    read_end = read.reference_end
    read_is_reverse = read.is_reverse
    read_mapping_quality = read.mapping_quality
    read_overflow, read_utr_overflow = _read_calculate_overlap(read, gene_start, gene_end, gene_is_reversed)

    key = (read_start, read_end, read_is_reverse)
    read_dict.setdefault(key, [0, 0, 0, 0])
    read_dict[key][0] += read_mapping_quality
    read_dict[key][1] += 1
    read_dict[key][2] += read_overflow
    read_dict[key][3] += read_utr_overflow


def _read_decode_key_value(key, value):
    (read_start, read_end, read_is_reverse) = key
    (read_mapping_quality_sum, read_count,
     read_overflow_sum, read_utr_overflow_sum) = value

    read_mapping_quality = read_mapping_quality_sum / read_count
    read_overflow = read_overflow_sum / read_count
    read_utr_overflow = read_utr_overflow_sum / read_count
    return (read_start, read_end, read_is_reverse, read_count,
            read_mapping_quality, read_overflow, read_utr_overflow)


def _read_dict_item_to_record(key, value):
    (read_start, read_end, read_is_reverse, read_count,
     read_mapping_quality, read_overflow, read_utr_overflow
     ) = _read_decode_key_value(key, value)

    return {
        'read_start': read_start,
        'read_end': read_end,
        'read_is_reverse': read_is_reverse,
        'read_mapping_quality': read_mapping_quality,
        'read_count': read_count,
        'read_overflow': read_overflow,
        'read_utr_overflow': read_utr_overflow,
    }


def _read_dict_item_to_count(key, value, gene_is_reversed, count_overflow):
    """ return tuple:
    number_of_reads, number_of_antisense_reads

    """
    (read_start, read_end, read_is_reverse, read_count,
     read_mapping_quality, read_overflow, read_utr_overflow
     ) = _read_decode_key_value(key, value)

    if not count_overflow and read_overflow:
        return 0, 0
    if read_is_reverse == gene_is_reversed:
        return read_count, 0
    else:
        return 0, read_count


# TODO: UTR - 100 bases before gene start
# TODO: overlap for partial overlap genes

def _reads_to_dict(samfile, contig,
                   gene_start, gene_end,
                   gene_is_reversed):
    # pysam indexing is 0 based, annotation is 1 based
    start = max(0, gene_start - 1)
    stop = gene_end

    if (stop is None) or (np.isnan(stop)):
        # last gene in the list
        return dict()

    if stop <= start:
        return dict()
    read_dict = dict()
    for read in samfile.fetch(contig=contig, start=start, stop=stop):
        # group reads by the start index
        print(read)
        _read_add_to_dict(read_dict, read, gene_start, gene_end, gene_is_reversed)
    return read_dict


def reads_per_gene(samfile, contig,
                   gene_start, gene_end,
                   gene_is_reversed):
    read_dict = _reads_to_dict(samfile, contig,
                               gene_start, gene_end,
                               gene_is_reversed)

    reads_list = [_read_dict_item_to_record(k, v) for k, v in read_dict.items()]

    return reads_list


def count_reads(samfile, contig,
                gene_start, gene_end,
                gene_is_reversed, count_overflow):
    read_dict = _reads_to_dict(samfile, contig,
                               gene_start, gene_end,
                               gene_is_reversed)

    count_list = [_read_dict_item_to_count(k, v, gene_is_reversed, count_overflow)
                  for k, v in read_dict.items()]
    if count_list == []:
        counts = [0, 0]
    else:
        counts = list(map(sum, zip(*count_list)))
    return pd.Series({'reads': counts[0], 'as_reads': counts[1]})


def _read_dict_item_to_cover(cover_df,
                             gene_start, gene_end, gene_is_reversed,
                             key, value):
    (read_start, read_end, read_is_reverse, read_count,
     read_mapping_quality, read_overflow, read_utr_overflow
     ) = _read_decode_key_value(key, value)
    if (gene_is_reversed == read_is_reverse):
        column = 'reads'
    else:
        column = 'reads_as'
    if read_overflow:
        column = f'overflow_{column}'
    start = max(read_start, gene_start)
    end = min(read_end, gene_end)
    cover_df.loc[range(start, end), column] = cover_df.loc[range(start, end), column] + read_count


def _select_gene_column(gene_type, gene_is_reversed, is_intergenic):
    if not gene_is_reversed:
        gene_column = f'{gene_type}_sense'
    else:
        gene_column = f'{gene_type}_as'
    if is_intergenic:
        gene_column = 'gene_inter'
    return gene_column


def cover_reads(samfile, contig,
                gene_start, gene_end,
                gene_is_reversed, gene_type, is_intergenic):
    read_dict = _reads_to_dict(samfile, contig,
                               gene_start, gene_end,
                               gene_is_reversed)

    cover_df_columns = [
        'gene_inter',
        'peg_sense', 'peg_as',
        'rna_sense', 'rna_as',
        'reads', 'reads_as',
        'overflow_reads', 'overflow_reads_as',
    ]
    cover_df = pd.DataFrame(index=range(gene_start, gene_end),
                            columns=cover_df_columns)

    cover_df.loc[:, 'reads'] = 0
    cover_df.loc[:, 'reads_as'] = 0
    cover_df.loc[:, 'overflow_reads'] = 0
    cover_df.loc[:, 'overflow_reads_as'] = 0
    # mark the gene location
    gene_column = _select_gene_column(gene_type, gene_is_reversed, is_intergenic)
    cover_df.loc[:, gene_column] = 1

    for k, v in read_dict.items():
        _read_dict_item_to_cover(cover_df,
                                 gene_start, gene_end, gene_is_reversed,
                                 k, v)

    return cover_df


def add_useful_columns_to_annotation_df(ann_df):
    ann_df['min_idx'] = ann_df[['start', 'stop']].min(axis=1)
    ann_df['max_idx'] = ann_df[['start', 'stop']].max(axis=1)
    ann_df['inter_stop_idx'] = ann_df.min_idx.shift(-1)
    ann_df['inter_length'] = ann_df.inter_stop_idx - ann_df.max_idx
    ann_df['gene_length'] = ann_df.max_idx - ann_df.min_idx
    ann_df['gene_is_reversed'] = (ann_df['strand'] == '-')
    ann_df['inter_length'] = ann_df.inter_length.fillna(0)
    return ann_df

def add_read_counts_to_annotation_df(ann_df, samfile, contig):

    count_func = lambda x: count_reads(samfile, contig,
            x['min_idx'], x['max_idx'], x['gene_is_reversed'], count_overflow=True)
    count_df = ann_df.apply(count_func, axis=1)
    inter_count_func = lambda x: count_reads(samfile, contig,
             x['max_idx'], x['inter_stop_idx'],
             gene_is_reversed=False, count_overflow=False)
    inter_count_df = ann_df.apply(inter_count_func, axis=1)
    df = ann_df.join(count_df)
    df = df.join(inter_count_df, rsuffix='_inter')
    return df



if __name__ == '__main__':
    import pprint

    pd.set_option('display.max_columns', 500)
    ann_df = pd.read_csv('MIT9313.txt', sep='\t')
    ann_df = add_useful_columns_to_annotation_df(ann_df)
    samfile = pysam.AlignmentFile("SRR3334788.sorted.bam", "rb")
    #lst = reads_per_gene(samfile, 'BX548175.1', 15, 20, False)
    #print(len(lst))
    #print(sum(r['read_count'] for r in lst))
    #pprint.pprint(lst)

    #print(count_reads(samfile, 'BX548175.1', 15, 20, gene_is_reversed=True, count_overflow=True))
    #print(count_reads(samfile, 'BX548175.1', 15, 20, gene_is_reversed=False, count_overflow=True))
    #print(count_reads(samfile, 'BX548175.1', 15, 20, gene_is_reversed=False, count_overflow=False))

    cover_df = cover_reads(samfile, 'BX548175.1', 15, 25, False, 'peg', False)
    print(cover_df)
