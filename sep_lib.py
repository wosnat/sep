#!/usr/bin/env python

import sys, os
import logging
from collections import defaultdict

import numpy as np
import pandas as pd

import pysam

class RnaAalignedRead:
    """
    class describing one read segment (grouping all reads that cover this segment)
    """
    def __init__(self):
        self.read_start = 0
        self.read_end = 0
        self.read_is_reverse = 0
        self.read_mapping_quality_sum = 0
        self.read_overflow_sum = 0
        self.read_utr_overflow_sum = 0
        self.read_count = 0

    def add_read(self, read :pysam.AlignedSegment , gene_start:int, gene_end:int, gene_is_reversed:bool):
        """
        add a new read to this segment
        :param read:  pysam.AlignedSegment
        :param gene_start: int start index beginning of reading frame). 1 based.
        :param gene_end: int end index (one ater the gene ends). 1 based.
        :param gene_is_reversed: bool is the gene on the reverse (antisense) strand
        :return: updates self
        """
        self.read_start = read.reference_start + 1
        self.read_end = read.reference_end
        self.read_is_reverse = read.is_reverse
        self.read_mapping_quality_sum += read.mapping_quality
        read_overflow, read_utr_overflow = \
            self._read_calculate_overlap(read, gene_start, gene_end, gene_is_reversed)
        self.read_overflow_sum += read_overflow
        self.read_utr_overflow_sum += read_utr_overflow
        self.read_count += 1

    def _read_calculate_overlap(self, read:pysam.AlignedSegment , gene_start:int, gene_end:int,
                                gene_is_reversed:bool):
        """
        calculate the overflow between the read and the gene. utr overflow: allow 100 bases before gene
        starts for utr (untranslated region)
        :param read:  pysam.AlignedSegment
        :param gene_start: int start index beginning of reading frame). 1 based.
        :param gene_end: int end index (one ater the gene ends). 1 based.
        :param gene_is_reversed: bool is the gene on the reverse (antisense) strand
        :return: tuple (int, int): (read overflow, read utr overflow)
        """
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

    @property
    def read_mapping_quality(self):
        return self.read_mapping_quality_sum / self.read_count

    @property
    def read_overflow(self):
        return self.read_overflow_sum / self.read_count

    @property
    def read_utr_overflow(self):
        return self.read_utr_overflow_sum / self.read_count

    @staticmethod
    def read_to_key(read: pysam.AlignedSegment):
        """
        :param read: samfile read
        :return: key for read dict
        """
        read_start = read.reference_start + 1
        read_end = read.reference_end
        read_is_reverse = read.is_reverse
        return (read_start, read_end, read_is_reverse)

    def to_record(self):
        """
        :return: current read as dict
        """
        return {
            'read_start': self.read_start,
            'read_end': self.read_end,
            'read_is_reverse': self.read_is_reverse,
            'read_mapping_quality': self.read_mapping_quality,
            'read_count': self.read_count,
            'read_overflow': self.read_overflow,
            'read_utr_overflow': self.read_utr_overflow,
        }


    def to_count(self, gene_is_reversed :bool, count_overflow: bool):
        """
        count number of reads covering this gene
        :param gene_is_reversed: bool is the gene on the reverse (antisense) strand
        :param count_overflow: whether to add partial cover reads to the count
        :return:  tuple: number_of_reads, number_of_antisense_reads
        """
        if not count_overflow and self.read_overflow:
            return 0, 0
        if self.read_is_reverse == gene_is_reversed:
            return self.read_count, 0
        else:
            return 0, self.read_count

    def to_cover(self, cover_df: pd.DataFrame,
                gene_start:int, gene_end:int, gene_is_reversed:bool):
        """
        compute cover per nucleotede in the gene
        :param cover_df: result of this function
        :param gene_start: int start index beginning of reading frame). 1 based.
        :param gene_end: int end index (one ater the gene ends). 1 based.
        :param gene_is_reversed: bool is the gene on the reverse (antisense) strand
        :return: nothing. return by reference of cover_df
        """
        if (gene_is_reversed == self.read_is_reverse):
            column = 'reads'
        else:
            column = 'reads_as'
        if self.read_overflow:
            column = f'overflow_{column}'
        start = max(self.read_start, gene_start)
        end = min(self.read_end, gene_end)
        cover_df.loc[range(start, end), column] = cover_df.loc[range(start, end), column] + self.read_count


class RnaReads:
    def __init__(self, samfile: pysam.AlignmentFile, contig: str,
                 gene_start:int, gene_end:int, gene_is_reversed:bool,
                 gene_type:str):
        """
        list reads (grouped by start,end, direction
        :param samfile: handle to pysam samfile
        :param contig: config name in the samfile
        :param gene_start: int start index beginning of reading frame). 1 based.
        :param gene_end: int end index (one ater the gene ends). 1 based.
        :param gene_is_reversed: bool is the gene on the reverse (antisense) strand
        :param gene_type: type (usually peg, rna)
        """
        self.reads_dict = defaultdict(RnaAalignedRead)
        self.gene_start = gene_start
        self.gene_end = gene_end
        if (gene_end is None) or (np.isnan(gene_end)):
            # last gene in the list
            self.gene_end = samfile.get_reference_length(contig)

        self.gene_type = gene_type
        self.gene_is_reversed = gene_is_reversed

        # pysam indexing is 0 based, annotation is 1 based
        samstart = max(0, gene_start - 1)
        samstop = self.gene_end

        if samstop <= samstart:
            return
        for read in samfile.fetch(contig=contig, start=samstart, stop=samstop):
            # group reads by the start index
            self.add_read(read)


    def add_read(self,read: pysam.AlignedSegment):
        """
        add read to as part of the constructor
        :param read: pysam read
        """
        logging.debug(read)
        key = RnaAalignedRead.read_to_key(read)
        self.reads_dict[key].add_read(read, self.gene_start, self.gene_end, self.gene_is_reversed)


    def get_record_list(self):
        """
        :return: a list of grouped reads (list of dicts)
        """
        return [r.to_record() for r in self.reads_dict.values()]

    def get_counts(self, count_overflow: bool):
        """

        :param count_overflow: should partially covered reads be counted
        :return: tuple (number of sense reads, number of antisense reads)
        """
        if not self.reads_dict:
            return  [0, 0]
        else:
            count_list = (r.to_count(self.gene_is_reversed, count_overflow) for r in self.reads_dict.values())
            counts = list(map(sum, zip(*count_list)))
            return counts

    def _cover_select_gene_column(self, is_intergenic:bool):
        if is_intergenic:
            gene_column = 'gene_inter'
        if not self.gene_is_reversed:
            gene_column = f'{self.gene_type}_sense'
        else:
            gene_column = f'{self.gene_type}_as'
        return gene_column

    def _cover_init_df(self, is_intergenic:bool):
        cover_df_columns = [
            'gene_inter',
            'peg_sense', 'peg_as',
            'rna_sense', 'rna_as',
            'reads', 'reads_as',
            'overflow_reads', 'overflow_reads_as',
        ]
        start = self.gene_start
        end = self.gene_end
        if self.gene_end <= self.gene_start:
            end = self.gene_start


        cover_df = pd.DataFrame(index=range(start, end),
                                columns=cover_df_columns)

        cover_df.loc[:, 'reads'] = 0
        cover_df.loc[:, 'reads_as'] = 0
        cover_df.loc[:, 'overflow_reads'] = 0
        cover_df.loc[:, 'overflow_reads_as'] = 0
        # mark the gene location
        gene_column = self._cover_select_gene_column(is_intergenic)
        cover_df.loc[:, gene_column] = 1
        return cover_df

    def get_cover_reads(self, is_intergenic:bool):
        """
        :param is_intergenic: is this intergenic (region between genes
        :return: df with number of covered reads per location in the gene
        """
        cover_df = self._cover_init_df(is_intergenic)
        for r in self.reads_dict.values():
            r.to_cover(cover_df, self.gene_start, self.gene_end, self.gene_is_reversed)
        return cover_df


############################################################3

def reads_per_gene(samfile: pysam.AlignmentFile, contig: str,
                 gene_start:int, gene_end:int, gene_is_reversed:bool,
                 gene_type:str):
    """

    :param samfile: handle to pysam samfile
    :param contig: config name in the samfile
    :param gene_start: int start index beginning of reading frame). 1 based.
    :param gene_end: int end index (one ater the gene ends). 1 based.
    :param gene_is_reversed: bool is the gene on the reverse (antisense) strand
    :param gene_type: type (usually peg, rna)
    :return: list of dicts per grouped read
    """
    read_dict = RnaReads(samfile, contig, gene_start, gene_end, gene_is_reversed, gene_type)
    return read_dict.get_record_list()


def count_reads(samfile: pysam.AlignmentFile, contig: str,
                 gene_start:int, gene_end:int, gene_is_reversed:bool,
                 gene_type:str, count_overflow):
    """
    :param samfile: handle to pysam samfile
    :param contig: config name in the samfile
    :param gene_start: int start index beginning of reading frame). 1 based.
    :param gene_end: int end index (one ater the gene ends). 1 based.
    :param gene_is_reversed: bool is the gene on the reverse (antisense) strand
    :param gene_type: type (usually peg, rna)
    :param count_overflow: should partially covered reads be counted
    :return: tuple (number of sense reads, number of antisense reads)
    """
    read_dict = RnaReads(samfile, contig, gene_start, gene_end, gene_is_reversed, gene_type)

    counts = read_dict.get_counts(count_overflow)
    return pd.Series({'reads': counts[0], 'as_reads': counts[1]})

def cover_reads(samfile: pysam.AlignmentFile, contig: str,
                 gene_start:int, gene_end:int, gene_is_reversed:bool,
                 gene_type:str, is_intergenic:bool):
    """
    :param samfile: handle to pysam samfile
    :param contig: config name in the samfile
    :param gene_start: int start index beginning of reading frame). 1 based.
    :param gene_end: int end index (one ater the gene ends). 1 based.
    :param gene_is_reversed: bool is the gene on the reverse (antisense) strand
    :param gene_type: type (usually peg, rna)
    :param is_intergenic: is this intergenic (region between genes
    :return: df with number of covered reads per location in the gene
    """
    read_dict = RnaReads(samfile, contig, gene_start, gene_end, gene_is_reversed, gene_type)
    cover_df = read_dict.get_cover_reads(is_intergenic)
    return cover_df



def add_useful_columns_to_annotation_df(ann_df:pd.DataFrame):
    """
    :param ann_df: dataframe based on rast annotation file
    :return: ann_df with additiona fields
    """
    ann_df['min_idx'] = ann_df[['start', 'stop']].min(axis=1)
    ann_df['max_idx'] = ann_df[['start', 'stop']].max(axis=1)
    ann_df['inter_stop_idx'] = ann_df.min_idx.shift(-1)
    ann_df['inter_length'] = ann_df.inter_stop_idx - ann_df.max_idx
    ann_df['gene_length'] = ann_df.max_idx - ann_df.min_idx
    ann_df['gene_is_reversed'] = (ann_df['strand'] == '-')
    ann_df['inter_length'] = ann_df.inter_length.fillna(0)
    return ann_df

def add_read_counts_to_annotation_df(ann_df:pd.DataFrame, samfile: pysam.AlignmentFile, contig:str):
    """
    :param ann_df: dataframe based on rast annotation file. assumes add_useful_columns_to_annotation_df()
    :param samfile: handle to pysam samfile
    :param contig: config name in the samfile
    :return: annotation df with additional fields :
        reads - reads in the direction of this gene
        reads_as - reads in the reverse direction
        reads_inter - reads in the forward direction covering the intergenic region after this gene
        reads_as_inter - reads in the reverse direction covering the intergenic region after this gene
    """
    count_func = lambda x: count_reads(samfile, contig,
            x['min_idx'], x['max_idx'], x['gene_is_reversed'], x['type'], count_overflow=True)
    inter_count_func = lambda x: count_reads(samfile, contig,
             x['max_idx'], x['inter_stop_idx'],
             gene_is_reversed=False, gene_type='intergenic', count_overflow=False)

    count_df = ann_df.apply(count_func, axis=1)
    inter_count_df = ann_df.apply(inter_count_func, axis=1)
    df = ann_df.join(count_df)
    df = df.join(inter_count_df, rsuffix='_inter')
    return df

def create_cover_df(ann_df:pd.DataFrame, samfile: pysam.AlignmentFile, contig:str):
    """
    :param ann_df: dataframe based on rast annotation file. assumes add_useful_columns_to_annotation_df()
    :param samfile: handle to pysam samfile
    :param contig: config name in the samfile
    :return: dataframe with row per location in the reference genome
    """
    cover_func = lambda x: cover_reads(samfile, contig,
                       x.iloc[0]['min_idx'], x.iloc[0]['max_idx'],
                       x.iloc[0]['gene_is_reversed'],
                       x.iloc[0]['type'], is_intergenic=False)
    inter_cover_func = lambda x: cover_reads(samfile, contig,
                       x.iloc[0]['max_idx'], x.iloc[0]['inter_stop_idx'],
                       gene_is_reversed=False,
                       gene_type='inter', is_intergenic=True)
    group_cols = ['contig_id', 'gene_id', 'type']
    cover_df = ann_df.groupby(group_cols).apply(cover_func)
    inter_count_df = ann_df.groupby(group_cols).apply(inter_cover_func)
    cover_df = pd.concat([cover_df, inter_cover_df])
    cover_df.reset_index(inplace=True)
    cover_df.rename(columns={'level_3': 'location'}, inplace=True)
    for c in ['gene_inter', 'peg_sense', 'peg_as', 'rna_sense', 'rna_as', 'reads', 'reads_as']:
        cover_df.loc[:,c] = pd.to_numeric(cover_df[c])


    return cover_df


if __name__ == '__main__':
    import pprint

    pd.set_option('display.max_columns', 500)
    ann_df = pd.read_csv('data/MIT9313.txt', sep='\t')
    ann_df = add_useful_columns_to_annotation_df(ann_df)
    #samfile = pysam.AlignmentFile("data/SRR3334788.sorted.bam", "rb")
    samfile = pysam.AlignmentFile("data/mit9313_SRR3334787.sorted.bam", "rb")
    lst = reads_per_gene(samfile, 'BX548175.1', 15, 20, False,'peg')
    print(len(lst))
    print(sum(r['read_count'] for r in lst))
    pprint.pprint(lst)

    #print(count_reads(samfile, 'BX548175.1', 15, 20, gene_is_reversed=True, count_overflow=True))
    print(count_reads(samfile, 'BX548175.1', 15, 20, gene_is_reversed=False, gene_type='peg' ,count_overflow=True))
    #print(count_reads(samfile, 'BX548175.1', 15, 20, gene_is_reversed=False, count_overflow=False))

    cover_df = cover_reads(samfile, 'BX548175.1', 15, 25, False, 'peg', False)
    print(cover_df)
