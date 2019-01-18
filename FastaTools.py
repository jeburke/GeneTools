import sys
import os
from subprocess import check_output
from subprocess import call
from subprocess import Popen
import math
import numpy as np
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneUtility
import SPTools as SP
from collections import OrderedDict
import csv


#Transcript dictionary: keys are transcript, values are [start, end, strand, chromosome, CDS start, CDS end]

def write_transcript_fasta(transcript_dict, fasta_dict, prefix='transcripts', sense=True, spliced=False):
    seq_dict = {}
    for transcript, values in transcript_dict.iteritems():
        start = values[0]
        end = values[1]
        strand = values[2]
        chrom = values[3]
        CDS_start_list = values[4]
        CDS_end_list = values[5]
        
        if spliced is False:
            seq = fasta_dict[chrom][start-1:end]
            if strand == '-':
                seq = SP.reverse_complement(seq)

        elif spliced is True:
            seq = ''
            for n in range(len(CDS_start_list)):
                if strand == '+':
                    seq = seq+fasta_dict[chrom][CDS_start_list[n]-1:CDS_end_list[n]]
                elif strand == '-':
                    new_seq = fasta_dict[chrom][CDS_start_list[n]-1:CDS_end_list[n]]
                    new_seq = SP.reverse_complement(new_seq)
                    seq = seq+new_seq
        
        if sense is False:
            seq = SP.reverse_complement(seq)
            
        seq_dict[transcript] = seq
        
    with open('{}.fa'.format(prefix), 'w') as fout:
        for transcript, seq in seq_dict.iteritems():
            fout.write('>'+transcript+'\n')
            fout.write(seq+'\n')
    
    return seq_dict
            
def write_intergenic_fasta(transcript_dict, fasta_dict, bps_us=0, bps_ds=0, all_intergenic=True, prefix='intergenic_transcripts'):
    seq_dict = {}
    if all_intergenic is False:
        for transcript, values in transcript_dict.iteritems():
            start = values[0]
            end = values[1]
            strand = values[2]
            chrom = values[3]
            
            if bps_us > 0:
                if strand == '+':
                    seq_us_sense = fasta_dict[chrom][start-bps_us:start]
                elif strand == '-':
                    seq_us_sense = fasta_dict[chrom][end:end+bps_us]
                    seq_us_sense = SP.reverse_complement(seq_us_sense)
                seq_us_antisense = SP.reverse_complement(seq_us_sense)
                seq_dict[transcript+'_us_sense'] = seq_us_sense
                seq_dict[transcript+'_us_antisense'] = seq_us_antisense
            
            if bps_ds > 0:
                if strand == '+':
                    seq_ds_sense = fasta_dict[chrom][end:bps_ds+end]
                elif strand == '-':
                    seq_ds_sense = fasta_dict[chrom][start-bps_ds:start]
                    seq_ds_sense = SP.reverse_complement(seq_ds_sense)
                seq_ds_antisense = SP.reverse_complement(seq_ds_sense)
                seq_dict[transcript+'_ds_sense'] = seq_ds_sense
                seq_dict[transcript+'_ds_antisense'] = seq_ds_antisense
    
    elif all_intergenic is True:
        chroms = fasta_dict.keys()
        for chrom in chroms:
            chrom_transcripts = dict((k, transcript_dict[k]) for k in transcript_dict if transcript_dict[k][3] == chrom)
            chr_txs_df = pd.DataFrame.from_dict(chrom_transcripts, orient='index')
            chr_txs_df.sort_values([0], inplace=True)
            sorted_transcripts = chr_txs_df.index.tolist()

            n = 0
            for n in range(len(sorted_transcripts)-1):
                transcript = sorted_transcripts[n]
                next_transcript = sorted_transcripts[n+1]
                transcript_end = chr_txs_df[1][transcript]
                next_start = chr_txs_df[0][next_transcript]
                if next_start > transcript_end:
                    seq_plus = fasta_dict[chrom][transcript_end:next_start]
                    seq_dict[transcript+'_'+next_transcript+'_plus'] = seq_plus
                    seq_dict[transcript+'_'+next_transcript+'_minus'] = SP.reverse_complement(seq_plus)
                else:
                    print 'Overlapping transcripts:'
                    print transcript
                    print next_transcript
                    
    with open('{}.fa'.format(prefix), 'w') as fout:
        for transcript, seq in seq_dict.iteritems():
            fout.write('>'+transcript+'\n')
            fout.write(seq+'\n')
    
    return seq_dict
        
def write_intron_fasta(transcript_dict, fasta_dict, prefix='introns', sense=True):
    seq_dict = {}
    for transcript, values in transcript_dict.iteritems():
        start = values[0]
        end = values[1]
        strand = values[2]
        chrom = values[3]
        CDS_start_list = values[4]
        CDS_end_list = values[5]

        for n in range(len(CDS_start_list)-1):
            if strand == '+':
                seq = fasta_dict[chrom][CDS_end_list[n]:CDS_start_list[n+1]-1]
            elif strand == '-':
                intron = len(CDS_start_list)-n-1
                seq = fasta_dict[chrom][CDS_end_list[intron]:CDS_start_list[intron-1]-1]
                seq = SP.reverse_complement(seq)
        
            if sense is False:
                seq = SP.reverse_complement(seq)
            
            seq_dict[transcript+'_'+str(n)] = seq
        
    with open('{}.fa'.format(prefix), 'w') as fout:
        for transcript, seq in seq_dict.iteritems():
            fout.write('>'+transcript+'\n')
            fout.write(seq+'\n')
    return seq_dict
        
    
def read_fasta(fasta_file):
    cds_dict = {}
    nts = ['A','T','C','G','N']
    with open(fasta_file, 'r') as fin:
        for line in fin:
            if line.startswith('>'):
                cds_dict[line.strip()] = ''
                cds = line.strip()
            elif line[0] in nts:
                cds_dict[cds] = cds_dict[cds]+line.strip()
    return cds_dict


def split_cds(cds_dict, end, number_nts):
    split_dict = {}
    for cds, seq in cds_dict.iteritems():
        if len(seq) <= number_nts:
            split_dict[cds] = seq
            continue
        if end == 'five_prime':
            split_dict[cds.split(' ')[0]+'_5'] = seq[:number_nts+1]
            split_dict[cds.split(' ')[0]+'_3'] = seq[number_nts+1:]
        elif end == 'three_prime':
            cds_len = len(seq)
            split_dict[cds.split(' ')[0]+'_5'] = seq[:(cds_len-number_nts)]
            split_dict[cds.split(' ')[0]+'_3'] = seq[(cds_len-number_nts):]
    return split_dict

def write_new_fasta(split_dict, fasta_file):
    new_name = fasta_file.split('.')[0]+'_split.fa'
    with open(new_name, 'w') as fout:
        for cds, seq in split_dict.iteritems():
            fout.write(cds+'\n')
            fout.write(seq+'\n')