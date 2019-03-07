import sys
import os
from subprocess import check_output
import math
import numpy as np
import pandas as pd
from collections import OrderedDict
import csv
import json

def find_organism_files(organism):
    ''' usage: organism, gff3, fa_dict, bowtie_index = find_organism_files(organism)'''
    if 'crypto' in organism:
        gff3 = '/home/jordan/GENOMES/CNA3_FINAL_CALLGENES_1_gobs.gff3'
        fa = '/home/jordan/GENOMES/H99_fa.json'
        bowtie_index = '/home/jordan/GENOMES/Crypto_for_gobs'
        organism = None
    elif 'pombe' in organism:
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        fa = '/home/jordan/GENOMES/POMBE/Sp_fasta_dict.json'
        bowtie_index = '/home/jordan/GENOMES/POMBE/Spombe'
        organism = 'pombe'
    elif 'cerev' in organism:
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        fa = '/home/jordan/GENOMES/S288C/S288C_fasta_dict.json'
        bowtie_index = '/home/jordan/GENOMES/S288C/S288C'
        organism = None
        
    with open(fa) as f:
        fa_dict = json.load(f)
        
    return organism, gff3, fa_dict, bowtie_index

def build_transcript_dict(gff3_file, organism=None):
    '''Creates a dictionary from gff3 file where keys are transcript names and values are [start, end, strand, chromosome, CDS_starts(list), CDS_ends(list]
    Note: change organism to 'pombe' if using pombe. Otherwise no need to indicate
    Additional note: See the tx_info function for retrieving gene information
    
    Parameters
    ----------
    gff3_file : str
              gff3 for your organism - tested on H99 C. neoformans, S288C S. cerevisiae, S. pombe and C. albicans
    organism : str, default `None`
              Change only if using pombe - then should be 'pombe'
    
    Returns
    -------
    tx_dict : File (comma separated) with annotation appended to each row.
              '''
    
    with open(gff3_file,"r") as gff3:
        transcript_dict = {}
        for line in gff3:
            columns = line.split('\t')
            
            if organism == 'pombe' and len(columns) > 1:
                chr_rom = columns[0]
                rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
                chrom = rom_lat[chr_rom]
                transcript_types = ['transcript','pseudogene','rRNA','snoRNA','tRNA','snRNA']
                if columns[2] in transcript_types:
                    if columns[8].split(':')[0].split('=')[1] == 'gene': continue
                    transcript = columns[8].split(';')[0].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], chrom, [], []]
                elif columns[2] == 'CDS':
                    transcript = columns[8].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
                        
            if len(columns) == 9 and organism is None:
                if columns[2] == "mRNA" or columns[2] == "snoRNA_gene" or columns[2] == "tRNA_gene":
                    transcript = columns[8]
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"): transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    #Transcript dictionary: keys are transcript, values are [start, end, strand, chromosome, CDS start, CDS end]
                    if transcript not in transcript_dict:
                        transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], columns[0], [], []]
                    else:
                        transcript_dict[transcript][0] = int(columns[3])
                        transcript_dict[transcript][1] = int(columns[4])
                        transcript_dict[transcript][2] = columns[6]
                        transcript_dict[transcript][3] = columns[0]
                elif columns[2] == "CDS":
                    transcript = columns[8].split("=")[1].split(".")[0].split(';')[0]
                    if 'mRNA' in transcript: transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    if transcript not in transcript_dict:
                        strand = columns[6]
                        chrom = columns[0]
                        transcript_dict[transcript] = [0,0,strand,chrom,[],[]]
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
    
    for tx in transcript_dict:
        if transcript_dict[tx][0] == 0:
            transcript_dict[tx][0] = transcript_dict[tx][4][0]
            transcript_dict[tx][1] = transcript_dict[tx][5][0]
    
    transcript_dict = OrderedDict(sorted(transcript_dict.items()))
    return transcript_dict


def seq_simple(chrom, start, end, strand, fasta_dict):
    if type(fasta_dict) == str:
        with open(fasta_dict, 'r') as f:
            fasta_dict = json.load(f)
    seq = fasta_dict[chrom][start:end+1]
    if strand == '-':
        seq = SP.reverse_complement(seq)
    return seq
            
def get_peak_sequence(input_file, fasta_file, gff3_file, window=1000):
    '''Makes a fasta file of peak sequences based on an input file.
    Input file columns - 1: transcript, 2: chromosome, 3: peak center
    Remember to save the input file as an MS-DOS CSV file if exporting from Excel
    Note: retrieves sequence
    
    Parameters
    ----------
    input_file : str
            CSV file - see above
    fasta_file : str
            .json dictionary of chromosome sequences or fasta file (.json will load faster)
    gff3_file : str
            gff3 file for your organism
    window : int, default 1000
            Size of sequence to retrieve (peak boundaries are window/2 on either side of peak summit)
            
    Outputs
    ------
    peak_fasta : fasta file with all peak sequences
    '''
    
    tx_dict = SP.build_transcript_dict(gff3_file)
    if type(fasta_file) == dict:
        fa_dict = fasta_file
    else:
        if fasta_file.endswith('json'):
            with open(fasta_file) as f:
                fa_dict = json.load(f)
        else:
            fa_dict = SP.make_fasta_dict(fasta_file)
    seq_list = []
    no_tx_n = 1
    with open(input_file,'r') as csv_file:
        f = csv.reader(csv_file, dialect=csv.excel)
        for row in f:
            tx_list = row[0].split(',')
            for tx in tx_list:
                tx = tx+'T0'
                if tx.startswith('3P'): tx = tx.split('3P')[1]
                
                chrom = row[1]
                if not chrom.startswith('chr'):
                    chrom = 'chr'+str(chrom)
                try:
                    center = int(row[2])
                    start = center-window/2
                    end = center+window/2
                    if tx in tx_dict:
                        strand = tx_dict[tx][2]
                    else:
                        print tx+" not in GFF3 file"
                        strand = '+'
                        tx = chrom+':'+str(center)
                    seq = seq_simple(chrom, start, end, strand, fa_dict)
                    seq_list.append((tx,seq))
                        
                except ValueError:
                    pass
    with open('{0}_peak_sequences.fa'.format(input_file.split('/')[-1].split('.')[0]),'w') as fout:
        for tx, seq in seq_list:
            fout.write('>'+tx+'\n')
            fout.write(seq+'\n')
    return seq_list


# This function reads my sorted-by-gene bedgraph format (output of SPTools.build_bedgraph_dict). Make sure the same
# transcript_dictionary is used for both function calls. In this case I have modified the dictionary to include 300 bp
# downstrame of the stop codon for each transcript

def read_bg_dict(bedgraph_dict_output, transcript_dict):
    bg_dict = {}
    count = 0
    dtype = [('coord', int), ('height', float)]
    with open(bedgraph_dict_output,'r') as f:
        n = 0
        for line in f:
            n += 1
            if len(line) > 1:
                
                #Read the transcript line
                if n%3 == 1:
                    tx = line.strip()
                    count += 1
                    #print count
                    if tx[-2] != 'T':
                        tx = tx+'T0'
                
                #Read the coordinate line
                elif n%3 == 2:
                    coords  = map(int, line.strip().split('\t'))
                
                #Read the values line
                elif n%3 == 0:
                    heights = map(float, line.strip().split('\t'))

                    all_coords = set(range(transcript_dict[tx][0],transcript_dict[tx][1]))
                    missing = all_coords.difference(coords)
                    coords = coords + list(missing)

                    #Fill in missing coordinates with zeros
                    zero_fill = [0]*len(missing)
                    heights = heights + zero_fill
                    
                    #Create a pandas series with all coordinates and sort so zeros are inserted appropriately
                    entry = pd.Series(heights, index=coords)
                    entry.sort_index(inplace=True)

                    bg_dict[tx] = entry
    return bg_dict

def seq_file_from_df(df, column_name, file_name):
    with open(file_name,'w') as fout:
        seq_list = df[column_name].tolist()
        for seq in seq_list:
            fout.write(seq+'\n')

def get_peak_sequence3(input_file, fasta_file, gff3_file, gene_list,window=1000):
    '''Makes a fasta file of peak sequences based on an input file.
    Input file columns - 1: transcript, 2: chromosome, 3: peak center
    Remember to save the input file as an MS-DOS CSV file if exporting from Excel
    Note: retrieves sequence
    
    Parameters
    ----------
    input_file : str
            CSV file - see above
    fasta_file : str
            .json dictionary of chromosome sequences or fasta file (.json will load faster)
    gff3_file : str
            gff3 file for your organism
    gene_list: str
             CSV file 
    window : int, default 1000
            Size of sequence to retrieve (peak boundaries are window/2 on either side of peak summit)
            
    Outputs
    ------
    peak_fasta : fasta file with all peak sequences
    '''

    tx_dict = SP.build_transcript_dict(gff3_file)
    if type(fasta_file) == dict:
        fa_dict = fasta_file
    else:
        if fasta_file.endswith('json'):
            with open(fasta_file) as f:
                fa_dict = json.load(f)
        else:
            fa_dict = SP.make_fasta_dict(fasta_file)
    seq_list = []
    no_tx_n = 1
    with open(input_file,'r') as csv_file:
        f = csv.reader(csv_file, dialect=csv.excel)
        for row in f:
            tx_list = row[0].split(',')
            for tx in tx_list:
                tx = tx+'T0'
                if tx.startswith('3P'): tx = tx.split('3P')[1]

                chrom = row[1]
                if not chrom.startswith('chr'):
                    chrom = 'chr'+str(chrom)
                try:
                    center = int(row[2])
                    start = center-window/2
                    end = center+window/2
                    if tx in tx_dict:
                        strand = tx_dict[tx][2]
                    else:
                        print tx+" not in GFF3 file"
                        strand = '+'
                        tx = chrom+':'+str(center)
                    seq = seq_simple(chrom, start, end, strand, fa_dict)
                    seq_list.append((tx,seq))

                except ValueError:
                    pass
    genes_of_interest=gene_list.split("\n")

    with open('{0}_peak_sequences.fa'.format(input_file.split('/')[-1].split('.')[0]),'w') as fout:
        for tx, seq in seq_list:
            if tx in genes_of_interest:
                fout.write('>'+tx+'\n')
                fout.write(seq+'\n')
    return seq_list
                     
