import sys
import os
from subprocess import check_output
import math
import numpy as np
import pandas as pd
script_path = os.path.dirname(os.path.realpath(__file__)).split('GeneTools')[0]
sys.path.append(script_path)
import GeneTools as GT
from collections import OrderedDict
import csv
import pysam
from matplotlib import pyplot as plt

def count_reads_in_window(bam, chrom, start, end, strand):
    '''Counts reads in a given window on one strand - assumes reads are from cDNA
    
    Parameters
    ----------
    bam : str
         Bam file generated by Bowtie or STAR
    chrom : str
         chromosome name (needs to match references in bam)
    start : int
         start of the window
    end : int
         end of the window (must be larger than start)
    strand : str
         "+" or "-"
    
    Returns
    ------
    read_count : int
         number of reads aligned to window'''
    
    if type(bam) == str:
        bam = pysam.Samfile(bam)
    bam_iter = bam.fetch(chrom, start, end)
    read_count = 0
    for read in bam_iter:
        if strand == "+":
            if read.is_reverse:
                read_count += 1
        elif strand == "-":
            if not read.is_reverse:
                read_count +=1
    return read_count

def tx_info(tx, tx_dict):
    '''Retrieves information on a given transcript from tx_dict generated by GeneTools.build_transcript_dict
    
    Parameters
    ----------
    tx : str
              transcript name - should end in something like 'T0' or '.1' in S. pombe
    tx_dict : OrderedDict
              transcript dictionary generated by GeneTools.build_transcript_dict
    
    Returns
    -------
    start : int
            start of transcript (always < end)
    end : int
            end of transcript (always > start)
    chrom : str
            chromosome name
    strand : str
            "+" or "-"
    CDS_start : int
            Start of the coding sequence (always < CDS_end)
    CDS_end : int
            End of the coding sequence (always > CDS_start)
    exons : list of tuples
            List of exon (start,end) - numbers are ordered lowest to highest, not based on transcript direction.
    '''
    
    chrom = tx_dict[tx][3]
    strand = tx_dict[tx][2]
    start = tx_dict[tx][0]
    end = tx_dict[tx][1]
    exons = None
    if len(tx_dict[tx]) > 4 and len(tx_dict[tx][4]) > 0:
        CDS_start = min(tx_dict[tx][4])
        CDS_end = max(tx_dict[tx][5])
        exons = zip(tx_dict[tx][4],tx_dict[tx][5])
    else:
        CDS_start = start
        CDS_end = end
        exons = [(start,end)]
        
    return start, end, chrom, strand, CDS_start, CDS_end, exons

def list_splice_sites(gff3_file, chromosome="All", gene_list=None, organism=None):
    '''Function to build dictionary of splice sites in each transcript.
    Parameters
    ----------
    gff3_file : str
                location and name of gff3 file
    chromosome : str, default "All"
               If not "All", will only list transcripts on the selected chromosome
    gene_list : list, default ``None``
               List of genes to limit the dictionary. Useful for optimizing performance of downstream functions.
    organism : str, default ``None``
                only necessary for S. pombe - then use 'pombe'
                
    Returns
    --------
    splice_site_dict : dict
                    Dictionary where transcript names are keys and values are list as follows:
                    [[intron starts], [intron stops]]
    intron flag : bool
                    True if introns rather than exons are defined in the gff3 file.
                    '''
        
    transcript_dict = build_transcript_dict(gff3_file, organism=organism)
    
    splice_site_dict = {}
    n = 1
    with open(gff3_file,"r") as fin:
        for line in fin:
            columns = re.split(r'\t+', line.strip())

            if organism == 'pombe' and len(columns)>1:
                intron_flag = False
                if columns[2] == 'exon':
                    chr_rom = columns[0]
                    rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
                    chrom = rom_lat[chr_rom]
                    transcript = columns[8].split(':')[0].split('=')[1]

                    if transcript not in splice_site_dict:
                        splice_site_dict[transcript] = [[],[],chrom]
                    if columns[6] == "+":
                        splice_site_dict[transcript][0].append(int(columns[4])-1)
                        splice_site_dict[transcript][1].append(int(columns[3])-2)
                    elif columns[6] == "-":
                        splice_site_dict[transcript][0].append(int(columns[3])-1)
                        splice_site_dict[transcript][1].append(int(columns[4]))

            if len(columns) > 1 and organism is None:
                if columns[2] == "mRNA" or columns[2] == "snoRNA_gene" or columns[2] == "tRNA_gene":
                    transcript = columns[8].strip()
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"):
                        transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T':
                        transcript = transcript+'T0'
                    if transcript not in splice_site_dict:
                        splice_site_dict[transcript] = [[],[],columns[0]]

                elif columns[2] == "exon":
                    intron_flag=False
                    transcript = columns[8].strip()
                    transcript = transcript[-12:]
                    if transcript[-2] != 'T':
                        transcript = transcript+'T0'
                    if gene_list is None:
                        if columns[6] == "+":
                            splice_site_dict[transcript][0].append(int(columns[4])-1)
                            splice_site_dict[transcript][1].append(int(columns[3])-2)
                        if columns[6] == "-":
                            splice_site_dict[transcript][0].append(int(columns[3])-1)
                            splice_site_dict[transcript][1].append(int(columns[4]))
                    else:
                        if transcript in gene_list:
                            if columns[6] == "+":
                                splice_site_dict[transcript][0].append(int(columns[4])-1)
                                splice_site_dict[transcript][1].append(int(columns[3])-2)
                            if columns[6] == "-":
                                splice_site_dict[transcript][0].append(int(columns[3])-1)
                                splice_site_dict[transcript][1].append(int(columns[4]))

                #For organisms where introns are annotated instead of exons (e.g. S. cerevisiae)
                elif "intron" in columns[2]: 
                    intron_flag=True
                    transcript = columns[8].strip()
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"):
                        transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T':
                        transcript = transcript+'T0'
                    if transcript not in splice_site_dict:
                        splice_site_dict[transcript] = [[],[],columns[0]]
                    if gene_list is None:
                        if columns[6] == "+":
                            splice_site_dict[transcript][0].append(int(columns[3])-1)
                            splice_site_dict[transcript][1].append(int(columns[4]))
                        elif columns[6] == "-":
                            splice_site_dict[transcript][0].append(int(columns[4]))
                            splice_site_dict[transcript][1].append(int(columns[3])-1)
                    else:
                        if transcript in gene_list:
                            if columns[6] == "+":
                                splice_site_dict[transcript][0].append(int(columns[3])-2)
                                splice_site_dict[transcript][1].append(int(columns[4])-1)
                            elif columns[6] == "-":
                                splice_site_dict[transcript][0].append(int(columns[4]))
                                splice_site_dict[transcript][1].append(int(columns[3])-1)
    
    
    #Trim to entries in gene list
    if gene_list is not None:
        splice_site_dict = {transcript: splice_site_dict[transcript] for transcript in gene_list}
        
    if chromosome != "All":
        splice_site_dict = dict([(transcript, coords) for transcript, coords in splice_site_dict.iteritems() if coords[2] == chromosome])
    
    if intron_flag is False:
        for transcript, sites in splice_site_dict.iteritems():
            sites5 = sorted(sites[0])
            sites3 = sorted(sites[1])
            if transcript_dict[transcript][2] == "+":
                sites5 = sites5[:-1]
                sites3 = sites3[1:]
            elif transcript_dict[transcript][2] == "-":
                sites5 = sites5[1:]
                sites3 = sites3[:-1]
            splice_site_dict[transcript] = [sites5, sites3]
        
    return (splice_site_dict, intron_flag)    

def collapse_ss_dict(splice_site_dict):
    '''Function to list introns by gene rather than isoform - removes redundant introns
    
    Parameters
    ----------
    splice_site_dict : dict
                splice site dictionary generated by list_splice_sites
                
    Returns
    --------
    splice_by_gene : dict
                    Dictionary where gene names (no "T0" or ".1" at the end) are keys and values are set as follows:
                    ([intron starts], [intron stops])
                    '''
    ss_by_gene = {}
    transcripts = splice_site_dict.keys()
    transcripts = [x[:-2] for x in transcripts]
    genes = list(set(transcripts))
    for gene in genes:
        for transcript, sites in splice_site_dict.iteritems():
            if gene in transcript:
                if gene not in ss_by_gene:
                    ss_by_gene[gene] = set()
                ss_by_gene[gene].update(set(zip(sites[0],sites[1])))
    return ss_by_gene


def count_aligned_reads(bam_file):
    '''Counts aligned reads in bam file reads using Samtools
    
    Parameters
    ----------
    bam_file : str
            bam file from Bowtie or STAR
    
    Returns
    ------
    total : float
         Million aligned reads'''
    
    total = check_output(['samtools','view','-F 0x04','-c',bam_file]).strip()
    total = float(total)/1000000.
    return total

### Builds a pandas series from a bam file for a set of coordinates. Index is genome coordinate and value is number of reads
def generate_read_series(bam_iterator, chrom, start, end, strand, baseline=0):
    '''Generates a pandas series that is effectively a bedgraph for a given genome region - index is location on chromosome and 
    values are the number of reads at that position. This object can be used for plotting and is easily normalized, smoothed, etc.
    
    Parameters
    ----------
    bam_iterator : pysam.libcsamfile.Samfile or pysam.libcalignmentfile.IteratorRowRegion
         Bam file opened with pysam.Samfile
         For improved speed - perform fetch on bam file to generate a pysam.libcalignmentfile.IteratorRowRegion object
    chrom : str
         chromosome name (needs to match references in bam)
    start : int
         start of the window
    end : int
         end of the window (must be larger than start)
    strand : str
         "+" or "-"
    baseline : int or float, default 0
         fill value for the baseline if no reads at a given position
    
    Returns
    ------
    s : pandas.core.series.Series
         read counts at each position in genome (see above)'''
    
    s = pd.Series(baseline, index=range(start-25, end+25))
    for read in bam_iterator:
        if read.is_reverse and strand == '+':
            pos = read.reference_start
            if pos not in s.index:
                s[pos] = 0
            s[pos] += 1
        elif not read.is_reverse and strand == '-':
            pos = read.reference_start
            if pos not in s.index:
                s[pos] = 0
            s[pos] = s[pos]+1
    s = s.dropna()
    #s = s[s > 0]
    s = s.sort_index()
    return s

### Build a dictionary of read series based on a bam file and transcript dictionary
def map_all_transcripts(gff3, bam_file):
    organism=None
    if 'pombe' in gff3:
        organism = 'pombe'
    tx_dict = GT.build_transcript_dict(gff3, organism=organism)
    
    bam = pysam.Samfile(bam_file)
    series_dict = {}
    for tx in tx_dict:
        start, end, chrom, strand, CDS_start, CDS_end, exons = tx_info(tx, tx_dict)
        series_dict[tx] = generate_read_series(bam, chrom, start, end, strand)
    return series_dict

def plot_transcripts(series_dict, tx_list):
    if type(tx_list) != list:
        tx_list = [tx_list]
    
    for tx in tx_list:
        fig = plt.figure(figsize=(12,6))
        ax = fig.add_subplot(111)
        ax.bar(series_dict[tx].index, series_dict[tx], width, color='darkslateblue')

        plt.show()
        plt.clf()
    
def count_PE_reads(open_bam, chrom, start, end, strand, both_strands=False, count_junctions=False):
    if type(open_bam) != pysam.libcsamfile.Samfile:
        open_bam = pysam.Samfile(open_bam)
    
    iterator = open_bam.fetch(chrom, start-150, end+150)
    
    count = 0
    for read in iterator:
        if both_strands is False:
            # For reads that start or end in the intron
            if read.reference_start in range(start, end) or read.reference_end in range(start, end):
                if not read.is_reverse and strand == '+': count += 1
                elif read.is_reverse and strand == '-': count += 1

            # For reads that span the intron
            elif read.reference_start <= start and read.reference_end >= end:
                intron = False
                if count_junctions is False:
                    # Check to see if the read contains a junction
                    if len(read.cigartuples) == 3 and read.cigartuples[1][0] == 3:
                        intron = True

                if intron is False:
                    if not read.is_reverse and strand == '+': count += 1
                    elif read.is_reverse and strand == '-': count += 1
                            
        else:
            # Don't need to worry about strand or read1 vs. read2, otherwise same as above
            if read.reference_start in range(start, end) or read.reference_end in range(start, end):
                count += 1
            elif read.reference_start <= start and read.reference_end >= end:
                intron = False
                if count_junctions is False:
                    # Check to see if the read contains a junction
                    if len(read.cigartuples) == 3 and read.cigartuples[1][0] == 3:
                        intron = True
                if intron is False:
                    count += 1
                    
    return count

def PE_intron_retention_from_annotation(bam_list, organism, both_strands=False, count_junctions=False):
    if 'crypto' in organism.lower():
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
        organism=None
    elif 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        organism='pombe'
    elif 'cerev' in organism.lower():
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        organism=None
        
    tx_dict = GT.build_transcript_dict(gff3, organism=organism)
    ss_dict, flag = list_splice_sites(gff3, organism=organism)
    ss_dict = collapse_ss_dict(ss_dict)
    #ss_dict = {k:v for k, v in ss_dict.items() if k in ss_dict.keys()[:100]}
    
    open_bams = {}
    data_dict = {}
    for bam in bam_list:
        open_bams[bam] = pysam.Samfile(bam)
        name = bam.split('/')[-1].split('_sorted.bam')[0]
        print name
        data_dict[bam] = {bam:name, 'reads in transcript':[], 'reads in intron':[]}
        
    column_dict = {'transcript':[],'intron start':[],'intron end':[],'transcript size':[],'chromosome':[],'strand':[]}
    
    for tx, splice_sites in ss_dict.iteritems():
        print tx
        # Get information for overall transcript
        if organism == 'pombe': iso = tx+'.1'
        else: iso = tx+'T0'
        start, end, chrom, strand, CDS_start, CDS_end, exons = GT.tx_info(iso, tx_dict)

        tx_counts = {}
        for bam, open_bam in open_bams.iteritems():
            # Count reads in transcript
            tx_counts[bam] = count_PE_reads(open_bam, chrom, start, end, strand, both_strands=both_strands, count_junctions=True)
        
        # Iterate over all annotated introns in transcript
        for five, three in splice_sites:
            column_dict['transcript'].append(tx)
            column_dict['intron start'].append(five)
            column_dict['intron end'].append(three)
            column_dict['transcript size'].append(end-start)
            column_dict['chromosome'].append(chrom)
            column_dict['strand'].append(strand)
            
            for bam, open_bam in open_bams.iteritems():
                if strand == '+':
                    intron_counts = count_PE_reads(open_bam, chrom, five, three, strand, both_strands=both_strands, count_junctions=count_junctions)
                if strand == '-':
                    intron_counts = count_PE_reads(open_bam, chrom, three, five, strand, both_strands=both_strands, count_junctions=count_junctions)
                
                data_dict[bam]['reads in transcript'].append(tx_counts[bam])
                data_dict[bam]['reads in intron'].append(intron_counts)
                
    df = pd.DataFrame(columns=column_dict.keys(), index=range(len(column_dict['transcript'])))
    for col, info in column_dict.iteritems():
        df[col] = info
    df['intron size'] = (df['intron start']-df['intron end']).apply(abs)
    
    for bam, data in data_dict.iteritems():
        df[data[bam]+': reads in transcript'] = data['reads in transcript']
        df[data[bam]+': reads in intron'] = data['reads in intron']
        df[data[bam]+': intron retention'] = ((df[data[bam]+': reads in intron']/float(sum(df[data[bam]+': reads in intron'])))/
                                              (df[data[bam]+': reads in transcript']/float(sum(df[data[bam]+': reads in transcript']))))
        
    return df

def PE_fragment_size(bam_file):
    '''Calculates average and standard deviation of insert fragment sizes from paired end data. Necessary for GEO deposition
    
    Parameters
    ----------
    bam_file : str
            bam file from Bowtie or STAR from paired end data
    
    Output
    ------
    Prints the average and standard deviation of the library fragment size'''
    
    bam = pysam.Samfile(bam_file)
    sizes = []
    try:
        reads = bam.fetch('chr1',1000,20000)
    except ValueError:
        reads = bam.fetch('I',1000,20000)
    for read in reads:
        if read.is_paired:
            try:
                mate = bam.mate(read)
                if read.reference_name == mate.reference_name:
                    if read.is_reverse and not mate.is_reverse:
                        size = read.reference_end-mate.reference_start
                        if size > 0 and size < 5000:
                            sizes.append(size)
                    elif not read.is_reverse and mate.is_reverse:
                        size = mate.reference_end-read.reference_start
                        if size > 0 and size < 5000:
                            sizes.append(size)
            except ValueError:
                pass
    print "Average fragment size: "+str(np.mean(sizes))
    print "Standard deviation: "+str(np.std(sizes))

def convert_bed_to_gff3(bed, save=True):
    names = ['chromosome','start','end','name','bitscore','strand','w_start','w_end','rgb']
    df = pd.read_csv(bed, sep='\t',names=names)
    
    gff3_df = df[['chromosome','start','end','strand','name']]
    
    if str(gff3_df.loc[0,'start']) == 'nan':
        gff3_df = gff3_df.loc[1:,:]
    
    gff3_df.loc[:,'source'] = bed
    gff3_df.loc[:,'type'] = 'gene'
    gff3_df.loc[:,'x'] = '.'
    gff3_df.loc[:,'y'] = '.'
    gff3_df.loc[:,'name'] = gff3_df['name'].apply(str)
    gff3_df.loc[:,'ID'] = ['ID='+x+';' for x in gff3_df['name']]
    gff3_df.loc[:,'start'] = gff3_df['start'].apply(int)
    gff3_df.loc[:,'end'] = gff3_df['end'].apply(int)
    
    gff3_df = gff3_df[['chromosome','source','type','start','end','x','strand','y','ID']]
    
    if save:
        gff3_name = bed.split('.bed')[0]+'.gff3'
        gff3_df.to_csv(bed.split('.bed')[0]+'.gff3', sep='\t', header=False, index=False)
        return gff3_name
    else:
        return gff3_df
    
    
def count_reads_from_gff3(bam_list, gff3, stranded='no', feature_type='gene', rpkm=True, csv=None, bed=False):
    total_reads = {}
    
    if bed:
        gff3 = convert_bed_to_gff3(gff3)
    
    for bam in bam_list:
        total_reads[bam] = GT.count_aligned_reads(bam)
        print '\n'+bam
        print total_reads[bam]
        args = 'htseq-count -q -f bam -s {0} -t {1} -i ID {2} {3}'.format(stranded, feature_type, bam, gff3)
        print args
        
        out = check_output(args.split(' '))
        with open(bam.split('_sorted')[0]+'.htseq','w') as fout:
            fout.write(out)
    
    htseq_df = None
    if rpkm:
        transcripts = GT.populate_transcripts(gff3, gff3_class=feature_type)
    for bam in bam_list:
        htseq = bam.split('_sorted')[0]+'.htseq'
        name = htseq.split('.htseq')[0]
        
        if htseq_df is None:
            htseq_df = pd.read_csv(htseq, sep='\t', index_col=0, names=[name])
            if rpkm:
                tx_lengths = []
                for ix, r in htseq_df.iterrows():
                    tx_lengths.append(transcripts[ix].length/1000.)
                htseq_df.loc[:,'Transcript length (kb)'] = tx_lengths
        else:
            new_htseq_df = pd.read_csv(htseq, sep='\t', index_col=0, names=[name])
            htseq_df = htseq_df.merge(new_htseq_df, how='outer', right_index=True, left_index=True)
        
        htseq_df.loc[:,name] = htseq_df[name].divide(total_reads[bam])
        if rpkm:
            htseq_df.loc[:,name] = htseq_df[name]/htseq_df['Transcript length (kb)']
    
    if csv is not None:
        htseq_df.to_csv(csv)
        return None
    else:
        return htseq_df
            
            
            
        