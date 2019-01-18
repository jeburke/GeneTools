import os
import sys
import pysam
import pandas as pd
import numpy as np
from scipy import stats
import json
import pickle

script_path = os.path.dirname(os.path.realpath(__file__)).split('GeneTools')[0]
sys.path.append(script_path)
import GeneTools as GT
    
class Intron:
    '''Intron object containing information about the intron location and size
    
    Parameters
    ----------
    chromosome_or_str : str
        If using the chromosome, must also provide start, end and strand
        If using a string, it must be a known string format containing the same information (see below)
    start : int, default `None`
        Beginning of intron on the chromosome (5 prime splice site on + strand, 3 prime splice site on - strand)
    end : int, default `None`
        End of intron on the chromosome
    strand : str, default `None`
        Strand of transcript containing intron
    str_format : str, default "JUM"
        Format of the string that contains the intron information
        "JUM" format is chrom_strand_start_end
        "Jordan" format is chrom:start-end,strand
        No others are currently supported
    organism : str or GT.Organism object, default crypto
        If using str - options are crypto, pombe, candida and cerevisiae
        - will create object from default files every time an Intron instance is initialized
        For better perforance - first initialize GT.organism and pass in the object
        
    Examples
    -------
    Constructing an intron from separate items
    >>> intron = Intron("chr1",611500,611563,"+")
    >>> intron.length
        63
    
    Constructing an intron from a string
    >>> intron = Intron("chr1_+_611500_611563", str_format="JUM")
    >>> intron.length
        63
        
    Providing an organism object
    >>> crypto = GT.Organism('Cryptococcus neoformans H99', 
                           '/home/jordan/GENOMES/H99_fa.json', 
                           '/home/jordan/GENOMES/CNA3_FINAL_CALLGENES_1_gobs.gff3')
    >>> intron = Intron("chr1",611500,611563,"+", organism=crypto)
    >>> intron.organism.name
        Cryptococcus neoformans H99
    
    '''
    def __init__(intron, chromosome_or_str, start=None, end=None, strand=None, str_format=None, organism='crypto'):
        if str_format is None:
            chromosome = chromosome_or_str
            if start is None or end is None or strand is None:
                raise ValueError('Must provide start, end and strand or indicate str_format')
        else:
            chromosome, start, end, strand = intron_from_string(chromosome_or_str, str_format=str_format)
            
        intron.chromosome = chromosome
        intron.start = int(start)
        intron.end = int(end)
        intron.strand = strand
        intron.length = abs(end-start)
        
        if type(organism) == str:
            organism = GT.autoload_organism(organism)
        intron.organism = organism
            
    def transcripts(intron, as_string=False):
        ''' Find transcripts that may contain this intron (same strand and location)
        
        Parameters
        ----------
        as_string : bool, default `False`
            If True, returns lists of transcripts as comma separated strings
        
        Returns
        -------
        tx_list : list or str
            list of transcripts if as_string is False
            comma separated string of transcripts if as_string is True
        '''
        ann_df = intron.organism.annotation
        
        if intron.organism.name == 'Schizosaccharomyces pombe':
            transcripts = ann_df[ann_df['type'] == 'transcript']
        else:
            transcripts = ann_df[ann_df['type'] == 'mRNA']
        matches = transcripts[(transcripts['chromosome'] == intron.chromosome) & (transcripts['strand'] == intron.strand)]
        matches = matches[(matches['start'] < intron.start) & (matches['end'] > intron.end)]
        if len(matches) > 0:
            if intron.organism.name == 'Schizosaccharomyces pombe':
                tx_list = [x.lstrip('ID=transcript:').split(';')[0] for x in matches['ID']] 
            else:
                tx_list = [x.lstrip('ID=').split(';')[0] for x in matches['ID']]
                if 'mRNA' in tx_list[0]:
                    tx_list = [x.rstrip('_mRNA') for x in tx_list]
        else:
            tx_list = None
        if as_string:
            if tx_list is not None:
                return ','.join(tx_list)
            else:
                return ''
        else:
            return tx_list

    def annotated(intron):
        ''' Determines whether the intron is in the current annotation
        
        Returns
        -------
        ann : bool
            If True, intron is annotated
            If False, intron is not annotated
        '''
        intron_df = build_intron_df(intron.organism.gff3, transcript_list=intron.transcripts())
        if intron.start in list(intron_df['start']) and intron.end in list(intron_df['end']):
            ann = True
        else:
            ann = False
        return ann
    
    def position(intron):
        ''' Finds the position of the intron in the transcript

        Returns
        -------
        position : int
            Zero indexed position of the intron within the transcript
            I.e. the first intron will be 0 etc.
        '''
        
        intron_df = build_intron_df(intron.organism.gff3, transcript_list=intron.transcripts())
        
        starts = list(intron_df['start']).append(intron.start)
        starts = sorted(list(set(starts)))
        
        position = starts.index(intron.start)
        return position
    
    def is_first(intron):
        ''' Determines whether the intron is the first intron in the transcript
        
        Returns
        -------
        ann : bool
            If True, intron is the first intron
        '''
        
        intron_df = build_intron_df(intron.organism.gff3, transcript_list=intron.transcripts())
        
        starts = list(intron_df['start'])
        starts.append(intron.start)
        starts = sorted(list(set(starts)))
        
        first = False
        if starts.index(intron.start) == 0:
            first = True
        return first
    
    def is_last(intron):
        ''' Determines whether the intron is the last intron in the transcript
        
        Returns
        -------
        ann : bool
            If True, intron is the last intron
        '''
        intron_df = build_intron_df(intron.organism.gff3, transcript_list=intron.transcripts())
        
        starts = list(intron_df['start'])
        starts.append(intron.start)
        starts = sorted(list(set(starts)))
        
        last = False
        if starts.index(intron.start) == len(starts)-1:
            last = True
        return last
    
    def is_only(intron):
        ''' Determines whether the intron is the only intron in the transcript
        
        Returns
        -------
        ann : bool
            If True, intron is the only intron
        '''
        intron_df = build_intron_df(intron.organism.gff3, transcript_list=intron.transcripts())
        
        starts = list(intron_df['start'])
        starts.append(intron.start)
        only = False
        if len(set(starts)) == 1:
            only = True
        return only
    
    def UTR_5prime(intron):
        ''' Determines whether the intron is in the 5 prime UTR
        
        Returns
        -------
        ann : bool
            If True, intron is in the 5 prime UTR
        '''
        transcript_dict = GT.populate_transcripts(intron.organism.gff3, organism=intron.organism)
        UTR = False
        for iso in intron.transcripts():
            tx = transcript_dict[iso]
            UTR_start, UTR_end = tx.UTR_5prime()
            if intron.start in range(UTR_start, UTR_end):
                UTR = True
        return UTR
    
    def UTR_3prime(intron):
        ''' Determines whether the intron is in the 3 prime UTR
        
        Returns
        -------
        ann : bool
            If True, intron is in the 3 prime UTR
        '''
        transcript_dict = GT.populate_transcripts(intron.organism.gff3, organism=intron.organism)
        UTR = False
        for iso in intron.transcripts():
            tx = transcript_dict[iso]
            UTR_start, UTR_end = tx.UTR_3prime()
            if intron.start in range(UTR_start, UTR_end):
                UTR = True
        return UTR
        
        
    ################################################
    ## Methods for obtaining sequence and scoring ##
    ################################################
    
    def sequence(intron, seq_range=(-2,2)):
        ''' Retrieve the intron sequence
        
        Parameters
        ----------
        seq_range : tuple, default (-2,2)
            Range of sequence to retrieve relative to start and end of intron. 
            Default retrieves 2 nt upstream and 2 nt downstream.
        
        Returns
        -------
        seq : str
            The sequence of the intron in the specified range
        '''
        
        fa_dict = intron.organism.sequence
        
        seq = fa_dict[intron.chromosome][intron.start+seq_range[0]-1:intron.end+seq_range[1]]
        if intron.strand == '-':
            seq = reverse_complement(seq)
        return seq
    
    def score5p(intron, PSSM_txt_file, position=(-2,6), quiet=True):
        ''' Score the 5 prime end of the intron
        
        Parameters
        ----------
        PSSM_txt_file : str
            Text file created from a numpy array containing the PSSM scores from annotated introns.
            This file can be created using the build_consensus_matrix function in this module.
        position : tuple, default (-2,6)
            Range of sequence to retrieve relative to start of intron
            Default retrieves 2 nt upstream and 6 nt downstream (e.g. TGGTAAGT)
        quiet : bool, default `True`
            Do not print the 5 prime splice site sequence if True 
            (use quiet when iterating through large numbers of introns)
        
        Returns
        -------
        score_5prime : float
            The score of the 5 prime splice site
            Higher scores are closer to consensus. A score of 0 is random.
        '''
        
        PSSM = np.loadtxt(PSSM_txt_file)
        base_dict = {"A":0, "C":1, "T":2, "G":3}
        seq = intron.sequence(seq_range=(position[0], position[0]*-1))
        seq5 = seq[:position[1]-position[0]]
        if not quiet:
            print 'Sequence: '+seq5
        score_5prime = 0
        for a, base in enumerate(seq5):
            score_5prime += PSSM[base_dict[base],a]
        return score_5prime
    
    def score3p(intron, PSSM_txt_file, position=(-6,2), quiet=True):
        ''' Score the 3 prime end of the intron
        
        Parameters
        ----------
        PSSM_txt_file : str
            Text file created from a numpy array containing the PSSM scores from annotated introns.
            This file can be created using the build_consensus_matrix function in this module.
        position : tuple, default (-6,2)
            Range of sequence to retrieve relative to start of intron
            Default retrieves 2 nt upstream and 6 nt downstream (e.g. TATTAGAC)
        quiet : bool, default `True`
            Do not print the 3 prime splice site sequence if True 
            (use quiet when iterating through large numbers of introns)
        
        Returns
        -------
        score_3prime : float
            The score of the 3 prime splice site
            Higher scores are closer to consensus. A score of 0 is random.
        '''
        PSSM = np.loadtxt(PSSM_txt_file)
        base_dict = {"A":0, "C":1, "T":2, "G":3}
        
        seq = intron.sequence(seq_range=(position[1]*-1, position[1]))
        seq3 = seq[position[0]-position[1]:]
        if not quiet:
            print 'Sequence: '+seq3
        score_3prime = 0
        for b, base in enumerate(seq3):
            score_3prime += PSSM[base_dict[base],b]
        return score_3prime
            
    def percent_py(intron, fraction=0.3):
        ''' Calculate the percent pyrimidine in the last 30% of the intron
        
        Parameters
        ----------
        fraction : float, default 0.3
            Fraction of the intron over which to calculate percent Py. Uses the last 30% by default.
        
        Returns
        -------
        perc_py : float
            Percent pyrimidine content in the specified fraction of the intron 
        '''
        seq = intron.sequence()
        start_ix = len(seq)-int(len(seq)*fraction)
        py_count = 0
        for n, base in enumerate(seq[start_ix:]):
            if base == 'T' or base == 'C':
                py_count += 1
        perc_py = py_count/float(n+1)*100
        return perc_py
        
        
    ############################################
    ## Methods for finding branch information ##
    ############################################
    
    def branch(intron, branch_db, guess=True):
        ''' Retrieve branch information if available
        
        Parameters
        ----------
        branch_db : str
            Pickle or csv file containing available branch information for introns in the organism.
            Must contain the following columns: [chromosome, start, end, strand, branch, sequence].
            Chromosome, start, end and strand refer to the intron.
            Branch is the position of the branch in the chromosome.
            Sequence is the 5 nt surrounding the branch where the branch A is position 4 (e.g. CTAAC)
        guess : bool, default `True`
            If no branch was experimentally detected for the intron, guess based on the most 
            common branch sequences in the database.
            
        Returns
        -------
        branch : int
            position of the branch in the chromosome
        br_seq : str
            sequence of the branch (e.g. CTAAC)
        br_dist : int
            distance from the branch to the 3 prime splice site
        py_content : float
            fraction pyrimidine between the branch and the 3 prime splice site
        '''
        if branch_db.endswith('pickle') or branch_db.endswith('.db'):
            branch_db = pd.read_pickle(branch_db)
        elif branch_db.endswith('csv'):
            branch_db = pd.read_csv(branch_db, index_col=0)
        else:
            print "Unrecognized branch database file type"
            return None

        seq_freq = []
        for seq in set(branch_db['sequence']):
            seq_freq.append((seq, len(branch_db[branch_db['sequence'] == seq])))

        seq_freq = sorted(seq_freq, key=lambda x: x[1])
        seq_freq.reverse()
        branches = [x[0] for x in seq_freq if x[1] > 5]
        
        intron_seq = intron.sequence()
        intron_br = branch_db[(branch_db['start'] == intron.start) & 
                              (branch_db['chromosome'] == intron.chromosome) &
                              (branch_db['end'] == intron.end) &
                              (branch_db['strand'] == intron.strand)]

        if len(intron_br) == 1:
            branch = intron_br['branch']
            br_seq = intron_br['sequence']
            if intron.strand == '+':
                br_dist = intron_br['end']-branch
            elif intron.strand == '-':
                br_dist = branch-intron_br['start']
        else:
            branch = None
            br_dist = np.NaN
            br_seq = None
            for br in branches:
                ix = intron_seq.rfind(br)
                if ix > 5:
                    if intron.strand == '+':
                        branch = intron.start+ix
                    elif intron.strand == '-':
                        branch = intron.end-ix
                    br_seq = br
                    br_dist = len(intron_seq)-ix
                    break
        if branch is not None:
            py_count = 0
            br_3p_seq = intron_seq[len(intron_seq)-br_dist:]
            for n, base in enumerate(br_3p_seq):
                if base == 'C' or base == 'T':
                    py_count += 1
            py_content = py_count/float(n+1)

        else:
            py_content = np.NaN
        return branch, br_seq, br_dist, py_content
    
    
    ################################
    ## Methods for counting reads ##
    ################################
    
    def count_5pss_reads(intron, bams, rpm=False, library_direction='reverse'):
        ''' Count the reads (or calculate RPM) spanning the 5 prime splice site
        
        Parameters
        ----------
        bams : str or list
            Single or list of sorted, indexed bam files for analysis
        rpm : bool, default `False`
            Normalize to total aligned reads (for comparing between samples)
        library_direction : str, default reverse
            Direction of first read relative to genome (usually reverse for RNA-seq)
        
        Returns
        -------
        counts : pandas.core.series.Series
            Pandas series where the index is the name of the bam files and the value is the read count or RPM
            Can be treated like a dictionary, but easier to perform math or populate a DataFrame
        '''
        counts = span_read_counter(bams, intron.chromosome, intron.start, intron.strand, rpm, library_direction)
        return counts

    def count_3pss_reads(intron, bams, rpm=False, library_direction='reverse'):
        ''' Count the reads (or calculate RPM) spanning the 3 prime splice site
        
        Parameters
        ----------
        bams : str or list
            Single or list of sorted, indexed bam files for analysis
        rpm : bool, default `False`
            Normalize to total aligned reads (for comparing between samples)
        library_direction : str, default reverse
            Direction of first read relative to genome (usually reverse for RNA-seq)
        
        Returns
        -------
        counts : pandas.core.series.Series
            Pandas series where the index is the name of the bam files and the value is the read count or RPM
            Can be treated like a dictionary, but easier to perform math or populate a DataFrame
        '''
        counts = span_read_counter(bams, intron.chromosome, intron.end, intron.strand, rpm, library_direction)
        return counts
    
    def count_intronic_reads(intron, bams, rpkm=False, library_direction='reverse'):
        ''' Count the reads (or calculate RPKM) within the intron.
            Will exclude reads with junctions that match the intron, but will 
            count reads that have other junctions (i.e. nested splicing events).
        
        Parameters
        ----------
        bams : str or list
            Single or list of sorted, indexed bam files for analysis
        rpkm : bool, default `False`
            Normalize to total aligned reads (for comparing between samples) and intron size
        library_direction : str, default reverse
            Direction of first read relative to genome (usually reverse for RNA-seq)
        
        Returns
        -------
        counts : pandas.core.series.Series
            Pandas series where the index is the name of the bam files and the value is the read count or RPM
            Can be treated like a dictionary, but easier to perform math or populate a DataFrame
        '''
        counts = simple_read_counter(bams, intron.chromosome, intron.start, intron.end, intron.strand, rpkm, library_direction)
        return counts
    
    def count_transcript_reads(intron, bams, rpkm=False, library_direction='reverse'):
        ''' Count the reads (or calculate RPKM) within the transcript that contains the intron.
        
        Parameters
        ----------
        bams : str or list
            Single or list of sorted, indexed bam files for analysis
        rpkm : bool, default `False`
            Normalize to total aligned reads (for comparing between samples) and transcript size
        library_direction : str, default reverse
            Direction of first read relative to genome (usually reverse for RNA-seq)
        
        Returns
        -------
        counts : dict of pandas.core.series.Series
            Dictionary of pandas series where the key is the transcript name.
            Values are andas series where the index is the name of the bam files and the value is the read count or RPM.
            Can be treated like a dictionary, but easier to perform math or populate a DataFrame.
        '''
        transcript_dict = GT.populate_transcripts(intron.organism.gff3, organism=intron.organism)
        tx_counts = {}
        for iso in intron.transcripts():
            tx = transcript_dict[iso]
            tx_counts[iso] = simple_read_counter(bams, tx.chromosome, tx.start, tx.end, tx.strand, rpkm, library_direction)
        return tx_counts
    
    def count_junction_reads(intron, bams, rpm=False, library_direction='reverse'):
        ''' Count the reads with junctions that match the intron.
        
        Parameters
        ----------
        bams : str or list
            Single or list of sorted, indexed bam files for analysis
        rpm : bool, default `False`
            Normalize to total aligned reads (for comparing between samples)
        library_direction : str, default reverse
            Direction of first read relative to genome (usually reverse for RNA-seq)
        
        Returns
        -------
        counts : pandas.core.series.Series
            Pandas series where the index is the name of the bam files and the value is the read count or RPM
            Can be treated like a dictionary, but easier to perform math or populate a DataFrame
        '''
        counts = pd.Series(index=bams)
        for bam in bams:
            read_count = 0
            reads = pysam.Samfile(bam).fetch(intron.chromosome, intron.start-200, intron.end+200)
            for read in reads:
                # Filter for reads with junctions (can be 3 or 5 entries in cigartuples)
                if len(read.cigartuples) >= 3 and read.cigartuples[1][0] == 3:
                    junction_start = read.reference_start + read.cigartuples[0][1]+1
                    junction_end = junction_start + read.cigartuples[1][1]-1
                    
                    # Make sure the junctions match the intron
                    if intron.start == junction_start and intron.end == junction_end:
                        if library_direction == 'reverse':
                            if read.is_reverse and intron.strand == '+':
                                read_count += 1
                            elif not read.is_reverse and intron.strand == '-':
                                read_count += 1
                        else:
                            if not read.is_reverse and intron.strand == '+':
                                read_count += 1
                            elif read.is_reverse and intron.strand == '-':
                                read_count += 1
                    
                    # Check if there is a second junction in the read that might match the intron
                    elif len(read.cigartuples) == 5:
                        junction_start2 = junction_end + read.cigartuples[2][1]
                        junction_end2 = junction_start2 + read.cigartuples[3][1]
                        if intron.start == junction_start2 and intron.end == junction_end2:
                            if library_direction == 'reverse':
                                if read.is_reverse and intron.strand == '+':
                                    read_count += 1
                                elif not read.is_reverse and intron.strand == '-':
                                    read_count += 1
                            else:
                                if not read.is_reverse and intron.strand == '+':
                                    read_count += 1
                                elif read.is_reverse and intron.strand == '-':
                                    read_count += 1
            counts[bam] = read_count
        return counts              
            
###############################################
## Functions for populating intron instances ##
###############################################

def build_intron_df(gff3, transcript_list=None):
    if 'pombe' in gff3: pombe = True
    else: pombe = False   
    
    ann_df = GT.read_gff3(gff3)
    if 'intron' in list(ann_df['type']):
        print "Introns in gff3"
        intron_df = ann_df[ann_df['type'] == 'intron']
        intron_df = intron_df[['chromosome','start','end','strand','ID']]
        intron_df.loc[:,'transcript'] = intron_df['ID'].str.split('Parent=').str[1].str.split(';').str[0]
        if 'cerev' in gff3:
            intron_df.loc[:,'transcript'] = intron_df.loc['transcript'].str.split('_mRNA').str[0]
        if transcript_list is not None:
            intron_df = intron_df['transcript'].isin(transcript_list)
    
    elif 'exon' in list(ann_df['type']):
        exon_df = ann_df[['chromosome','start','type','end','strand','ID']]
        exon_df = exon_df[exon_df['type'] == 'exon']
        exon_df.loc[:,'transcript'] = exon_df['ID'].str.split('Parent=').str[1].str.split(';').str[0]
        if pombe:
            exon_df.loc[:,'transcript'] = exon_df['transcript'].str.split('transcript:').str[1]
        if transcript_list is not None:
            exon_df = exon_df[exon_df['transcript'].isin(transcript_list)]
        intron_df = pd.DataFrame(columns = exon_df.columns)
        df_dict = {'chromosome':[],'transcript':[],'start':[],'end':[],'strand':[]}
        for tx in set(exon_df['transcript']):
            tx_df = exon_df[exon_df['transcript'] == tx]
            tx_df = tx_df.sort_values('start').reset_index(drop=True)
            for n in tx_df.index[:-1]:
                try:
                    df_dict['transcript'].append(tx)
                    df_dict['chromosome'].append(tx_df.loc[n,'chromosome'])
                    df_dict['strand'].append(tx_df.loc[n,'strand'])
                    df_dict['start'].append(tx_df.loc[n,'end']+1)
                    df_dict['end'].append(tx_df.loc[n+1,'start']-1)
                except KeyError:
                    pass
        for col, data in df_dict.iteritems():
            intron_df.loc[:,col] = data
    
    intron_df['start'] = intron_df['start'].apply(int)
    intron_df['end'] = intron_df['end'].apply(int)
    return intron_df
    
def intron_from_string(intron_str, str_format="JUM"):
    if str_format == "JUM":
        chromosome = intron_str.split('_')[0]
        strand = intron_str.split('_')[1]
        start = int(intron_str.split('_')[2])
        end = int(intron_str.split('_')[3].strip())
    elif str_format == "Jordan":
        chromosome = intron_str.split(':')[0]
        strand = intron_str.split(',')[-1].strip()
        start = int(intron_str.split(':')[-1].split(',')[0].split('-')[0])
        end = int(intron_str.split(':')[-1].split(',')[0].split('-')[1])
    return chromosome, start, end, strand

def introns_from_gff3(organism, transcript_list=None, by_transcript=False, sort=True, pickle_name=None):
    '''Function to build a list of intron objects in the genome
    
    Parameters
    ----------
    gff3_file : str
        location and name of gff3 file
    transcript_list : list, default `None`
        List of genes to limit the dictionary. Useful for optimizing performance of downstream functions.
    by_transcript : bool, default `False`
        Return a dictionary of introns in each transcript, rather than a simple list
    sort : bool, default `True`
        Sort introns first by chromosome then by start position. 
        If False will be in order of occurence in gff3
                
    Returns
    --------
    introns : list of intron objects if by_transcript is False
              dict of intron objects with transcripts as keys if by_transcript is True
                    '''
    if type(organism) == str:
        organism = GT.autoload_organism(organism)
 
    intron_df = build_intron_df(organism.gff3, transcript_list=transcript_list)
    
    #first group introns to remove redundancy
    intron_df = intron_df.sort_values('transcript')
    intron_df = intron_df.drop_duplicates(subset=['chromosome','start','end'])
    
    #then go through all introns and create instances
    if not by_transcript:
        introns = []
        for ix, r in intron_df.iterrows():
            introns.append(Intron(r['chromosome'],r['start'],r['end'],r['strand'], organism=organism))
        if sort:
            introns.sort(key=lambda intron: (intron.chromosome, intron.start))
    else:
        introns = {}
        for tx in set(intron_df['transcript']):
            tx_df = intron_df[intron_df['transcript'] == tx]
            introns[tx] = []
            for ix, r in tx_df.iterrows():
                introns[tx].append(Intron(r['chromosome'],r['start'],r['end'],r['strand'],organism=organism))
            if sort:
                introns[tx].sort(key=lambda intron: intron.start)
    if pickle_name is not None:
        with open(pickle_name, 'w') as fout:
            pickle.dump(introns, fout)
    return introns

##########################################
## Functions for working with sequences ##
##########################################
    
def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N', 'Y':'R', 'R':'Y'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
        return complement(s[::-1])

def gc_content(fa_dict):
    base_dict = {0:"A", 1:"C", 2:"T", 3:"G"}
    tally = np.zeros([1,4])[0]
    for chrom, seq in fa_dict.iteritems():
        for n in range(len(tally)):
            tally[n] += seq.count(base_dict[n])
    total = float(sum(tally))   
    nucleotide_prob = tally/total
    return nucleotide_prob

def build_consensus_matrix(gff3, fa, out_name, seq_type='intron', position=('5prime',-2,6)):
    if fa.endswith('.json'):
        with open(fa) as f: fa_dict = json.load(f)
    else:
        fa_dict = make_fasta_json(fasta_file)
    
    all_transcripts = GT.populate_transcripts(gff3)
    
    nuc_prob = gc_content(fa_dict)
    base_dict = {"A":0, "C":1, "T":2, "G":3}
    
    # 1st row is A counts, second row is C, third row is T, fourth row is G.
    array_len = position[2]-position[1]
    PSSM = np.zeros([4,array_len])

    sequences = []
    print "Loading sequences..."
    if seq_type == 'intron':
        all_introns = introns_from_gff3(gff3)
        for intron in all_introns:
            if position[0] == '5prime':
                seq = intron.sequence(fa_dict, seq_range=(position[1]-1, position[1]*-1))
            elif position[0] == '3prime':
                seq = intron.sequence(fa_dict, seq_range=(position[2]*-1-1, position[2]))
            else:
                print "unrecognized sequence position type"
                return None
            sequences.append(seq)
    
    elif seq_type == 'transcript':
        for tx in all_transcripts:
            seq = tx.sequence(fa_dict, seq_range=seq_range)
            sequences.append(seq)
            
    print "Calculating base composition..."
    for seq in sequences:
        if position[0] == '5prime':
            seq = seq[:position[2]-position[1]]
        elif position[0] == '3prime':
            seq = seq[position[1]-position[2]:]

        for a, base in enumerate(seq):
            PSSM[base_dict[base],a] += 1

    for a in range(PSSM.shape[0]):
        for b in range(PSSM.shape[1]):
            if PSSM[a,b] == 0: PSSM[a,b] += 1
            PSSM[a,b] = np.log2((PSSM[a,b]/float(len(sequences)))/nuc_prob[a])
    
    float_formatter = lambda x: "%.1f" % x
    np.set_printoptions(formatter={'float_kind':float_formatter})
    print PSSM
    
    np.savetxt(out_name, PSSM)
    
    
##################################
## Functions for counting reads ##
##################################
            
def span_read_counter(bams, chromosome, coordinate, strand, rpm, library_direction):
    if type(bams) == str:
        bams = [bams]
    counts = pd.Series(index=bams)
    for bam in bams:
        ob = pysam.Samfile(bam)
        read_iter = ob.fetch(chromosome, coordinate-200, coordinate+200)
        read_count = 0
        for read in read_iter:
            if len(read.cigartuples) < 3:
                if read.reference_start < coordinate and read.reference_end > coordinate:
                    if library_direction == 'reverse':
                        if strand == '+' and read.is_reverse:
                            read_count += 1
                        elif strand == '-' and not read.is_reverse:
                            read_count += 1
                    else:
                        if strand == '+' and not read.is_reverse:
                            read_count += 1
                        elif strand == '-' and read.is_reverse:
                            read_count += 1
        counts[bam] = read_count
        
    if rpm:
        for bam in counts.keys():
            scale = GT.count_aligned_reads(bam)
            counts[bam] = counts[bam]*scale
    return counts

def simple_read_counter(bams, chromosome, start, end, strand, rpkm, library_direction):
    if type(bams) == str:
        bams = [bams]
    counts = pd.Series(index=bams)
    for bam in bams:
        ob = pysam.Samfile(bam)
        read_iter = ob.fetch(chromosome, start-200, end+200)
        read_count = 0
        for read in read_iter:
            if read.reference_end >= start and read.reference_start <= end:
                if library_direction == 'reverse':
                    if strand == '+' and read.is_reverse:
                        read_count += 1
                    elif strand == '-' and not read.is_reverse:
                        read_count += 1
                else:
                    if strand == '+' and not read.is_reverse:
                        read_count += 1
                    elif strand == '-' and read.is_reverse:
                        read_count += 1
        counts[bam] = read_count
    if rpkm:
        for bam in counts.keys():
            total = GT.count_aligned_reads(bam)
            counts[bam] = (counts[bam]/total)/((end-start)/1000)
    return counts        
        



    