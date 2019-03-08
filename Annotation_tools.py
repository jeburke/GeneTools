import pandas as pd
import json
import os
import sys
script_path = os.path.dirname(os.path.realpath(__file__)).split('GeneTools')[0]
sys.path.append(script_path)
import GeneTools as GT

class Organism:
    def __init__(self, name, fasta_file, gff3):
        self.name = name
        self.fasta_file = fasta_file
        self.gff3 = gff3
        if fasta_file.endswith('.json'):
            with open(fasta_file) as f: fa = json.load(f)
        else:
            fa = make_fasta_json(f)    
        self.sequence = fa
        self.chromosome_lengths = {k:len(v) for k,v in fa.items()}
        self.annotation = read_gff3(gff3)

class Transcript:
    def __init__(transcript, name, chromosome, start, end, strand, organism='crypto'):
        transcript.name = name
        transcript.chromosome = chromosome
        transcript.start = start
        transcript.end = end
        transcript.strand = strand
        transcript.length = end-start
        if type(organism) == str:
            organism = autoload_organism(organism)
        transcript.organism = organism
        
    def CDS(transcript):
        df = transcript.organism.annotation
        CDS_df = df[df['type'] == 'CDS']
        CDS_df.loc[:,'transcript'] = CDS_df['ID'].str.split('Parent=').str[1].str.split(';').str[0]
        if transcript.organism.name == 'Schizosaccharomyces pombe':
            CDS_df.loc[:,'transcript'] = CDS_df['transcript'].str.split('transcript:').str[-1]
        elif transcript.organism.name == 'Saccharomyces cerevisiae S288C':
            CDS_df.loc[:,'transcript'] = CDS_df['transcript'].str.split('_mRNA').str[0]
       
        tx_CDS = CDS_df[CDS_df['transcript'] == transcript.name]
        CDS_start = min(tx_CDS['start'])
        CDS_end = max(tx_CDS['end'])
        
        return CDS_start, CDS_end
    
    def UTR_5prime(transcript):
        CDS_start, CDS_end = transcript.CDS()
        if transcript.strand == '+':
            return transcript.start, CDS_start
        elif transcript.strand == '-':
            return CDS_end, transcript.end
    
    def UTR_3prime(transcript):
        CDS_start, CDS_end = transcript.CDS()
        if transcript.strand == '+':
            return CDS_end, transcript.end
        elif transcript.strand == '-':
            return transcript.start, CDS_start
        
    def sequence(transcript, seq_range=(0,0)):
        fa_dict = transcript.organism.sequence
        seq = fa_dict[transcript.chromosome][transcript.start+seq_range[0]:transcript.end+seq_range[1]]
        if transcript.strand == '-':
            seq = GT.reverse_complement(seq)
        return seq
        
    def transcript_count_reads(transcript, bams, library_direction='reverse'):
        counts = GT.simple_read_counter(bams, transcript.chromosome, transcript.start, transcript.end, transcript.strand, library_direction=library_direction)
        return counts
    
    def UTR_5prime_count_reads(transcript, bams, rpkm=False, library_direction='reverse'):
        UTR = transcript.UTR_5prime()
        counts = GT.simple_read_counter(bams, transcript.chromosome, UTR[0], UTR[1], transcript.strand, rpkm=rpkm, library_direction=library_direction)
        return counts
    

#####################################################
## Tools for reading annotation files and creating ##
## instances of Transcript and Organism            ##
#####################################################

def read_gff3(gff3):
    names = ['chromosome','source','type','start','end','x','strand','y','ID']
    ann_df = pd.read_csv(gff3, sep='\t', names=names)
    return ann_df

def make_fasta_json(fa, write=True):
    '''Convert fasta file to a dictionary stored in json format
    
    Parameters
    ----------
    fa : fasta file (e.g. genome fasta downloaded from fungidb)
    write : bool, default `True`
        write the json file - if False will just return the dictionary
    '''
    
    fa_dict = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                contig = line.lstrip('>').strip()
                fa_dict[contig] = ''
            else:
                fa_dict[contig] = fa_dict[contig] + line.strip()
    if write:
        fa_name = fa.split('.fa')[0]+'_fa.json'
        with open(fa_name, 'w') as fout: json.dump(fa_dict, fout)
        return fa_name
    else:
        return fa_dict
    
def populate_transcript_df(gff3, gff3_class='mRNA', organism='crypto'):
    ann_df = GT.read_gff3(gff3)
    transcripts = ann_df[ann_df['type'] == gff3_class]
    if len(transcripts) == 0:
        print "gff3_class not found in gff3 file!"
        return None
            
    transcripts.loc[:,'transcript'] =transcripts['ID'].str.split('ID=').str[1].str.split(';').str[0]
    transcripts = transcripts.reset_index(drop=True)
    if '_mRNA' in transcripts.loc[0,'transcript']:
        transcripts.loc[:,'transcript'] =transcripts.loc['transcript'].str.split('_mRNA').str[0]
    if 'transcript:' in transcripts.loc[0,'transcript']:
        transcripts.loc[:,'transcript'] =transcripts['transcript'].str.split('transcript:').str[1]
    return transcripts

def populate_transcripts(gff3, gff3_class='mRNA', organism='crypto'):
    ann_df = GT.read_gff3(gff3)
    transcripts = ann_df[ann_df['type'] == gff3_class]
    if len(transcripts) == 0:
        print "gff3_class not found in gff3 file!"
        return None

    transcripts.loc[:,'transcript'] =transcripts['ID'].str.split('ID=').str[1].str.split(';').str[0]
    transcripts = transcripts.reset_index(drop=True)
    if '_mRNA' in transcripts.loc[0,'transcript']:
        transcripts.loc[:,'transcript'] =transcripts.loc['transcript'].str.split('_mRNA').str[0]
    if 'transcript:' in transcripts.loc[0,'transcript']:
        transcripts.loc[:,'transcript'] =transcripts['transcript'].str.split('transcript:').str[1]

    tx_dict = {}
    for ix, r in transcripts.iterrows():
        tx_dict[r['transcript']] = GT.Transcript(r['transcript'], r['chromosome'], r['start'], r['end'], r['strand'], organism=organism)

    return tx_dict

def autoload_organism(organism):
    if 'crypto' in organism.lower() or 'neoformans' in organism.lower() :
        org_obj = Organism('Cryptococcus neoformans H99', 
                           '/home/jordan/GENOMES/H99_fa.json', 
                           '/home/jordan/GENOMES/CNA3_FINAL_CALLGENES_1_gobs.gff3')
    elif 'pombe' in organism.lower():
        org_obj = Organism('Schizosaccharomyces pombe',
                           '/home/jordan/GENOMES/POMBE/Sp_fasta_dict.json',
                           '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3')
    elif 'cerev' in organism.lower() or 'S288C' in organism:
        org_obj = Organism('Saccharomyces cerevisiae S288C',
                           '/home/jordan/GENOMES/S288C/S288C_fasta_dict.json',
                           '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3')
    elif 'candida' in organism.lower() or 'albicans' in organism.lower():
        org_obj = Organism('Candida albicans',
                           '/home/jordan/GENOMES/C_albicans_fa.json',
                           '/home/jordan/GENOMES/C_albicans_SC5314_version_A21-s02-m09-r10_features.gff')
    else:
        raise ValueError('Organism currently not supported, please create instance of Organism manually')
    
    return org_obj
                           

#############################################################
## Tools for adding homologs and orthologs to spreadsheets ##
#############################################################

def crypto_annotate(text_file, sep=',', multi_index=False):
    '''Annotates a spreadsheet or text file where the gene ID (CNAG) is the first column.
    Uses Hiten's hand annotation and S. cerevisiae and S. pombe orthologs determined by Yi.
    
    Parameters
    ----------
    text_file : str
              Your spreadsheet or text file
    sep : str, default ','
              Column separator character. Change to '\t' if tab separated.
    
    Output
    -------
    csv file : File (comma separated) with annotation appended to each row.
              '''
    
    yi = pd.read_csv('/home/jordan/GENOMES/Crypto_homologs_Spombe_S288C.csv', index_col=0)
    yi.index = [x[:-2] for x in yi.index]
    yi = yi.drop_duplicates()
    hiten = pd.read_csv('/home/chomer/Code_Base/AnnotateGenes/HitensCryptoAnnot_new.txt', sep='\t', index_col=0)
    
    interpro_df = pd.read_pickle('/home/jordan/GENOMES/InterproDomains.pickle')
    
    if multi_index is True:
        yi.columns = pd.MultiIndex.from_product(['Annotation',yi.columns])
        hiten.columns = pd.MultiIndex.from_product(['Annotation',hiten.columns])
        interpro_df.columns = pd.MultiIndex.from_product(['Annotation',interpro_df.columns])
        df = pd.read_csv(text_file, sep=sep, index_col=0, header=[0,1])
    else:
        df = pd.read_csv(text_file, sep=sep, index_col=0)
        
    if len(df.index) == 0:
        print "Empty spreadsheet, cannot annotate"
        return

    if df.index[0][-2] == 'T':
        df.index = [x[:-2] for x in df.index]
        
    
    df = df.merge(hiten, right_index=True, left_index=True, how='left')
    df = df.merge(yi, right_index=True, left_index=True, how='left')
    df = df.merge(interpro_df, right_index=True, left_index=True, how='left')
    df = df.groupby(df.index).first()
    
    print text_file.split('.')[-2]+'_annotated.csv'
    df.to_csv(text_file.split('.')[-2]+'_annotated.csv')
    
def pombe_annotate(text_file, sep=',', multi_index=False):
    '''Annotates a spreadsheet or text file where the gene ID is the first column.
    Uses Pombase annotation.
    
    Parameters
    ----------
    text_file : str
              Your spreadsheet or text file
    sep : str, default ','
              Column separator character. Change to '\t' if tab separated.
    
    Output
    -------
    csv file : File (comma separated) with annotation appended to each row.
              '''
    
    pombase = pd.read_csv('/home/jordan/GENOMES/POMBE/gene_association.pombase', sep='\t', skiprows=42, index_col=0, usecols=[1,2,9,10,11])
    
    if multi_index is True:
        df = pd.read_csv(text_file, sep=sep, index_col=0, header=[0,1])
    else:
        df = pd.read_csv(text_file, sep=sep, index_col=0)
        
    if df.index[0][-2] == '.':
        df.index = [x[:-2] for x in df.index]
    
    if df.index[0].startswith('gene'):
        df.index = [x.split('gene:')[1] for x in df.index]
        
    df = df.merge(pombase, right_index=True, left_index=True, how='left')
    df = df.drop_duplicates()
    
    print text_file.split('.csv')[0]+'_annotated.csv'
    df.to_csv(text_file.split('.csv')[0]+'_annotated.csv')