import sys
import os
from subprocess import check_output
from subprocess import call
from subprocess import Popen
from multiprocessing import Pool
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

def count_aligned_reads_bedgraph(bam):
    command = 'samtools view -F 0x904 -c {0}'.format(bam)
    aligned_reads = check_output(command.split(), shell=False)
    aligned_reads = float(aligned_reads)/1000000
    print '\n'+bam
    print "{0} million reads".format("%.2f" % aligned_reads)
    return {bam:1/aligned_reads}

def run_genomeCoverageBed(command):
    with open('error.log','a') as fout:
        process = Popen(command, shell=True, stdout=fout, stderr=fout)
    ret_code = process.wait()
    return ret_code

def generate_scaled_bedgraphs2(directory, untagged, organism='crypto', start_only=False, stranded=False, threads=1, file_provided=False, expand=False):
    if 'crypto' in organism.lower():
        genome = '/home/jordan/GENOMES/crypto_for_bedgraph.genome'
    elif 'cerev' in organism.lower():
        genome = '/home/jordan/GENOMES/S288C/S288C_for_bedgraph.genome'
    elif 'pombe' in organism.lower():
        genome = '/home/jordan/GENOMES/POMBE/Sp_for_bg.genome'
    elif 'albicans' in organism.lower() or 'candida' in organism.lower():
        genome = '/home/jordan/GENOMES/C_albicans_for_bg.genome'
    else:
        genome = organism
    
    bam_list = []
    untagged_other_dir = False
    if not file_provided:
        for file in os.listdir(directory):
            if file.lower().endswith("sorted.bam") or file.endswith('sortedByCoord.out.bam'):
                print file
                bam_list.append(directory+file)
    else:
        bam_list.append(directory)
        base_dir = directory.split('/')[:-1]
        base_dir = '/'.join(base_dir)+'/'
        untagged_bams = [base_dir+x for x in os.listdir(base_dir) if untagged in x and x.endswith('.bam')]
        if len(untagged_bams) == 1:
            untagged = untagged_bams[0]
        elif len(untagged_bams) > 1:
            print "Too many matches for untagged"
            return None
        else:
            untagged = untagged
            untagged_other_dir = True
        bam_list.append(untagged)
    
    p = Pool(threads)
    totals = {}
    entries = p.map(count_aligned_reads_bedgraph, bam_list)
    for x in entries:
        totals.update(x)
    
    commands = {}
    if expand:
        flag = '-d'
    else:
        flag = '-bga'
    for n, bam in enumerate(bam_list):
        command = 'genomeCoverageBed -ibam {0} -g {1} {2} -scale {3} '.format(bam, genome, flag, "%.2f" % totals[bam])
        commands[bam] = command

    if start_only is True:
        commands = [x+'-5 ' for x in comands]
    if stranded is True:
        stranded_cmds = {}
        for bam, command in commands.iteritems():
            stranded_cmds[bam] = []
            stranded_cmds[bam].append(command+'-strand + ')
            stranded_cmds[bam].append(command+'-strand - ')
        commands = stranded_cmds
    
    final_cmds = []
    for bam, command in commands.iteritems():
        if type(command) == list:
            if untagged_other_dir and bam == untagged:
                new_name = base_dir + untagged.split('/')[-1].split('.bam')[0]
                final_cmds.append(command[0]+'> {0}_plus.bedgraph'.format(new_name))
                final_cmds.append(command[1]+'> {0}_minus.bedgraph'.format(new_name))
            else:
                final_cmds.append(command[0]+'> {0}_plus.bedgraph'.format(bam.split('.bam')[0]))
                final_cmds.append(command[1]+'> {0}_minus.bedgraph'.format(bam.split('.bam')[0]))
        else:
            if untagged_other_dir and bam == untagged:
                new_name = base_dir + untagged.split('/')[-1].split('.bam')[0]
                final_cmds.append(command+'> {0}.bedgraph'.format(new_name))
            else:
                final_cmds.append(command+'> {0}.bedgraph'.format(bam.split('.bam')[0]))
    
    p = Pool(threads)
    codes = p.map(run_genomeCoverageBed, final_cmds)
    
    return codes
        
## Deprecated
def generate_scaled_bedgraphs(directory, organism='crypto', start_only=False, stranded=False, file_provided=False):
    if 'crypto' in organism.lower():
        genome = '/home/jordan/GENOMES/crypto_for_bedgraph.genome'
    elif 'cerev' in organism.lower():
        genome = '/home/jordan/GENOMES/S288C/S288C_for_bedgraph.genome'
    elif 'pombe' in organism.lower():
        genome = '/home/jordan/GENOMES/POMBE/Sp_for_bg.genome'
    elif 'albicans' in organism.lower() or 'candida' in organism.lower():
        genome = '/home/jordan/GENOMES/C_albicans_for_bg.genome'
    
    bam_list = []
    if not file_provided:
        for file in os.listdir(directory):
            if file.lower().endswith("sorted.bam"):
                bam_list.append(directory+file)
    else:
        bam_list.append(directory)
            
    total_aligned = []
    for bam in bam_list:
        print bam
        command = 'samtools view -F 0x904 -c {0}'.format(bam)
        aligned_reads = check_output(command.split(), shell=False)
        total_aligned.append(aligned_reads)
        print "Total aligned reads in "+bam
        print aligned_reads
     
    total_aligned = [float(x)/1000000 for x in total_aligned]
    
    for n in range(len(bam_list)):
        out = bam_list[n].split('/')[-1].split('.')[0]
        if stranded is True:
            if start_only is False:
                command1 = 'genomeCoverageBed -ibam {0} -g {1} -bga -strand + -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            elif start_only is True:
                command1 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand + -5 -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            print command1
            bg1 = check_output(command1.split(), shell=False)
            with open('{0}_plus.bedgraph'.format(out),'w') as fout:
                fout.write(bg1)
            if start_only is False:
                command2 = 'genomeCoverageBed -ibam {0} -g {1} -bga -strand - -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            elif start_only is True:
                command2 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand - -5 -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            bg2 = check_output(command2.split(), shell=False)
            with open('{0}_minus.bedgraph'.format(out),'w') as fout:
                fout.write(bg2)
        else:
            if start_only is False:
                command = 'genomeCoverageBed -ibam {0} -g {1} -bga -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            else:
                command = 'genomeCoverageBed -ibam {0} -g {1} -bg -5 -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            print command
            bg = check_output(command.split(), shell=False)
            with open('{0}.bedgraph'.format(out),'w') as fout:
                fout.write(bg)            

def list_bedgraphs(directory):
    plus_list = []
    minus_list = []
    for file in os.listdir(directory):
        if file.lower().endswith("plus.bedgraph"):
            plus_list.append(directory+file)
        elif file.lower().endswith("minus.bedgraph"):
            minus_list.append(directory+file)
    plus_list.sort()
    minus_list.sort()
    bedgraphs = zip(plus_list,minus_list)
    return bedgraphs

############################################################
## Read bedgraph and sort by transcript into dictionary   ##
############################################################

def build_bedgraph_dict(transcript_dict, bedgraph_file):
    '''Function for sorting bedgraph files by gene.
    
    Parameters
    ----------
    bedgraph_file : str
                    Bedgraph file
    transcript_dict : dict
                    Transcript dict generated by build_transcript_dict (in SeqTools module)
                    
    Output
    -------
    sorted bedgraph : file
            Bedgraph file sorted by gene (ends in _by_gene.bedgraph).'''
    
    print datetime.now()
    bedgraph_dict = {}
    transcript_by_chr = {}
    for transcript, coords in transcript_dict.iteritems():
        chromosome = coords[3]
        bedgraph_dict[transcript] = [[],[]]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(transcript)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(transcript)
    
    with open(bedgraph_file, "r") as bedgraph:
        for line in bedgraph:
            columns = re.split(r'\t', line)
            bed_chr = columns[0].strip()
            rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
            if bed_chr in rom_lat:
                bed_chr = rom_lat[bed_chr]
            bed_position = int(columns[1])
            bed_peak = float(columns[3])
            
            if bed_chr in transcript_by_chr:
                transcript_list = transcript_by_chr[bed_chr]
            for transcript in transcript_list:  
                
                #Dictionary for bedgraph. Values will be [list of genomic positions][reads starting at that position]
                if bed_chr == transcript_dict[transcript][3].strip() and bed_position > transcript_dict[transcript][0] and bed_position < transcript_dict[transcript][1]:
                    bedgraph_dict[transcript][0].append(bed_position)
                    bedgraph_dict[transcript][1].append(bed_peak)
   
    with open("{0}_by_gene.bedgraph".format(bedgraph_file.split("/")[-1].split(".")[0]), "a") as fout:
        for transcript, values in bedgraph_dict.iteritems():
            fout.write(transcript+"\n")
            coord_list = map(str, bedgraph_dict[transcript][0])
            coord_line = "\t".join(coord_list)
            fout.write(coord_line+"\n")
            count_list = map(str, bedgraph_dict[transcript][1])
            count_line = "\t".join(count_list)
            fout.write(count_line+"\n")
                                
    bedgraph_dict = collections.OrderedDict(sorted(bedgraph_dict.items()))

    print datetime.now()

def read_sorted_bedgraph(bedgraph_dict_output, transcript_dict, organism=None):
    '''Function for loading sorted bedgraph files.
    
    Parameters
    ----------
    bedgraph_dict_output : str
                    File output by build_bedgraph_dict (ends in _by_gene.bedgraph)
    transcript_dict : dict
                    Transcript dict generated by build_transcript_dict (in SeqTools module)
    organism : str, default ``None``
                    change to 'pombe' if working with S. pombe
                    
    Returns
    -------
    bg_dict : dict
            Dictionary where keys are transcript name and values are pandas series of bedgraph data'''
    
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
                    if tx[-2] != 'T' and organism != 'pombe':
                        tx = tx+'T0'
                
                #Read the coordinate line
                elif n%3 == 2:
                    coords  = map(int, line.strip().split('\t'))
                
                #Read the values line
                elif n%3 == 0:
                    heights = map(float, line.strip().split('\t'))

                    if tx not in transcript_dict:
                        pass
                    else:
                        all_coords = set(range(min(coords),max(coords)))
                        missing = all_coords.difference(coords)
                        coords = coords + list(missing)

                        #Fill in missing coordinates with zeros
                        zero_fill = [0]*len(missing)
                        heights = heights + zero_fill

                        #Create a pandas series with all coordinates and sort so zeros are inserted appropriately
                        entry = pd.Series(heights, index=coords)
                        entry.sort_index(inplace=True)

                        selected_range = range(transcript_dict[tx][0],transcript_dict[tx][1])
                        entry = entry[entry.index.isin(selected_range)]
                        
                        bg_dict[tx] = entry
    return bg_dict

def decollapse_bedgraph(bedgraph):
    new_bedgraph = bedgraph.split('.bedgraph')[0]+'_full.bedgraph'
    counter1 = 0
    counter2 = 0
    with open(bedgraph) as f:
        with open(new_bedgraph,'w') as fout:
            for line in f:
                counter1 += 1
                data = line.split('\t')
                chrom = data[0]
                start = int(data[1])
                end = int(data[2])
                value = data[3].strip()

                if end-start != 1:
                    n_lines = end-start
                    for n in range(n_lines):
                        counter2 += 1
                        new_start = str(start+n)
                        new_end = str(start+n+1)
                        new_line = '\t'.join([chrom, new_start, new_end, value+'\n'])
                        fout.write(new_line)
                else:
                    counter2 += 1
                    fout.write(line)
    #print counter1
    #print counter2
    
def collapse_bedgraph(bedgraph):
    new_bedgraph = bedgraph+'.tmp'
    counter1 = 0
    counter2 = 0
    with open(bedgraph) as f:
        with open(new_bedgraph, 'w') as fout:
            for line in f:
                data = line.split('\t')
                chrom = data[0]
                start = data[1]
                end = data[2]
                value = data[3].strip()
                if counter1 == 0:
                    prev_chrom = chrom
                    prev_value = value
                    block_start = start
                    prev_end = end
                else:
                    if prev_chrom != chrom:
                        new_line = [prev_chrom, block_start, prev_end, prev_value+'\n']
                        new_line = '\t'.join(new_line)
                        fout.write(new_line)
                        counter2 += 1
                        
                        prev_chrom = chrom
                        prev_value = value
                        prev_end = end
                        block_start = start
                    elif prev_chrom == chrom and prev_value == value:
                        pass
                    elif prev_chrom == chrom and prev_value != value:
                        new_line = [chrom, block_start, end, prev_value+'\n']
                        new_line = '\t'.join(new_line)
                        fout.write(new_line)
                        counter2 += 1
                        
                        prev_value = value
                        block_start = start
                        prev_end = end

                counter1 += 1
    print counter1
    print counter2
    
    os.remove(bedgraph)
    os.rename(bedgraph+'.tmp', bedgraph)

def bedgraph_reader(bedgraph, chromosomes=None):
    df = pd.read_csv(bedgraph, sep='\t', header=None, names=['chromosome','start','end','RPM'])
    
    if chromosomes is not None:
        df = df[df['chromosome'].isin(chromosomes)]
        
    df.index = df['chromosome'].str.cat(df['start'].apply(str),sep=':')
    
    return df



def write_bedgraph(dataframe, name):
    dataframe.to_csv(name, index=False, header=False, sep='\t')

def bedgraph_reader2(bedgraph, chromosomes=None, write=False):
    ''' For expanded begraphs (made with -a flag)'''
    df = pd.read_csv(bedgraph, sep='\t', header=None, names=['chromosome','start','RPM'])
    df.loc[:,'end'] = df['start']+1
    df = df[['chromosome','start','end','RPM']]
    
    if chromosomes is not None:
        df = df[df['chromosome'].isin(chromosomes)]
        
    df.index = df['chromosome'].str.cat(df['start'].apply(str),sep=':')
    if write:
        write_bedgraph(df, bedgraph)
    
    return df
    
def combine_stranded_bedgraph(directory, file_provided=False):
    bg_pairs = []
    if not file_provided:
        for file in os.listdir(directory):
            if file.endswith('plus.bedgraph'):
                bg_pairs.append((file,file.split('plus.bedgraph')[0]+'minus.bedgraph'))
    else:
        name1 = directory.split('.bam')[0]+'plus.bedgraph'
        name2 = directory.split('.bam')[0]+'minus.bedgraph'
        bg_pairs.append((name1, name2))

    for pair in bg_pairs:
        name = pair[0].split('.bedgraph')[0].split('plus')[0]
        if not name.endswith('_'):
            name = name+'_'
        plus = bedgraph_reader(pair[0])
        minus = bedgraph_reader(pair[1])
        minus[3] = minus[3].multiply(-1)

        new = plus.append(minus)
        new = new.sort_values([0,1])
        write_bedgraph(new, name)
##
 
def normalize_bedgraph(tagged, untagged, smooth=False, last=False):
    tagged_RPM = bedgraph_reader2(tagged, write=True)
    untagged_RPM = bedgraph_reader2(untagged, write=last)

    total = tagged_RPM.merge(untagged_RPM, right_index=True, left_index=True, how='left')
    total.loc[:,'norm RPM'] = total['RPM_x']/total['RPM_y']
    
    normalized = total[['chromosome_x','start_x','end_x','norm RPM']]
    normalized = normalized.replace([np.inf,np.inf*-1],np.NaN).dropna(how='any')
    
    normalized.to_csv(tagged.split('.bedgraph')[0]+'_norm.bedgraph', sep='\t', index=False, header=False)
    
    if smooth is False:
        collapse_bedgraph(tagged.split('.bedgraph')[0]+'_norm.bedgraph')
    
def smooth_bedgraphs(bedgraph_list, window):
    for bedgraph in bedgraph_list:
        bg_df = bedgraph_reader(bedgraph)
        new_bg = pd.DataFrame(columns=bg_df.columns)
        for chrom in set(bg_df['chromosome']):
            chrom_df = bg_df[bg_df['chromosome'] == chrom]
            chrom_df = chrom_df.sort_values(['start'])
            new_intensities = chrom_df['RPM'].rolling(window=window, center=True).mean()
            chrom_df.loc[:,'RPM'] = new_intensities
            new_bg = new_bg.append(chrom_df.dropna(how='any'))
            
        new_bg = new_bg.sort_values(['chromosome','start'])
        new_bg.loc[:,'start'] = new_bg['start'].apply(int)
        new_bg.loc[:,'end'] = new_bg['end'].apply(int)
        new_bg.to_csv(bedgraph.split('.bedgraph')[0]+'_{0}bp_smooth.bedgraph'.format(str(window)), sep='\t', index=False, header=False)
        
        collapse_bedgraph(bedgraph.split('.bedgraph')[0]+'_{0}bp_smooth.bedgraph'.format(str(window)))
        
def background_subtraction(bedgraph):
    bg_df = bedgraph_reader2(bedgraph)
    new_bg = pd.DataFrame(columns=bg_df.columns)
    
    print "Calculating background..."
    for chrom in set(bg_df['chromosome']):
        chrom_df = bg_df[bg_df['chromosome'] == chrom]
        chrom_df = chrom_df.sort_values(['start'])
        
        new_intensities = chrom_df['RPM'].rolling(window=5000, center=True).mean()
        chrom_df.loc[:,'RPM'] = new_intensities
        chrom_df.iloc[:2500,2] = chrom_df.iloc[2500,2]
        chrom_df.iloc[-2499:,2] = chrom_df.iloc[-2500,2]
        new_bg = new_bg.append(chrom_df.dropna(how='any'))
    
    print "Adjusting values..."
    bg_df.loc[:,'RPM'] = bg_df['RPM'] - new_bg['RPM']
    bg_df.loc[:,'end'] = bg_df['start'] + 1
    
    final_bg = bg_df[['chromosome','start','end','RPM']]
    neg_index = final_bg[final_bg['RPM'] < 0].index
    final_bg.loc[neg_index,['RPM']] = 0
    
    final_bg.to_csv(bedgraph.split('.bedgraph')[0]+'_sub.bedgraph', sep='\t', header=False, index=False)
    
    print "Collapsing bedgraph..."
    collapse_bedgraph(bedgraph.split('.bedgraph')[0]+'_sub.bedgraph')