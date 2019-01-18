import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.append('/home/jordan/RNA-is-awesome/')
import GeneTools as GT
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/CHIP/')
sys.path.append('/home/jordan/RNA-is-awesome/CHIP/')
import Compare_RPKM
import pandas as pd
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import subprocess
import os
import pysam
import json
from collections import OrderedDict
import seaborn as sns
from copy import deepcopy
from multiprocessing import Pool
import pickle

def run(cmd, logfile):
    '''Function to open subprocess, wait until it finishes and write all output to the logfile'''
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile, stderr=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

def align_fastq_chip(directory, threads=1, organism=None, adaptor='GATCGGAAGA', bowtie_ix=None):
    '''Automatically aligns all fastq.gz files in a drectory using Bowtie'''
    if directory[-1] != '/':
        directory = directory+'/'
        
    if organism is None:
        bowtie_ix=bowtie_ix
    elif 'crypto' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/Crypto_for_gobs'
    elif 'cerev' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/S288C/S288C'
        if 'DBK' in organism:
            bowtie_ix = '/home/jordan/GENOMES/S288C/S288C_DBK'
    elif 'pombe' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/POMBE/Spombe'
    elif 'candida' in organism.lower() or 'albicans' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/C_albicans'
    else:
        print "Organism not recognized"
        return None
        
    fastq_list = []
    for file in os.listdir(directory):
        if file.lower().endswith("fastq.gz") or file.lower().endswith("fastq"):
            fastq_list.append(directory+file)
    
    if not os.path.exists(directory+'MULTI/'):
        os.makedirs(directory+'MULTI/')
    if not os.path.exists(directory+'UNIQUE/'):
        os.makedirs(directory+'UNIQUE/')
    
    ##Check to make sure tophat hasn't already completed on these files. Do some cleanup
    for fastq in fastq_list:
        bam = None
        sorted_bam = None
        sam = None
        for align_type in ('MULTI/','UNIQUE/'):
            for name in os.listdir(directory+align_type):    
                if fastq.split('.fastq')[0] in name and name.endswith('_sorted.bam'):
                    fastq_list.remove(directory+fastq)
                    print "Bowtie and Samtools view already completed on "+fastq
                    sorted_bam = name
                elif fastq.split('.fastq')[0] in name and name.endswith('.sam'):
                    fastq_list.remove(directory+fastq)
                    print "Bowtie already completed on "+fastq
                    sam = name
                elif fastq.split('.fastq')[0] in name and name.endswith('.bam'):
                    fastq_list.remove(directory+fastq)
                    print "Bowtie already completed on "+fastq
                    bam = name
            if sorted_bam is not None:
                if sam is not None:
                    os.remove(directory+align_type+sam)
                if bam is not None:
                    os.remove(directory+align_type+bam)
            elif sorted_bam is None and bam is not None:
                if sam is not None:
                    os.remove(directory+align_type+sam)
    
    if len(fastq_list) > 0:
        processes = []
        for n, fastq in enumerate(fastq_list):
            prefix = fastq.split('/')[-1].split('.fastq')[0]
            print prefix
            with open(prefix+'.log', 'w') as logfile:
                trim_name = prefix+'_trim.fastq'
                if trim_name not in os.listdir(directory):
                    print "Trimming reads..."
                    cutadapt_args = "cutadapt -a {0} -m 20 -o {1}_trim.fastq {2}".format(adaptor, prefix, fastq)
                    logfile.write('**********Cutadapt output**********\n')
                    processes.append(subprocess.Popen(cutadapt_args, shell=True, universal_newlines=True, stdout=logfile, stderr=logfile))
                    if (n+1) % threads == 0 or n+1 == len(fastq_list):
                        for p in processes:
                            p.wait()
                    
                    #ret_code = run(cutadapt_args, logfile)

        for n, fastq in enumerate(fastq_list):
            prefix = fastq.split('/')[-1].split('.fastq')[0]
            print prefix
            trim_name = prefix+'_trim.fastq'
            with open(prefix+'.log', 'a') as logfile:
                if len([x for x in os.listdir(directory+'MULTI/') if prefix in x]) == 0:
                    logfile.write('\n**********Bowtie multi alignment output**********\n')
                    print "Aligning allowing multiple alignments (random assignment)..."
                    
                    multi_args = "bowtie -p{0} -v2 -M1 --best --chunkmbs 1024 --un {1}_multi_un.fastq --max {1}_multi.fastq {2} -q {3} --sam {1}_multi.sam".format(str(threads), directory+'MULTI/'+prefix, bowtie_ix, trim_name)
                    ret_code = run(multi_args, logfile)
                
                if len([x for x in os.listdir(directory+'UNIQUE/') if prefix in x]) == 0:
                    logfile.write('\n**********Bowtie unique alignment output**********\n')
                    print "Aligning allowing only unique alignments...\n"
                    
                    unique_args = "bowtie -p{0} -v2 -m1 --chunkmbs 1024 --un {1}_unique_un.fastq {2} -q {3} --sam {1}_unique.sam".format(str(threads), directory+'UNIQUE/'+prefix, bowtie_ix, trim_name)
                    ret_code = run(unique_args, logfile)
    
    for subdir in ('MULTI/','UNIQUE/'):
        sam_files = [x for x in os.listdir(directory+subdir) if x.endswith('.sam')]
        processes = []
        for n, sam in enumerate(sam_files):
            name = sam.split('.sam')[0]
            args = "samtools view -Sbo {0}.bam {0}.sam".format(directory+subdir+name)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            
            if (n+1) % threads == 0 or n+1 == len(sam_files):
                for p in processes:
                    p.wait()
                    
        bam_files = [x for x in os.listdir(directory+subdir) if x.endswith('bam') and not x.endswith('_sorted.bam')]
        processes = []
        for n, bam in enumerate(bam_files):
            name = bam.split('.bam')[0]
            args = "samtools sort {0}.bam -o {0}_sorted.bam".format(directory+subdir+name)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            if (n+1) % threads == 0 or n+1 == len(bam_files):
                for p in processes:
                    p.wait()
    
        sorted_bams = [x for x in os.listdir(directory+subdir) if x.endswith('_sorted.bam')]
        processes = []
        for n, bam in enumerate(sorted_bams):
            print bam
            args = "samtools index {0}".format(directory+subdir+bam)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            if (n+1) % threads == 1 or n+1 == len(sorted_bams):
                for p in processes:
                    p.wait()

        # Remove processing intermediates
        #for sam in [x for x in os.listdir(directory+subdir) if x.endswith('.sam')]:
        #    os.remove(directory+subdir+sam)
        for bam in [x for x in os.listdir(directory+subdir) if x.endswith('.bam')]:
            if not bam.endswith('_sorted.bam'):
                os.remove(directory+subdir+bam)
    
    for trim in [x for x in os.listdir(directory) if x.endswith('trim.fastq')]:
        os.remove(directory+trim)

def main():
    alt_adaptor = None
    index = None
    organism = None
    for n, arg in enumerate(sys.argv):
        if arg == '-h' or arg == '--help':
            print "\nUsage:\npython ChIP_tools.py --directory fastq_directory --threads num_threads --organism crypto/pombe/cerevisiae/candida <--adaptor GATCGGAAGA> <--index bowtie_index_prefix>"
            print "Note: --adaptor argument is optional and will default to the one shown above\n --index must  be provided if not using the default files for your organism."
            return None
        
        elif arg == '--directory':
            directory = sys.argv[n+1]
        elif arg == '--threads':
            threads = int(sys.argv[n+1])
        elif arg == '--organism':
            organism = sys.argv[n+1]
        elif arg == '--adaptor':
            alt_adaptor = sys.argv[n+1]
        elif arg == '--index':
            index = sys.argv[n+1]
    
    if alt_adaptor is not None:
        adaptor = alt_adaptor
    else:
        adaptor = 'GATCGGAAGA'
        
    align_fastq_chip(directory, threads=threads, organism=organism, adaptor=adaptor, bowtie_ix=index)
    
if __name__ == "__main__":
    main()
        

#########################################################################
## Functions for making scatter plots of the centermeres and telomeres ##
## (Used in Sandra's paper)                                            ##
#########################################################################
        
def make_name(string):
    name = string.split('/')[-1].split('_sorted')[0]
    return name

def prep_bam(df, bam, tx_dict):
    name = bam.split('/')[-1].split('_sorted')[0]
    read_list = Compare_RPKM.count_reads_in_ChIP(tx_dict, bam)
    df[name] = read_list
    return df

def plot_min_max(lists, ax):
    all_min = 0.9*min([x for y in lists for x in y])
    all_max = 1.1*max([x for y in lists for x in y])
    ax.set_ylim([all_min,all_max])
    ax.set_xlim([all_min,all_max])
    ax.plot([all_min,all_max],[all_min,all_max],color='black')
    return ax

def ChIP_rpkm_scatter(WCE_bam, WT1_bam, WT2_bam, Mut1_bam, Mut2_bam, gff3, plot_name, Z_change=False, cen_tel=False):
    '''Plots RPKM as scatter plots from two different samples - can do promoters or just centromeres and telomeres
    
    Parameters
    ----------
    WCE_bam : str or list
            Whole cell extract bam file, or list of WCE bam files - must be in the same order as the WT and Mut bam files.
            Example of list: ['WCE_WT1_sorted.bam','WCE_WT2_sorted.bam','WCE_Mut1_sorted.bam','WCE_Mut2_sorted.bam']
    WT1_bam : str
            Wild type (or first condition) bam file - replicate 1
    WT2_bam : str
            Wild type (or first condition) bam file - replicate 2
    Mut1_bam : str
            Mutant (or second condition) bam file - replicate 1
    Mut2_bam : str
            Mutant (or second condition) bam file - replicate 2
    gff3 : str
            gff3 file - if plotting centromeres and telomeres, provide a custom gff3 file with their locations.
            Otherwise provide the standard gff3 file for your organism
    plot_name : str
            name to save plot under (will be a .eps file)
    Z_change : bool, default `False`
            Whether or not to evaluate outliers that change in the mutant based on Z score
    cen_tel : bool, default `False`
            Whether to plot the RPKM for the promoter (1 kb upstream)+ORF (False) or for centromeres/telomeres (True)
            
    Outputs
    ------
    scatter plot : eps file'''
    
    tx_dict = GT.build_transcript_dict(gff3)
    if cen_tel is False:
        tx_dict = Compare_RPKM.make_promoter_dict(tx_dict, '/home/jordan/GENOMES/H99_chrom_lengths.json')
    
    df = pd.DataFrame(index=tx_dict.keys())
    if type(WCE_bam) == list:
        wce_list = WCE_bam
    else:
        wce_list = [WCE_bam]
    chip_list = [WT1_bam, WT2_bam, Mut1_bam, Mut2_bam]
    
    for bam in chip_list:
        df = prep_bam(df, bam, tx_dict)
    #print df.columns
    
    for bam in wce_list:
        df = prep_bam(df, bam, tx_dict)
    #print df.columns
        
    for column in df:
        df[column] = pd.to_numeric(df[column])
    
    names = df.columns
    if len(wce_list) == 4:
        df[names[0]+' Normalized'] = df[names[0]]/df[names[4]]
        df[names[1]+' Normalized'] = df[names[1]]/df[names[5]]
        df[names[2]+' Normalized'] = df[names[2]]/df[names[6]]
        df[names[3]+' Normalized'] = df[names[3]]/df[names[7]]
    elif len(wce_list) == 1:
        df[names[0]+' Normalized'] = df[names[0]]/df[names[4]]
        df[names[1]+' Normalized'] = df[names[1]]/df[names[4]]
        df[names[2]+' Normalized'] = df[names[2]]/df[names[4]]
        df[names[3]+' Normalized'] = df[names[3]]/df[names[4]]

    df['Enrichment 1'] = df[names[2]+' Normalized']/df[names[0]+' Normalized']
    df['Enrichment 2'] = df[names[3]+' Normalized']/df[names[1]+' Normalized']
    
    df = df[(df[names[0]] > 0) & (df[names[1]] > 0)]
    
    for_plot = []
    for_plot2 = []
    Z_scores = []
    for column in df:
        if 'Normalized' in column:
            for_plot2.append(column)
            df['log2 RPKM '+column.split(' ')[0]] = df[column].apply(np.log2)
            for_plot.append('log2 RPKM '+column.split(' ')[0])
        elif 'Enrichment' in column:
            df['Z-score '+column] = pd.Series(stats.mstats.zscore(df[column]), index=df.index)
            Z_scores.append('Z-score '+column)
    
    ## make the plot
    f, ax = plt.subplots(2, 2, figsize=(10,10))
    for axis in ax:
        axis[0] = plot_min_max([df[for_plot[0]].tolist()+df[for_plot[1]].tolist()+df[for_plot[2]].tolist()+df[for_plot[3]].tolist()], axis[0])
        axis[1] = plot_min_max([df[for_plot[0]].tolist()+df[for_plot[1]].tolist()+df[for_plot[2]].tolist()+df[for_plot[3]].tolist()], axis[1])
    
    ax[0,0].set_xlabel(for_plot[0])
    ax[0,0].set_ylabel(for_plot[1])
    ax[0,0].plot(df[for_plot[0]],df[for_plot[1]],'o', alpha=0.8, color='0.4')
    print stats.pearsonr(df[for_plot[0]],df[for_plot[1]])
    
    ax[0,1].set_xlabel(for_plot[2])
    ax[0,1].set_ylabel(for_plot[3])
    ax[0,1].plot(df[for_plot[2]],df[for_plot[3]],'o', alpha=0.8, color='0.4')
    print stats.pearsonr(df[for_plot[2]],df[for_plot[3]])
    
    ax[1,0].set_xlabel(for_plot[0])
    ax[1,0].set_ylabel(for_plot[2])
    
    if Z_change is True:
        ax[1,0].plot(df[for_plot[0]],df[for_plot[2]],'o', alpha=0.5, color='0.4')
        dep_df1 = df[(df[Z_scores[0]] >= 2) | (df[Z_scores[0]] <= -2)]
        ax[1,0].plot(dep_df1[for_plot[0]],dep_df1[for_plot[2]],'o', alpha=0.6, color='darkorchid')
    
    ax[1,1].set_xlabel(for_plot[1])
    ax[1,1].set_ylabel(for_plot[3])
    
    if Z_change is True:
        ax[1,1].plot(df[for_plot[1]],df[for_plot[3]],'o', alpha=0.5, color='0.4')
        dep_df2 = df[(df[Z_scores[1]] >= 2) | (df[Z_scores[1]] <= -2)]
        ax[1,1].plot(dep_df2[for_plot[1]],dep_df2[for_plot[3]],'o', alpha=0.6, color='darkorchid')
    
        up = set(df.index)
        down = set(df.index)
        for Z_col in Z_scores:
            up = up.intersection(df[df[Z_col] >= 2].index)
            down = down.intersection(df[df[Z_col] <= -2].index)
        df['Z-change'] = 'Unchanged'
        df.loc[up, 'Z-change'] = 'Up'
        df.loc[down, 'Z-change'] = 'Down'
    
    df.to_csv(plot_name+'.csv')    

    if cen_tel is True:
        cen_df = df[df.index.str.contains('Cen')]
        tel_df = df[df.index.str.contains('tel')]
        ax[1,0].plot(cen_df[for_plot[0]], cen_df[for_plot[2]], 'o', alpha=0.8, color='crimson', label='Centromeres')
        print stats.pearsonr(cen_df[for_plot[0]], cen_df[for_plot[2]])
        ax[1,0].plot(tel_df[for_plot[0]], tel_df[for_plot[2]], 'o', alpha=0.8, color='mediumblue', label='Telomeres')
        ax[1,0].legend()
        
        ax[1,1].plot(cen_df[for_plot[1]], cen_df[for_plot[3]], 'o', alpha=0.8, color='crimson', label='Centromeres')
        print stats.pearsonr(cen_df[for_plot[1]], cen_df[for_plot[3]])
        ax[1,1].plot(tel_df[for_plot[1]], tel_df[for_plot[3]], 'o', alpha=0.8, color='mediumblue', label='Telomeres')
        ax[1,1].legend()

    f.savefig(plot_name+'.eps', format='eps')
    print "Saved as "+plot_name+'.eps'
    
    plt.show()
    plt.clf()
    
    if cen_tel is True:
        cen_df = cen_df.sort_index()
        tel_df = tel_df.sort_index()
        
        cen_df['WT avg'] = (cen_df[for_plot2[0]] + cen_df[for_plot2[1]]).divide(2)
        cen_df['WT range'] = (cen_df[for_plot2[0]] - cen_df[for_plot2[1]]).apply(abs).divide(2)
        cen_df['Mut avg'] = (cen_df[for_plot2[2]] + cen_df[for_plot2[3]]).divide(2)
        cen_df['Mut range'] = (cen_df[for_plot2[2]] - cen_df[for_plot2[3]]).apply(abs).divide(2)
        
        tel_df['WT avg'] = (tel_df[for_plot2[0]] + tel_df[for_plot2[1]]).divide(2)
        tel_df['WT range'] = (tel_df[for_plot2[0]] - tel_df[for_plot2[1]]).apply(abs).divide(2)
        tel_df['Mut avg'] = (tel_df[for_plot2[2]] + tel_df[for_plot2[3]]).divide(2)
        tel_df['Mut range'] = (tel_df[for_plot2[2]] - tel_df[for_plot2[3]]).apply(abs).divide(2)
        
        f2, ax2 = plt.subplots(2, 1, figsize=(10,10))
        width = 0.35
        
        ind_cen = np.arange(len(cen_df))
        ax2[0].bar(ind_cen, cen_df['WT avg'], width, color='mediumblue', yerr=cen_df['WT range'], label="WT")
        ax2[0].bar(ind_cen + width, cen_df['Mut avg'], width, color='crimson', yerr=cen_df['Mut range'], label="Mutant")
        
        ax2[0].set_ylabel('RPKM')
        ax2[0].set_title('Centromeres')
        ax2[0].set_xticks(ind_cen + width / 2)
        ax2[0].set_xticklabels(cen_df.index)
        ax2[0].legend()
        
        ind_tel = np.arange(len(tel_df))
        ax2[1].bar(ind_tel, tel_df['WT avg'], width, color='mediumblue', yerr=tel_df['WT range'], label="WT" )
        ax2[1].bar(ind_tel + width, tel_df['Mut avg'], width, color='crimson', yerr=tel_df['Mut range'], label="Mutant")
        
        ax2[1].set_ylabel('RPKM')
        ax2[1].set_title('Telomeres')
        ax2[1].set_xticks(ind_tel + width / 2)
        ax2[1].set_xticklabels(tel_df.index, rotation='vertical')
        ax2[1].legend()

        f2.savefig(plot_name+'_bar.eps', format='eps')
        print "Saved as "+plot_name+'_bar.eps'

        plt.show()
        plt.clf()
        
        
##################################################
## Functions for analyzing peaks called by MACS ##
##################################################

def add_transcript(df, gff3, organism=None):
    ''' Used by MACS_peak_RPKM_scatters'''
    tx_dict = GT.build_transcript_dict(gff3, organism=organism)
    transcripts = set()
    transcripts_for_df = []
    for ix, r in df.iterrows():
        r_tx = ''
        chrom_tx = {k:v for k, v in tx_dict.items() if v[3] == r['chr']}
        for tx, info in chrom_tx.iteritems():
            if tx[:-2] not in r_tx:
                if info[2] == '+':
                    if r['abs_summit'] in range(info[0]-500,info[1]):
                        r_tx = r_tx+tx[:-2]+','
                        transcripts.add(tx[:-2])
                if info[2] == '-':
                    if r['abs_summit'] in range(info[0],info[1]+500):
                        r_tx = r_tx+tx[:-2]+','
                        transcripts.add(tx[:-2])
        if r_tx == '':
            for tx, info in chrom_tx.iteritems():
                if tx[:-2] not in r_tx:
                    if info[2] == '+':
                        if r['abs_summit'] in range(info[0]-1000,info[1]):
                            r_tx = r_tx+tx[:-2]+','
                            transcripts.add(tx[:-2])
                    if info[2] == '-':
                        if r['abs_summit'] in range(info[0],info[1]+1000):
                            r_tx = r_tx+tx[:-2]+','
                            transcripts.add(tx[:-2])
        transcripts_for_df.append(r_tx)
    df.loc[:,'transcript'] = transcripts_for_df
    return df, transcripts
        
def compare_MACS_output(rep1_xls, rep2_xls, untagged_xls, organism, return_df=False, min_overlap=0.5, cutoff=2):
    ''' Used by MACS_peak_RPKM_scatters'''
    
    if 'crypto' in organism.lower():
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
        organism = None
    elif 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        organism = 'pombe'
    elif 'cerev' in organism.lower() or '288' in organism:
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        organism = None
        
    df1 = pd.read_csv(rep1_xls, sep='\t', skiprows=28, header=0)
    df1 = df1[df1['fold_enrichment'] >= cutoff]
    df2 = pd.read_csv(rep2_xls, sep='\t', skiprows=28, header=0)
    df2 = df2[df2['fold_enrichment'] >= cutoff]
    df_un = pd.read_csv(untagged_xls, sep='\t', skiprows=28, header=0)
    
    # Determine if peak is in each replicate
    rep = []
    enrich2 = []
    starts = []
    ends = []
    for ix, r in df1.iterrows():
        r_rep = False
        peak = range(r['start'],r['end'])
        
        # Check if there are peaks in the second replicate with the minimum overlap
        rep2_df = df2[df2['chr'] == r['chr']]
        rep2_matches = []
        for ix2, r2 in rep2_df.iterrows():
            peak_overlap = len(set(peak).intersection(range(r2['start'],r2['end'])))
            if peak_overlap/float(len(peak)) > min_overlap:
                rep2_matches.append(ix2)
        
        
        rep2_matches = rep2_df[rep2_df.index.isin(rep2_matches)]
        if len(rep2_matches) > 0:
            rep.append(True)
            r_rep = True
            
            # Determine the nearest absolute summit and record fold enrichment in replicate
            s = (rep2_matches['abs_summit']-r['abs_summit']).apply(abs)
            min_dist = s[s == min(s)].index[0]
            enrich2.append(rep2_matches.loc[min_dist,'fold_enrichment'])
            
            # Contract peak window to the reproducible region
            new_start = max(rep2_matches.loc[min_dist,'start'], r['start'])
            new_end = min(rep2_matches.loc[min_dist,'end'], r['end'])
            starts.append(new_start)
            ends.append(new_end)
        else:
            rep.append(False)
            enrich2.append(np.NaN)
            starts.append(r['start'])
            ends.append(r['end'])
            
    df1.loc[:,'In replicate'] = rep
    df1.loc[:,'fold_enrichment2'] = enrich2
    df1.loc[:,'start'] = starts
    df1.loc[:,'end'] = ends
    
    # Determine if peak is in untagged
    un = []
    for ix, r in df1.iterrows():
        peak = range(r['start'],r['end'])
        un_matches = df_un[(df_un['chr'] == r['chr']) & (df_un['abs_summit'].isin(peak))]
        if len(un_matches) > 0:
            un.append(True)
        else:
            un.append(False)
    df1.loc[:,'In untagged'] = un
    
    # Filter based on reproducibility and untagged
    print "\nPeaks in "+rep1_xls+":"
    print len(df1)
    print "Peaks in "+rep2_xls+":"
    print len(df2)
    filtered = df1[df1['In replicate'] == True]
    print "Peaks in both replicates:"
    print len(filtered)
    filtered = filtered[filtered['In untagged'] == False]
    print "Peaks not in untagged:"
    print len(filtered)
    filtered = filtered.drop(['In replicate','In untagged'], axis=1)
    
    # Remove redundant peaks
    filtered = filtered.sort_values(['fold_enrichment'], ascending=False)
    filtered = filtered.drop_duplicates(subset=['chr','start','end'])
    print "Peaks remaining after removing summits within the same peak region:"
    print len(filtered)
    
    # Add transcripts
    filtered, transcripts = add_transcript(filtered, gff3, organism=organism)
    
    filtered.to_csv(rep1_xls.split('/')[-1].split('.xls')[0]+'_comparison.csv', index=False)
    
    with open(rep1_xls.split('/')[-1].split('.xls')[0]+'_genes_with_peaks.txt', 'w') as fout:
        for transcript in transcripts:
            fout.write(transcript+'\n')
    if return_df:
        return filtered
    else:
        return None

def wt_v_mut_MACS(df1, df2, min_overlap=0.5):
    ''' Used by MACS_peak_RPKM_scatters'''
    enrich_mut1 = []
    enrich_mut2 = []
    for ix, r in df1.iterrows():
        match = False
        chrom_df = df2[df2['chr'] == r['chr']]
        range1 = set(range(r['start'],r['end']))
        for ix2, r2 in chrom_df.iterrows():
            range2 = range(r2['start'],r2['end'])
            if len(range1.intersection(range2))/float(len(range1)) >= min_overlap:
                enrich_mut1.append(r2['fold_enrichment'])
                enrich_mut2.append(r2['fold_enrichment2'])
                match = True
                break
        if not match:
            enrich_mut1.append(np.NaN)
            enrich_mut2.append(np.NaN)
    
    df1.loc[:,'fold_enrichment mut1'] = enrich_mut1
    df1.loc[:,'fold_enrichment mut2'] = enrich_mut2
    
    print len(df1)
    counter = 0
    for ix, r in df2.iterrows():
        match = False
        chrom_df = df1[df1['chr'] == r['chr']]
        range1 = set(range(r['start'],r['end']))
        for ix2 in chrom_df.index:
            range2 = range(chrom_df.loc[ix2,'start'], chrom_df.loc[ix2,'end'])
            if len(range1.intersection(range2))/float(len(range1)) >= min_overlap:
                match = True
                break
        if not match:
            new_row = df2.loc[ix,:]
            new_row = new_row.rename({'fold_enrichment': 'fold_enrichment mut1', 'fold_enrichment2': 'fold_enrichment mut2'})
            df1 = df1.append(new_row)
            df1 = df1.reset_index(drop=True)
            counter += 1
    print counter
    print len(df1)
    
    df1.index = df1['chr'].str.cat(df1['start'].apply(str),sep=':').str.cat(df1['end'].apply(str),sep='-')
            
    return df1    
    
def MACS_scatter_plots(df, xls_pair1, xls_pair2):
    ''' Used by MACS_peak_RPKM_scatters'''
    WT_a = [x for x in df.columns if x in xls_pair1[0]][0]
    WT_b = [x for x in df.columns if x in xls_pair1[1]][0]
    MUT_a = [x for x in df.columns if x in xls_pair2[0]][0]
    MUT_b = [x for x in df.columns if x in xls_pair2[1]][0]
    
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10,10))
    ax[0][0].scatter(df[WT_a].apply(np.log2), df[WT_b].apply(np.log2), s=15, color='0.3', alpha=0.5)
    ax[0][1].scatter(df[MUT_a].apply(np.log2), df[MUT_b].apply(np.log2), s=15, color='0.3', alpha=0.5)
    ax[1][0].scatter(df[WT_a].apply(np.log2), df[MUT_a].apply(np.log2), s=15, color='0.3', alpha=0.5)
    ax[1][1].scatter(df[WT_b].apply(np.log2), df[MUT_b].apply(np.log2), s=15, color='0.3', alpha=0.5)

    min_xy = min(ax[0][0].get_xlim()[0], ax[0][0].get_ylim()[0])
    max_xy = max(ax[0][0].get_xlim()[1], ax[0][0].get_ylim()[1])
    
    for subax in ax:
        for sub in subax:
            min_xy = min(min_xy, sub.get_xlim()[0], sub.get_ylim()[0])
            max_xy = max(max_xy, sub.get_xlim()[1], sub.get_ylim()[1])

    for subax in ax:
        for sub in subax:
            sub.set_xlim(min_xy,max_xy)
            sub.set_ylim(min_xy,max_xy)

    ax[0][0].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[0][1].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[1][0].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[1][1].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)

    ax[0][0].set_xlabel(xls_pair1[0].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)
    ax[0][0].set_ylabel(xls_pair1[1].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)

    ax[0][1].set_xlabel(xls_pair2[0].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)
    ax[0][1].set_ylabel(xls_pair2[1].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)

    ax[1][0].set_xlabel(xls_pair1[0].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)
    ax[1][0].set_ylabel(xls_pair2[0].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)

    ax[1][1].set_xlabel(xls_pair1[1].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)
    ax[1][1].set_ylabel(xls_pair2[1].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)

    fig.tight_layout()
    plt.show()
    plt.clf()
    return fig

def genome_tiles(chromosome_sizes, tile_size=1000):
    ''' Used by MACS_peak_RPKM_scatters'''
    names = []
    chroms = []
    start = []
    end = []
    with open(chromosome_sizes) as f:
        for line in f:
            chrom = line.split('\t')[0]
            size = int(line.split('\t')[1].strip())
            n = 0
            while n*tile_size < size-tile_size:
                chroms.append(chrom)
                names.append(chrom+'-'+str(n))
                start.append(n*tile_size)
                end.append(n*tile_size+tile_size)
                n += 1
            chroms.append(chrom)
            names.append(chrom+'-'+str(n))
            start.append(n*tile_size)
            end.append(size)
    df = pd.DataFrame(index=names)
    df.loc[:,'chromosome'] = chroms
    df.loc[:,'start'] = start
    df.loc[:,'end'] = end
    return df

def count_reads_in_tiles(bam_list, WCE, chromosome_sizes, tile_size=1000):
    ''' Used by MACS_peak_RPKM_scatters'''
    df = genome_tiles(chromosome_sizes)
    
    bam_list.append(WCE)
    names = []
    for bam in bam_list:
        bam_name = bam.split('/')[-1].split('_sorted.bam')[0]
        if bam == WCE:
            WCE_name = bam_name
        else:
            names.append(bam_name)
        open_bam = pysam.Samfile(bam)

        read_dens = []
        for ix, r in df.iterrows():
            read_counter = 0
            try:
                for read in open_bam.fetch(r['chromosome'],r['start'],r['end']):
                    read_counter += 1
                read_dens.append(read_counter/((r['end']-r['start'])/1000.))
            except ValueError:
                read_dens.append(np.NaN)
        
        df.loc[:,bam_name] = read_dens
        df.loc[:,bam_name] = df[bam_name].divide(sum(df[bam_name].dropna())/1000000.)
    for name in names:
        df.loc[:,name+' norm to WCE'] = df[name]/df[WCE_name]
    
    df = df.dropna(how='any')
    return df


def MACS_peak_RPKM_scatters(xls_pair1, xls_pair2, untagged_xls, bam_list, WCE_bam_list, name, organism='crypto', enrichment_cutoff=2, min_overlap=0.5, exclude_chrom=None, collapse_peaks=True, gene_list=None, adjust=True):
    
    '''This function does several things:
        1: It compares the MACS output from replicate ChIP-seq experiments and defines a minimal set of reproducible peaks.
               Importantly, it gives the peak summit as the largest peak in the region and trims the peak boundaries to the
               minimum width that is reproducible in both replicates
        2: It looks for  peaks that appear in the untagged sample and removes them.
        3: It compares mutant to wild type and finds which peaks are present on both or just one.
        4: It counts reads in the newly defined peak set to quantitatively compare between samples (and normalizes to WCE)
        5: It adjusts for amplification or efficiency bias by performing a linear fit to the background
        6: It assigns transcripts based on whether the peak is within or upstream of an ORF - tries within 500 bp first then
           expands search to 1 kb
        
    
    Parameters
    ----------
    xls_pair1 : tuple of str
            Wild type or condition 1 replicates - .xls output from MACS
    xls_pair2 : tuple of str
            Mutant or condition 2 replicates - .xls output from MACS
    untagged_xls : str
            Untagged or no antibody .xls output from MACS
    bam_list : list of str
            The bam files to use for re-quantitating the peak enrichments
    WCE_bam_list : list of str
            The corresponding whole cell extract bam files - make sure they match the length of the bam_list and correspond
            in position to the bam_list
    name : str
            Prefix for saving files
    organism : str, default 'crypto'
            change to 'pombe' or 'cerevisiae' if using with another organism
    min_overlap : float, default 0.5
            minimum overlap required for peaks to be considered reproducible (default is 50%)
    exclude_chrom : str or list, default `None`
            Chromosome(s) to exclude from analysis
    collapse_peaks : bool, default `True`
            Whether or not peaks are collapsed by transcript in the 'by_transcript.csv' output file
            
    Outputs
    ------
    Comparison csv files : csv files containing the reproducible, reduced set of MACS peaks
    Genes with peaks : files with a list of genes that had at least one reproducible peak
    Scatter spreadsheet : The data the corresponds to the scatter plots - RPKM of each peak normalized to whole cell extract
    Scatter plot : pdf file with scatter plot
    By transcript spreadsheet : Summation of enrichment values by transcript (if collapse_peaks is True)
    '''
    
    if type(WCE_bam_list) == list:
        single_WCE = False
        if len(bam_list) != len(WCE_bam_list):
            if len(WCE_bam_list) != 1:
                print "Must provide only 1 WCE bam file or a matched WCE bam file for each sample in the bam_list!"
                return None
    elif type(WCE_bam_list) == str:
        single_WCE = True
        WCE = WCE_bam_list
        WCE_bam_list = []
        for bam in bam_list:
            WCE_bam_list.append(WCE)
    else:
        print "Whole cell extract is neither a list or a string!"
        return None
            
    # First run the comparison script to get peaks that align between replicates and to contract peak regions
    df_wt = compare_MACS_output(xls_pair1[0], xls_pair1[1], untagged_xls, organism, return_df=True, cutoff=enrichment_cutoff)
    df_mut = compare_MACS_output(xls_pair2[0], xls_pair2[1], untagged_xls, organism, return_df=True, cutoff=enrichment_cutoff)
    
    # Now need to run through and compare peaks between genotypes - see code for MACS enrichment scatters
    compare_df_A = wt_v_mut_MACS(df_wt, df_mut, min_overlap=min_overlap)
    compare_df_A.to_csv(name+'_MACS_wt_mut_merge.csv')
    
    # Generate dictionary of peaks to quantitate
    peak_dict = {}
    for ix, r in compare_df_A.iterrows():
        key = r['chr']+':'+str(r['start'])+'-'+str(r['end'])
        peak_dict[key] = [r['start'],r['end'],'+',r['chr'], r['transcript']]
    
    if exclude_chrom is not None:
        if type(exclude_chrom) == list:
            for chrom in exclude_chrom:
                peak_dict = {k:v for k,v in peak_dict.items() if v[3] != chrom}
        elif type(exclude_chrom) == str:
            print exclude_chrom
            print len(peak_dict)
            peak_dict = {k:v for k,v in peak_dict.items() if v[3] != exclude_chrom}
            print len(peak_dict)
    
    # Then need to count reads in each bam file
    print "\nCalculating RPKM in peak regions..." 
    data_dict = {}
    bam_names = []
    for bam_file in bam_list: 
        bam_name = bam_file.split('_sorted')[0].split('/')[-1]
        bam_names.append(bam_name)
        print bam_name
        data_dict[bam_name] = Compare_RPKM.count_reads_in_ChIP(peak_dict, bam_file)
        new_ix = data_dict[bam_name].index

    # Assemble dataframe from results
    data_df_columns = data_dict.keys()
    data_df = pd.DataFrame(index=new_ix, columns=data_df_columns)
    for col in data_df.columns:
        data_df.loc[:,col] = data_dict[col]
    in_wt = []
    in_mut = []
    for ix, r in data_df.iterrows():
        if str(compare_df_A.loc[ix,'fold_enrichment']) == 'nan':
            in_wt.append(False)
        else:
            in_wt.append(True)
        
        if str(compare_df_A.loc[ix,'fold_enrichment mut1']) == 'nan':
            in_mut.append(False)
        else:
            in_mut.append(True)
    data_df.loc[:,'In WT'] = in_wt
    data_df.loc[:,'In Mut'] = in_mut

    # Now get the background level and apply adjustment
    if adjust:
        print "\n Adjusting for biases..."
        chromosome_sizes = '/home/jordan/GENOMES/crypto_for_bedgraph.genome'
        tile_df = count_reads_in_tiles(bam_list, WCE_bam_list[0], chromosome_sizes)
        
        slopes = []
        for bam_name in bam_names:
            peak_min=None
            if bam_name in xls_pair1[0] or bam_name in xls_pair1[1]:
                peak_min = np.percentile(data_df[data_df['In WT'] == True][bam_name], 0.5)
            elif bam_name in xls_pair2[0] or bam_name in xls_pair2[1]:
                peak_min = np.percentile(data_df[data_df['In Mut'] == True][bam_name], 0.5)
            if peak_min is None:
                print "BAM file name does not match xls file name!"
                print bam_name
                return None

            s = tile_df.sort_values(bam_name+' norm to WCE').reset_index()[bam_name+' norm to WCE']
            top = max(s[(s <= peak_min)].index)
            slope = stats.linregress(tile_df.sort_values(bam_name).reset_index().iloc[1000:top].index,
                tile_df.sort_values(bam_name).reset_index().iloc[1000:top][bam_name])[0]*1000
            print bam_name
            print slope
            
            slopes.append(slope)
            #data_df.loc[:,bam_name] = data_df[bam_name].divide(slope)
            
        adj_slopes = [x/min(slopes) for x in slopes]
        for n, bam_name in enumerate(bam_names):
            data_df.loc[:,bam_name] = data_df[bam_name].divide(adj_slopes[n])

    # Count whole cell extract and normalize
    if not single_WCE:
        for n, bam_file in enumerate(WCE_bam_list):
            print bam_file
            WCE = Compare_RPKM.count_reads_in_ChIP(peak_dict, bam_file)
            data_df.loc[:,bam_names[n]] = data_df[bam_names[n]]/WCE
    else:
        print bam_file
        WCE = Compare_RPKM.count_reads_in_ChIP(peak_dict, WCE_bam_list[0])
        for bam in bam_names:
            data_df.loc[:,bam] = data_df[bam]/WCE
            
    # Add transcripts to spreadsheet
    tx_list = []
    for ix, r in data_df.iterrows():
        tx_list.append(peak_dict[ix][4])
        
    data_df.loc[:,'transcript'] = tx_list
    data_df.to_csv(name+'.csv')
    
    # Make additional spreadsheet organized by gene
    data_df2 = deepcopy(data_df)
    data_df2 = data_df2.reset_index()
    data_df2 = data_df2.rename(columns={"index":"peak coord"})

    tx_df = pd.DataFrame(columns=data_df2.columns)
    for ix, r in data_df2.iterrows():
        try:
            if ',' in r['transcript']:
                tx_list = r['transcript'].split(',')
                for tx in tx_list:
                    if len(tx) > 0:
                        new_r = deepcopy(r)
                        new_r['transcript'] = tx
                        tx_df = tx_df.append(new_r)
            else:
                tx_df = tx_df.append(r)
        except TypeError:
            new_r = deepcopy(r)
            new_r['transcript'] = "None"
            tx_df = tx_df.append(new_r)

    if not collapse_peaks:
        tx_df = tx_df.reset_index(drop=True)
        tx_df.to_csv(name+'_by_tx.csv')
    else:
        coll_tx_df = pd.DataFrame(columns=tx_df.columns)
        groups = tx_df.columns.to_series().groupby(tx_df.dtypes).groups
        data_cols = [v for k, v in groups.items() if k == 'float64']
        for tx in set(tx_df['transcript']):
            if tx != 'None':
                mini_df = tx_df[tx_df['transcript'] == tx]
                if len(mini_df) > 1:
                    mini_df = mini_df.reset_index(drop=True)
                    new_r = deepcopy(mini_df.iloc[0,:])
                    for col in data_cols:
                        new_r[col] = mini_df[col].apply(sum)
                    coll_tx_df = coll_tx_df.append(new_r)
                else:
                    coll_tx_df = coll_tx_df.append(mini_df)
        coll_tx_df.index = coll_tx_df['transcript']
        coll_tx_df = coll_tx_df.drop(['peak coord','transcript','In WT','In Mut'], axis=1)
        coll_tx_df.to_csv(name+'_by_tx.csv')
    
    # Make scatter plots - 0v1, 2v3, 0v2, 1v3
    fig = MACS_scatter_plots(data_df, xls_pair1, xls_pair2)
    fig.savefig(name+'.pdf',format='pdf',bbox_inches='tight')
    
    #return data_df 
    
def get_transcripts_from_peak_csv(csv):
    df = pd.read_csv(csv, index_col=0)

    transcripts = []
    for ix, r in df.iterrows():
        if type(r['transcript']) == str:
            r_tx = r['transcript'].split(',')
            transcripts = transcripts + r_tx
    transcripts = set([x for x in transcripts if len(x) > 0])
    print len(transcripts)

    with open(csv.split('.csv')[0]+'_transcripts.txt','w') as fout:
        for tx in transcripts:
            fout.write(tx+'\n')
            
def MACS_Z_score(csv, wt, mut):
    data_df = pd.read_csv(csv, index_col=0)
    wt_names = [x for x in data_df.columns if wt in x]
    mut_names = [x for x in data_df.columns if mut in x]

    s1a = (data_df[mut_names[0]]/data_df[wt_names[0]]).apply(np.log2)
    s1a_Z = pd.Series(stats.mstats.zscore(s1a), index=data_df.index)
    
    s1b = (data_df[mut_names[0]]/data_df[wt_names[1]]).apply(np.log2)
    s1b_Z = pd.Series(stats.mstats.zscore(s1b), index=data_df.index)

    s2a = (data_df[mut_names[1]]/data_df[wt_names[1]]).apply(np.log2)
    s2a_Z = pd.Series(stats.mstats.zscore(s2a), index=data_df.index)
    
    s2b = (data_df[mut_names[1]]/data_df[wt_names[0]]).apply(np.log2)
    s2b_Z = pd.Series(stats.mstats.zscore(s2b), index=data_df.index)
    
    up = set(s1a_Z[s1a_Z >= 2].index)
    for Z in [s1b_Z,s2a_Z,s2b_Z]:
        new_up_ix = set(Z[Z >= 2].index)
        up = up.intersection(new_up_ix)
        
    down = set(s1a_Z[s1a_Z <= -2].index)
    for Z in [s1b_Z,s2a_Z,s2b_Z]:
        new_down_ix = set(Z[Z <= -2].index)
        down = down.intersection(new_down_ix)

    #up = set(s1_Z[s1_Z >= 1.8].index).intersection(s2_Z[s2_Z >= 1.8].index)
    #print "Genes up in mutant: "+str(len(up))
    #down = set(s1_Z[s1_Z <= -1.8].index).intersection(s2_Z[s2_Z <= -1.8].index)
    #print "Genes down in mutant: "+str(len(down))
    
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10,10))
    ax[0][0].scatter(data_df[wt_names[0]].apply(np.log2), data_df[wt_names[1]].apply(np.log2), s=15, color='0.3', alpha=0.5)
    ax[0][1].scatter(data_df[mut_names[0]].apply(np.log2), data_df[mut_names[1]].apply(np.log2), s=15, color='0.3', alpha=0.5)

    up_df = data_df[data_df.index.isin(up)]
    up_df.to_csv(csv.split('.csv')[0]+'_up_in_mut.csv')
    with open(csv.split('.csv')[0]+'_up.txt', 'w') as fout:
        for ix, r in up_df.iterrows():
            try:
                for tx in r['transcript'].rstrip(',').split(','):
                    fout.write(tx+'\n')
            except AttributeError:
                pass
    GT.crypto_annotate(csv.split('.csv')[0]+'_up.txt')
    os.remove(csv.split('.csv')[0]+'_up.txt')
    
    down_df = data_df[data_df.index.isin(down)]
    down_df.to_csv(csv.split('.csv')[0]+'_down_in_mut.csv')
    with open(csv.split('.csv')[0]+'_down.txt', 'w') as fout:
        for ix, r in down_df.iterrows():
            try:
                for tx in r['transcript'].rstrip(',').split(','):
                    fout.write(tx+'\n')
            except AttributeError:
                pass
    GT.crypto_annotate(csv.split('.csv')[0]+'_down.txt')
    os.remove(csv.split('.csv')[0]+'_down.txt')
    
    other_df = data_df[(~data_df.index.isin(up)) & (~data_df.index.isin(down))]

    ax[1][0].scatter(other_df[wt_names[0]].apply(np.log2), other_df[mut_names[0]].apply(np.log2), s=15, color='0.7', alpha=0.5)
    ax[1][0].scatter(up_df[wt_names[0]].apply(np.log2), up_df[mut_names[0]].apply(np.log2), s=15, color='orangered')
    ax[1][0].scatter(down_df[wt_names[0]].apply(np.log2), down_df[mut_names[0]].apply(np.log2), s=15, color='royalblue')

    ax[1][1].scatter(other_df[wt_names[1]].apply(np.log2), other_df[mut_names[1]].apply(np.log2), s=15, color='0.7', alpha=0.5)
    ax[1][1].scatter(up_df[wt_names[1]].apply(np.log2), up_df[mut_names[1]].apply(np.log2), s=15, color='orangered')
    ax[1][1].scatter(down_df[wt_names[1]].apply(np.log2), down_df[mut_names[1]].apply(np.log2), s=15, color='royalblue')

    min_xy = min(ax[0][0].get_xlim()[0], ax[0][0].get_ylim()[0])
    max_xy = max(ax[0][0].get_xlim()[1], ax[0][0].get_ylim()[1])

    for subax in ax:
        for sub in subax:
            min_xy = min(min_xy, sub.get_xlim()[0], sub.get_ylim()[0])
            max_xy = max(max_xy, sub.get_xlim()[1], sub.get_ylim()[1])

    for subax in ax:
        for sub in subax:
            sub.set_xlim(min_xy,max_xy)
            sub.set_ylim(min_xy,max_xy)

    ax[0][0].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[0][1].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[1][0].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[1][1].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)

    ax[0][0].set_xlabel(wt_names[0]+'\n Enrichment over WCE', fontsize=14)
    ax[0][0].set_ylabel(wt_names[1]+'\n Enrichment over WCE', fontsize=14)

    ax[0][1].set_xlabel(mut_names[0]+'\n Enrichment over WCE', fontsize=14)
    ax[0][1].set_ylabel(mut_names[1]+'\n Enrichment over WCE', fontsize=14)

    ax[1][0].set_xlabel(wt_names[0]+'\n Enrichment over WCE', fontsize=14)
    ax[1][0].set_ylabel(mut_names[0]+'\n Enrichment over WCE', fontsize=14)

    ax[1][1].set_xlabel(wt_names[1]+'\n Enrichment over WCE', fontsize=14)
    ax[1][1].set_ylabel(mut_names[1]+'\n Enrichment over WCE', fontsize=14)

    fig.tight_layout()
    fig.savefig(csv.split('.csv')[0]+'_Z_scatter.pdf', format='pdf', bbox_inches='tight')
    plt.show()
    plt.clf()
    
def RNAseq_ChIP_overlap_volcano(DESeq2_csv, ChIP_gene_list):
    RNA = GT.load_DESeq2_results([DESeq2_csv])
    chip = pd.read_csv(ChIP_gene_list, index_col=0, header=None)
    RNA = RNA[RNA.index.isin(chip.index)]
    GT.volcano_plot(RNA)

def get_peak_sequence2(csv, gene_list=None, fa_dict_loc='/home/jordan/GENOMES/H99_fa.json'):
    if type(gene_list) == str:
        holder = []
        with open(gene_list) as f:
            for line in f:
                holder.append(line.strip())
        gene_list = holder
    
    with open(fa_dict_loc) as f: fa_dict = json.load(f)
    df = pd.read_csv(csv, index_col=0)
    
    seq_list = []
    seq_list2 = []
    for ix, r in df.iterrows():
        chrom = ix.split(':')[0]
        start = int(ix.split(':')[1].split('-')[0])
        end = int(ix.split(':')[1].split('-')[1])
        seq = GT.seq_simple(chrom, start, end, '+', fa_dict)
        try:
            transcripts = [x for x in r['transcript'].split(',') if '_1' not in x]
        except AttributeError:
            transcripts = ['']
        if gene_list is not None:
            other = set(transcripts).difference(gene_list)
            transcripts = set(transcripts).intersection(gene_list)
            
            if len(transcripts) > 0:
                seq_list.append((','.join(transcripts), seq))
            else:
                seq_list2.append((','.join(other), seq))
        else:
            seq_list.append((','.join(transcripts), seq))
    
    print csv.split('/')[-1].split('.csv')[0]+'_sequences.fasta'
    with open(csv.split('/')[-1].split('.csv')[0]+'_sequences.fasta', 'w') as fout:
        for tx, seq in seq_list:
            fout.write('>'+tx+'\n')
            fout.write(seq+'\n')
            
    if gene_list is not None:
        with open(csv.split('/')[-1].split('.csv')[0]+'_other_sequences.fasta', 'w') as fout:
            for tx, seq in seq_list2:
                fout.write('>'+tx+'\n')
                fout.write(seq+'\n')
                
def plot_peaks(bam_list, peak_csv, organism='crypto', colors=None, save_dir=None, log_scale=False, same_yaxis=True, rpm=True):
    organism, gff3, fa_dict, bowtie_index = GT.find_organism_files(organism)
    tx_dict = GT.build_transcript_dict(gff3, organism=organism)
    fix_info = {'I':'chr1','II':'chr2','III':'chr3','chr1':'I','chr2':'II','chr4':'IV','chr5':'V','chr6':'VI',
                'chr7':'VII','chr8':'VIII','chr9':'IX','chr10':'X','chr11':'XI','chr12':'XII','chr13':'XIII',
                'chr14':'XIV','chr15':'XV','chr16':'XVI','-':'+','+':'-','chr1':'I','chr2':'II','chr3':'III'}
    if organism == 'pombe':
        tx_suffix = '.1'
    else:
        tx_suffix = 'T0'
    
    if colors is None:
        colors = []
        for n in range(len(bam_list)): colors.append('k')
    
    open_bams = {}
    totals = {}
    for bam in bam_list:
        open_bams[bam] = pysam.Samfile(bam)
        if rpm:
            totals[bam] = GT.count_aligned_reads(bam)
        else:
            totals[bam] = 1
    
    peak_df = pd.read_csv(peak_csv, index_col=0)
    num_ax = len(bam_list)+1
    
    for ix, r in peak_df.iterrows():
        fig, ax = plt.subplots(num_ax, figsize=(10,num_ax), sharex=True)
        fig.subplots_adjust(hspace=0)
        
        chrom = ix.split(':')[0]
        start = int(ix.split(':')[1].split('-')[0])-1000
        end = int(ix.split(':')[1].split('-')[1])+1000
        
        max_y = 0
        junc_ymax = 0
        for n, bam in enumerate(bam_list):
            try:
                bam_iter = open_bams[bam].fetch(chrom, start, end)
            except ValueError:
                chrom = fix_info[chrom]
                bam_iter = open_bams[bam].fetch(chrom, start, end)
            s = GT.generate_read_series(bam_iter, chrom, start, end, '+')
            s.add(GT.generate_read_series(bam_iter, chrom, start, end, '-'))
            s = s.replace([np.NaN],0)
            s = s.rolling(100, center=True).mean()
            s = s.dropna()
            
            # Normalize to rpm
            s = s.divide(totals[bam])
            if log_scale:
                s = s.apply(np.log2)
  
            # Plot!
            ax[n].bar(s.index, s, linewidth=1, color=colors[n], edgecolor=colors[n], zorder=2)
            ax[n].tick_params(axis='both', which='major', labelsize=14)
            
            max_y = max([max_y,max(s)])
        
        # Make transcript boxes
        try:
            for tx in r['transcript'].split(','):
                title = r['transcript'].rstrip(',')
                if len(tx) > 0:
                    if tx+tx_suffix in tx_dict:
                        strand = GT.gene_patches(tx, tx_dict, ax[-1])
                        ax[-1].set_xlim(start, end)
                    else:
                        try:
                            new_tx = tx.split(' ')[0]
                            if new_tx[-2] == 'T' or new_tx[-2] == '.':
                                new_tx = new_tx[:-2]
                            strand = GT.gene_patches(new_tx, tx_dict, ax[-1])
                            ax[-1].set_xlim(start, end)
                        except KeyError:
                            print "Transcript unknown"
        except AttributeError:
            title = ix
                
        # Set x limits
        for n in range(len(bam_list)):
            ax[n].set_xlim(start, end)
            if same_yaxis:
                ax[n].set_ylim(0,max_y+0.1*max_y)

        # Label axes
        ax[0].set_ylabel('RPM', fontsize=16)
        ax[0].set_title(title, fontsize=16)
        plt.show()
        
        # Save if indicated
        if save_dir is not None:
            if not save_dir.endswith('/'): save_dir = save_dir+'/'
                
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
                
            fig.savefig(save_dir+'_'.join(title.split(','))+'.eps', format='eps')
            
        plt.clf()

###################################################################
## Tile functions for analyzing ChIP-seq data across chromosomes ##
###################################################################
                
def read_feature_gff3(gff3):
    feature_dict = {}
    with open(gff3) as f:
        for line in f:
            data = line.split('\t')
            try:
                chrom = data[0]
                start = int(data[3])
                end = int(data[4])
                name = data[8].split('ID=')[-1].strip()
            except IndexError:
                pass
            
            feature_dict[name] = [chrom, start, end]
    
    return feature_dict

def read_fa(fa):
    fa_dict = {}
    with open(fa) as f:
        for line in f:
            if line.startswith('>'):
                contig = line.lstrip('>').strip()
                fa_dict[contig] = ''
            else:
                fa_dict[contig] = fa_dict[contig]+line.strip()
    return fa_dict
                
def make_tile_df(directory, tile_size=5000, organism='crypto', gff3=None, fa=None):
    '''Creates a spreadsheet that contains RPKM values for tiles across the entire genome from bam files. Can also
    classify each tile based on a gff3 file if desired (e.g. centromeres and telomeres).
    
    Parameters
    ----------
    directory : str
            Directory containing bam files and index files
    tile_size : int, default 5000
            Size of tiles within each chromosome
    organism : str, default 'crypto'
            'crypto' or 'pombe' 
    gff3 : str, default `None`
            gff3 format file containing features. Tiles will be classified based on the boundaries described in this file
            
    Returns
    ------
    df : pandas.DataFrame
            Dataframe containing the RPKM calculations for each tile from each bam file in the directory'''
    
    # Find all bam files in directory
    bam_list = []
    if not directory.endswith('/'): directory = directory+'/'
    for file in os.listdir(directory):
        if file.endswith(('_sorted.bam')):
            bam_list.append(directory+file)
    
    # Load in tiles from dictionary
    if 'crypto' in organism.lower():
        if fa is None:
            with open('/home/jordan/GENOMES/H99_fa.json') as f: fa_dict = json.load(f)
        else:
            fa_dict = read_fa(fa)
        length_dict = {}
        for chrom, seq in fa_dict.iteritems():
            length_dict[chrom] = len(seq)
    elif 'pombe' in organism.lower():
        if fa is None:
            with open('/home/jordan/GENOMES/POMBE/Sp_fasta_dict.json') as f: fa_dict = json.load(f)
        else:
            fa_dict = read_fa(fa)
        length_dict = {}
        for chrom, seq in fa_dict.iteritems():
            length_dict[chrom] = len(seq)
    else:
        print "Organism not supported at this time"
        return None
    
    frag_list = []
    for chrom, length in length_dict.iteritems():
        n=0
        counter = 1
        while n < length-tile_size:
            frag_list.append([chrom+'-'+str(counter),n,n+tile_size-1,'+',chrom])
            n += tile_size
            counter += 1
        frag_list.append([chrom+'-'+str(counter),n,length,'+',chrom])
    
    all_df = pd.DataFrame(index=[x[0] for x in frag_list])
    all_df.loc[:,'chrom'] = [x[4] for x in frag_list]
    all_df.loc[:,'start'] = [x[1] for x in frag_list]
    all_df.loc[:,'end'] = [x[2] for x in frag_list]
    
    # Label categories
    if gff3 is not None:
        feature_dict = read_feature_gff3(gff3)
        
        features =[]
        for ix, r in all_df.iterrows():
            match = False
            for feat, info in feature_dict.iteritems():
                if r['chrom'] == info[0]:
                    if r['start'] in range(info[1],info[2]) or r['end'] in range(info[1],info[2]):
                        features.append(feat)
                        match = True
                        break
            if not match:
                features.append('')
        all_df.loc[:,'Feature'] = features
    
    
    # Open bam files for counting
    open_bams = {}
    for bam in bam_list:
        open_bams[bam] = pysam.Samfile(bam)

    fix_strand = {'+':'-','-':'+'}

    # Count reads in each tile
    count_dict = {}
    col_dict = {'name':[],'chrom':[],'start':[],'end':[]}
    for bam in bam_list:
        name = bam.split('/')[-1].split('_sorted')[0]
        print name
        count_dict[name] = []
        for frag, start, end, strand, chrom in frag_list:
            
            # Count reads on both strands
            count = GT.count_reads_in_window(open_bams[bam], chrom, start, end, strand)
            count += GT.count_reads_in_window(open_bams[bam], chrom, start, end, fix_strand[strand])
            
            # Divide by tile size in kb
            count = count/((end-start)/1000.)
            count_dict[name].append(count)

    # Populate dataframe
            
    for name, data in count_dict.iteritems():
        all_df[name] = data
        
        # Calculate RPKM
        bam_file = [x for x in bam_list if name in x]
        if len(bam_file) > 1:
            print "Redundant bam names, cannot calculate RPKM"
        else:
            bam_file = bam_file[0]
        total_reads = GT.count_aligned_reads(bam_file)
        all_df[name+' RPKM'] = [x/total_reads for x in data]

    # Remove weird chromosome
    #all_df = all_df[~all_df.index.str.contains('Caalf')]
    


    return all_df

def normalize_tiles_to_WCE(df, WCE_sample_dict):
    '''Adds normalized columns to the spreadsheet created by make_tile_df
    
    Parameters
    ----------
    df : pandas.DataFrame
            Dataframe generated by make_tile_df function (modifies this dataframe in place
    WCE_sample_dict : dict
            format - {WCE1:[sample1a,sample1b],WCE2:[sample2a,sample2b]}'''
    
    # format for WCE_sample_dict is {WCE_name:[list of bam files], WCE_name:[list of bam files]}
    
    RPKM_cols = [x for x in df.columns if x.endswith('RPKM')]
    
    for WCE, samples in WCE_sample_dict.iteritems():
        WCE_RPKM = [x for x in RPKM_cols if WCE in x][0]
        
        for sample in samples:
            sample_RPKM = [x for x in RPKM_cols if sample in x][0]
            df.loc[:,sample+' RPKM norm'] = df[sample_RPKM]/df[WCE_RPKM]
            
def tile_scatters(df, name1, name2, highlight=None, colors=['orangered','0.5']):
    '''Creates scatter plots from dataframe generated by make_tile_df and normalize_tiles_to_WCE
    
    Parameters
    ----------
    df : pandas.DataFrame
            Dataframe generated by make_tile_df function (modifies this dataframe in place
    name1 : str
            all or part of the name of the first sample 
    name2 : str
            all or part of the name of the second sample
    highlight : str, default `None`
            feature to highlight (e.g. 'Cen') - note string must be contained within the feature names provided in the gff3 file
    colors : list, default ['orangered','0.5']
            colors for scatter plot - first color is the highlight color
            
    Returns
    ------
    fig : matplotlib.figure
            matplotlib figure object that can be exported with fig.savefig(filename, format='pdf')'''
    
    cols = [x for x in df.columns if x.endswith('RPKM norm')]
    a = [x for x in cols if name1 in x][0]
    b = [x for x in cols if name2 in x][0]
    
    
    fig, ax = plt.subplots(figsize=(4,4))
    ax.scatter(df[a].apply(np.log2), df[b].apply(np.log2), s=15, alpha=0.5, label='Genome', color=colors[1], zorder=1)
    if highlight is not None:
        c = df[df['Feature'].str.contains(highlight)][a].apply(np.log2)
        d = df[df['Feature'].str.contains(highlight)][b].apply(np.log2)
        print len(c)
        ax.scatter(c, d, s=15, color=colors[0], label=highlight, zorder=2, alpha=0.5)

    min_xy = min(ax.get_xlim()[0], ax.get_ylim()[0])
    max_xy = max(ax.get_xlim()[1], ax.get_ylim()[1])

    ax.set_xlim(min_xy,max_xy)
    ax.set_ylim(min_xy,max_xy)

    ax.set_xlabel(a.split('_S')[0]+' (log2)', fontsize=14)
    ax.set_ylabel(b.split('_S')[0]+' (log2)', fontsize=14)
    ax.legend(fontsize=14)

    ax.plot((min_xy,max_xy),(min_xy,max_xy), '--', zorder=0, color='0.5')
    
    plt.show()
    return fig

def tile_boxplots(df, name1, name2, feature=None, color='0.5', yscale=None):
    '''Creates boxplots from dataframe generated by make_tile_df and normalize_tiles_to_WCE
    
    Parameters
    ----------
    df : pandas.DataFrame
            Dataframe generated by make_tile_df function (modifies this dataframe in place
    name1 : str
            all or part of the name of the first sample 
    name2 : str
            all or part of the name of the second sample
    feature : str, default `None`
            feature to highlight (e.g. 'Cen') - note string must be contained within the feature names provided in the gff3 file
    color : str, default '0.5'
            color of fill for boxplot
    yscale : tuple, default `None`
            Provide if autoscale is not sufficient on the y-axis (due to fliers). E.g. (0,5)
            
    Returns
    ------
    fig : matplotlib.figure
            matplotlib figure object that can be exported with fig.savefig(filename, format='pdf')'''
        
    cols = [x for x in df.columns if x.endswith('RPKM norm')]
    a = [x for x in cols if name1 in x][0]
    b = [x for x in cols if name2 in x][0]
    box_data = df[[a,b,'Feature']]

    if feature is not None:
        feat_box = box_data[box_data['Feature'].str.contains(feature)]
        oth_box = box_data[~box_data['Feature'].str.contains(feature)]

        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8,4), sharey=True)
        sns.boxplot(data=feat_box, ax=ax1, color=color)
        sns.boxplot(data=oth_box, ax=ax2, color='0.5')

        ax1.set_ylabel('Enrichment over WCE (RPKM)', fontsize=14)
        ax1.set_title(feature, fontsize=14)
        ax2.set_title('Genome', fontsize=14)
        
        if yscale is not None:
            ax1.set_ylim(yscale)
            
        print 'p-value for feature:'
        print stats.mannwhitneyu(feat_box[a],feat_box[b])[1]
        print 'p-value for remainder:'
        print stats.mannwhitneyu(oth_box[a],oth_box[b])[1]

    else:
        fig, ax = plt.subplots(figsize=(4,4))
        sns.boxplot(box_data, ax=ax, color=color)
        
        ax.set_ylabel('Enrichment over WCE (RPKM)', fontsize=14)
        
        if yscale is not None:
            ax.set_ylim(yscale)
    
    fig.tight_layout()
    plt.show()
    return fig

############################
## Metatelomere functions ##
############################

def tel_read_counts(open_bam, row, tel_len=60000):
    ''' Function for counting reads in telomeres. Called by generate_metatel_df by way of compile_telomere_series'''
    s = pd.Series(1, index=np.arange(tel_len))
    if row['name'].endswith('L'):
        reads = open_bam.fetch(row['chromosome'],0,tel_len)
    elif row['name'].endswith('R'):
        reads = open_bam.fetch(row['chromosome'],row['end']-tel_len,row['end'])

    for read in reads:
        if row['name'].endswith('L'):
            read_range = np.arange(read.reference_start, read.reference_end)
        elif row['name'].endswith('R'):
            read_start = row['end']-read.reference_end
            read_end = row['end']-read.reference_start
            #print read_start
            read_range = np.arange(read_start, read_end)
        else:
            print "Can't parse telomere name: "+r['name']+"... skipping"
            pass

        s[s.index.isin(read_range)] += 1
        #print s
    return s

def cen_read_counts(open_bam, row, cen_flank=30000):
    ''' Function for counting reads in centromeres. Called by generate_meta_df via compile_cen_series'''
    s = pd.Series(1, index=np.arange(cen_flank*-1, cen_flank))
    center = row['start']+(row['end']-row['start'])/2
    reads = open_bam.fetch(row['chromosome'], center-cen_flank, center+cen_flank)

    for read in reads:
        read_range = np.arange(read.reference_start-center, read.reference_end-center)
        s[s.index.isin(read_range)] += 1
    return s

def compile_cen_tel_series((bam, cen_tel_gff3, cen_tel_len, cen_or_tel)):
    '''Function for compiling centromere or telomere series. Called by generate_meta_df'''
    print bam
    cen_tel = pd.read_csv(cen_tel_gff3, sep='\t', header=None, 
                          names=['chromosome','name','type','start','end','x','strand','y','ID']).dropna()
    open_bam = pysam.Samfile(bam)
    s_dict = {}
    for ix, r in cen_tel.iterrows():
        if cen_or_tel == "telomere" and r['name'].lower().startswith('tel'):
            s_dict[r['name']] = tel_read_counts(open_bam, r, tel_len=cen_tel_len)
        elif cen_or_tel == 'centromere' and r['name'].lower().startswith('cen'):
            s_dict[r['name']] = cen_read_counts(open_bam, r, cen_flank=cen_tel_len)
    
    s_dict['total aligned reads'] = GT.count_aligned_reads(bam)
    
    return (bam, s_dict)

def cen_tel_bins(sample, chip_dict, wce_dict, cen_tel_len=60000, bin_size=5000, cen_or_tel='telomere'):
    bin_df = pd.DataFrame(columns=['name','sample','bin','RPKM'])

    bins = []
    if cen_or_tel == 'telomere':
        n=0
    elif cen_or_tel == 'centromere':
        n=cen_tel_len*-1
        
    while n <= cen_tel_len-bin_size:
        bins.append(np.arange(n,n+bin_size-1))
        n += bin_size
    if cen_tel_len-n >= 100:
        bins.append(np.arange(n,cen_tel_len))
    
    for reg, s in chip_dict.items():
        if reg != 'total aligned reads':
            for n, b in enumerate(bins):
                bin_s_chip = s[s.index.isin(b)]
                wce_s = wce_dict[reg]
                bin_s_wce = wce_s[wce_s.index.isin(b)]
                RPKM = sum(bin_s_chip)/(chip_dict['total aligned reads']) / sum(bin_s_wce)/(wce_dict['total aligned reads'])
                new_row = pd.Series([reg, sample, n, RPKM], index=['name','sample','bin','RPKM'], name=n)
                bin_df = bin_df.append(new_row)
    return bin_df

def create_meta_cen_tel(chip_dict, wce_dict, chip_name, cen_tel_len=60000, cen_or_tel='telomere'):
    '''Function for creating metacentromere or metatelomere from region series (generated by compile_cen_tel_series).
    Called by generate_meta_df.'''
    if cen_or_tel == 'telomere':
        avg = pd.Series(0, index=np.arange(cen_tel_len))
    elif cen_or_tel == 'centromere':
        avg = pd.Series(0, index=np.arange(cen_tel_len*-1,cen_tel_len))
    
    for reg, s in chip_dict.items():
        if reg != 'total aligned reads':
            s = s.divide(chip_dict['total aligned reads']) / wce_dict[reg].divide(wce_dict['total aligned reads'])
            avg = avg+s
    avg = avg.divide(len(chip_dict)-1)
    avg.name = chip_name
    return avg

def generate_meta_df(chip_wce_pairs, out_name, cen_or_tel, threads=4, tel_len=60000, cen_flank=30000, bin_size=5000, cen_tel_gff3='/home/jordan/GENOMES/Cen_Tel.gff3.txt'):
    '''Creates spreadsheet with metatelomeres for each ChIP and WCE pair provided. 
    Should work for any organism as long as gff3 file contig names match bam file contig names.
    Note: Second column in gff3 file must indicate which rows are telomeres. 
          Telomere names must start with 'tel' or 'Tel' and end with 'R' or 'L' indicating right or left.
    
    Parameters
    ----------
    chip_wce_pairs : list of tuples
            Provide file location of ChIP and WCE sorted, indexed bam files in pairs
            Example: [('chip1_sorted.bam','wce1_sorted.bam'),('chip2_sorted.bam','wce2_sorted.bam'),('chip3_sorted.bam','wce3_sorted.bam')]
    out_name : str
            Indicate name to give to csv file
    cen_or_tel : str
            'centromere' or 'telomere'
    threads : int, default 4
            Number of threads to use for parallel processing
    tel_len : int, default 60000
            Maximum size of telomere, 60000 is based on previous descriptions of telomeres in C. neoformans
    cen_flank : int, default 30000
            Maximum size of one half of centromere, 30000 is based on previous descriptions of telomeres in C. neoformans
    bin_size : int, default 5000
            Bin size for determining difference between sections of telomeres or centromeres. E.g. a bin size of 5000
            with a tel_len of 60000 gives 12 bins
    cen_tel_gff3 : str, default '/home/jordan/GENOMES/Cen_Tel.gff3.txt'
            GFF3 file indicating telomere locations. See note above for formatting
            
    Returns
    ------
    df : pandas.DataFrame
            DataFrame object with metatelomeres for all ChIP bam files
    bins : pandas.DataFrame
            DataFrame object containing RPKM values for each telomere for each sample in indicated bins
    Output
    ------
    csv_file : Comma separated table containing metatelomeres for all ChIP bam files
    '''
    
    if cen_or_tel == 'centromere':
        reg_len = cen_flank
        df = pd.DataFrame(index=np.arange(cen_flank*-1, cen_flank))
    elif cen_or_tel == 'telomere':
        reg_len = tel_len
        df = pd.DataFrame(index=np.arange(tel_len))
    else:
        print "Must indicate 'centromere' or 'telomere'"
        return None
    

    wce_tel_dicts = {}
    
    bam_set = set()
    for chip, wce in chip_wce_pairs:
        bam_set.add((wce, cen_tel_gff3, reg_len, cen_or_tel))
        bam_set.add((chip, cen_tel_gff3, reg_len, cen_or_tel))
    bam_list = list(bam_set)
    if len(bam_list) < threads: threads = len(bam_list)
    p = Pool(threads)
    series = p.map(compile_cen_tel_series, bam_list)
    
    chip_ser = {x[0]:x[1] for x in series if x[0] in zip(*chip_wce_pairs)[0]}
    wce_ser = {x[0]:x[1] for x in series if x[0] in zip(*chip_wce_pairs)[1]}
    
    with open(out_name+'.pickle', 'w') as fout:
        pickle.dump((chip_ser, wce_ser), fout)
    
    bins = None
    for chip, chip_dict in chip_ser.items():
        chip_name = '_'.join(chip.split('/')[-1].split('_')[:-1])
        print chip_name
        n = [x[0] for x in chip_wce_pairs].index(chip)
        wce_name = [x[1] for x in chip_wce_pairs][n]
        wce_dict = wce_ser[wce_name]
        s = create_meta_cen_tel(chip_dict, wce_dict, chip_name, cen_tel_len=reg_len, cen_or_tel=cen_or_tel)
        df.loc[:,s.name] = s
        
        if bins is None:
            bins = cen_tel_bins(chip_name, chip_dict, wce_dict, cen_tel_len=reg_len, bin_size=bin_size, cen_or_tel=cen_or_tel)
        else:
            bins = bins.append(cen_tel_bins(chip_name, chip_dict, wce_dict, cen_tel_len=reg_len, bin_size=bin_size, cen_or_tel=cen_or_tel))
    
    bins = bins.reset_index(drop=True)
    bins.to_csv(out_name+'_bins.csv')
    
    df.to_csv(out_name+'.csv')
    return df, bins

def plot_metaregion(df, samples, out_name, colors=['#222f3e','#ee5253','#5f27cd','#10ac84'], smooth_window=1000):
    '''Plots metatelomere generated from generate_metatel_df. 
    
    Parameters
    ----------
    df : pandas.DataFrame
            DataFrame object generated by generate_megatel_df function.
            Note, to create this object from the csv file, first import pandas, then use df = pandas.read_csv('my_metateomere.csv')
    samples : list of strings
            list of sample names (must match df columns, so take a look at your DataFrame first)    
    out_name : str
            Indicate name to give to pdf of plot
    colors : list of strings, default ['#222f3e','#ee5253','#5f27cd','#10ac84']
            List of colors. Must be at least as long as samples
            Use either matplotlib color names (https://i.stack.imgur.com/fMx2j.png) or flatui colors (https://flatuicolors.com/)
    smooth_window : int, default 1000
            Size of window in bp over which to smooth metatelomere. 

    Output
    ------
    pdf : PDF formatted figure of your metatelomere plot
    '''
    
    sns.set(font_scale=1.2)
    sns.set_style('whitegrid')
    fig, ax = plt.subplots(1)
    for n, sample in enumerate(samples):
        ax.plot(df.index, df[sample].rolling(smooth_window).mean(), color=colors[n])
    ax.set_ylabel('Fold enrichment over WCE')
    ax.set_xlabel('Distance from end (bp)')
    
    plt.show()
    fig.savefig(out_name+'.pdf', format='pdf', bbox_inches='tight')
    plt.clf()
    
def plot_bin_diff(bins, samples, out_name, colors=['#222f3e','#ee5253','#5f27cd','#10ac84']):
    '''Plots telomere or centromere bins per sample and reports p-values for pairwise sample comparisons (T-test) 
    
    Parameters
    ----------
    bins : pandas.DataFrame
            DataFrame object generated by generate_megatel_df function.
            Note, to create this object from the csv file, first import pandas, then use df = pandas.read_csv('my_metateomere_bins.csv')
    samples : list of strings
            list of sample names (must match df columns, so take a look at your DataFrame first)    
    out_name : str
            Indicate name to give to pdf of plot

    Output
    ------
    pdf : PDF formatted barplot of each pairwise bin comparison between samples
    csv : CSV file containing p-values for each bin comparison pairwise between samples
    '''
    bins = bins[bins['sample'].isin(samples)]
    sns.set_style('darkgrid')
    sns.set(font_scale=1.2)
    
    df = pd.DataFrame(columns=['sample 1','sample 2','bin','p-value'])
    
    sample_pairs = []
    for n in range(len(samples)-1):
        for m in range(1,len(samples)-n):
            sample_pairs.append((samples[n], samples[n+m]))
        
    for s1, s2, in sample_pairs:
        sub_bins = bins[bins['sample'].isin([s1,s2])]
        sns.barplot(x='bin', y='RPKM', data=sub_bins, hue='sample', 
                palette=sns.color_palette(colors), saturation=0.9)
    
        plt.savefig(out_name+'_'+s1+'_'+s2+'_bin_comparison.pdf', format='pdf', bbox_inches='tight')
        plt.show()
        plt.clf()
        
        for b in set(bins['bin']):
            sub_bins1 = bins[(bins['bin'] == b) & (bins['sample'] == s1)]
            sub_bins2 = bins[(bins['bin'] == b) & (bins['sample'] == s2)]
            p = stats.ttest_ind(sub_bins1['RPKM'],sub_bins2['RPKM'])[1]
            row = pd.Series([s1, s2, str(b), p], index=['sample 1','sample 2','bin','p-value'])
            df = df.append(row, ignore_index=True)
    
    df.to_csv(out_name+'_bin_comparison.csv')
    #return df