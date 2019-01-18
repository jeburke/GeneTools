import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.append('/home/jordan/RNA-is-awesome/')
import SPTools as SP
import pysam
from subprocess import check_output
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import pysam
from subprocess import check_output
import os
import numpy as np

def read_genes_from_csv(csv, header=True):
    gene_list = []
    with open(csv) as f:
        lines = f.readlines()
        if header is True:
            first = lines[0]
        else: first = None
        
        print len(lines)
        for line in lines:
            if line != first and len(line) > 2:
                data = line.split(',')
                if len(data) > 1:
                    gene = data[0]
                else:
                    data2 = data[0].split('\t')
                    if len(data2) > 1:
                        gene = data2[0]
                    else:
                        gene = data2[0].strip()
                try:
                    if gene[-2] == 'T' or gene[-2] == '.':
                        gene = gene[:-2]
                    gene_list.append(gene)
                except IndexError:
                    print "Gene name less than 2 characters. Check formatting!"
    print len(gene_list)
    return gene_list
        
def list_from_arg(arg, length):
    if type(arg) != list:
        value = arg
        arg = []
        for n in range(length):
            arg.append(value)
    return arg
    
def get_junctions(open_bam, chrom, start, end, strand, baseline=0):
    iterator = open_bam.fetch(chrom, start, end)
    intron_dict = {}
    
    for read in iterator:
        intron = False
        if len(read.cigartuples) == 3 and read.cigartuples[1][0] == 3:
            intron = True
        if intron is True:
            
            intron_start = read.reference_start + read.cigartuples[0][1]
            intron_end = intron_start + read.cigartuples[1][1] - 1
            middle = (intron_end-intron_start)/2 + intron_start
            
            if (intron_start, middle, intron_end) not in intron_dict:
                intron_dict[(intron_start, middle, intron_end)] = [0.1,0,0.1]
            intron_dict[(intron_start, middle, intron_end)][1] += 1
    
    return intron_dict

    
def gene_patches(tx, tx_dict, ax, arrow=False):
    iso_list = [x for x in tx_dict if tx in x]
    if len(iso_list) == 0:
        return None
    
    for n, iso in enumerate(iso_list):
        start, end, strand, CDS_start, CDS_end, exons, chrom = SP.tx_info(iso, tx_dict)
        if arrow is False:
            tx_patch = patches.Rectangle((start,0.8-n*0.15),end-start,0.04,edgecolor='0.1',facecolor='0.1')
            ax.add_patch(tx_patch)
        else:
            if strand == '+':
                ax.arrow(start, 0.9, end-start-0.02*(end-start), 0, linewidth=2, head_width=0.1, 
                         head_length=0.02*(end-start), fc='k', ec='k')
            elif strand == '-':
                ax.arrow(end, 0.9, start-end-0.02*(start-end), 0, linewidth=2, head_width=0.1, 
                         head_length=0.02*(end-start), fc='k', ec='k')

        if exons is not None:
            exon_patches = []
            for exon_start, exon_stop in exons:
                exon_patches.append(patches.Rectangle((exon_start, 0.775-n*0.15), exon_stop-exon_start, 0.10,
                                                      edgecolor='0.1',facecolor='0.1'))
            for patch in exon_patches:
                ax.add_patch(patch)
        else:
            CDS_patch = patches.Rectangle((CDS_start, 0.75-n*0.15),CDS_end-CDS_start, 0.10, edgecolor='0.1', facecolor='0.1')
            ax.add_patch(CDS_patch)
        ax.get_yaxis().set_ticks([])
    return strand  
    
def igv_plots_general(bam_list, gene_list, organism, colors=None, names=None, save_dir=None, 
                      unstranded=False, end_only=False, same_yaxis=False, specific_range=None, transcript_direction=True,
                     log_scale=False, rpm=True, PE=False, plot_junctions=False):
    '''Usage:
    Parameters
    ----------
    bam_list : list, bam files in order of plotting (top to bottom)
    gene_list : list of transcripts to plot (should be genes not transcript isoforms)
            if dataframe passed instead of list, will plot introns (must have intron information in datafame)
    organism : str, pombe or crypto
    colors : list, default `None`
            list of colors to use, same length as bam_list, check matplotlib documentation for valid color names
    names : list, default `None`
            list of sample names to use instead of bam file names. Same length as bam_files
    save_dir : str, default `None`
            directory to save eps files. If None, does not save files
    unstranded : bool, default `False`
            Use True for ChIP or DNA sequencing data (or unstranded RNAseq)
    end_only : bool or list, default `False`
            Whether to plot only the ends of reads. If different for each bam, make a list of bools same length as bam_list
    same_yaxis : bool, default `False`
            Whether all samples should be plotted on the same axis after normalizing to total number of aligned reads
    specific_range : str, default `None`
            Options: ('end', window)
                     ('start', window)
                     ([coordinate], window)
    transcript_direction : bool, default `True`
            If True, will plot in the direction of transcription, not in the direction of the DNA
    '''
    
    # Get all organism information (annotation etc.)
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    fix_info = {'I':'chr1','II':'chr2','III':'chr3','chr1':'I','chr2':'II','chr4':'IV','chr5':'V','chr6':'VI',
                'chr7':'VII','chr8':'VIII','chr9':'IX','chr10':'X','chr11':'XI','chr12':'XII','chr13':'XIII',
                'chr14':'XIV','chr15':'XV','chr16':'XVI','-':'+','+':'-','chr1':'I','chr2':'II','chr3':'III'}
    if organism == 'pombe':
        tx_suffix = '.1'
    else:
        tx_suffix = 'T0'
    
    # Set up range parameters if specific range is indicated
    if specific_range is not None:
        window = int(specific_range[1])
        new_tx_dict = {}
        for gene in gene_list:
            info = tx_dict[gene+tx_suffix]
            if specific_range[0] == 'end':
                if info[2] == '+':
                    start = info[1]-window
                    end = info[1]+window
                else:
                    start = info[0]-window
                    end = info[0]+window
            elif specific_range[0] == 'start':
                if info[2] == '-':
                    start = info[1]-window
                    end = info[1]+window
                else:
                    start = info[0]-window
                    end = info[0]+window    
            else:
                start = int(specific_range[0])-window
                end = int(specific_range[0])+window           
            new_tx_dict[gene+tx_suffix] = [start, end, info[2], info[3]]
    else:
        new_tx_dict = tx_dict
                
    # Open bam files and count reads if rpm is True
    open_bams = {}
    total_list = []
    for bam in bam_list:
        open_bams[bam] = pysam.Samfile(bam)
        if rpm is True:
            total = check_output(['samtools','view','-F 0x04','-c',bam]).strip()
            total = float(total)/1000000.
            total_list.append(total)
        else:
            total_list.append(1.)
    
    # Expand optional arguments to lists if necessary
    colors = list_from_arg(colors, len(bam_list))
    end_only = list_from_arg(end_only, len(bam_list))
    log_scale = list_from_arg(log_scale, len(bam_list))
    unstranded = list_from_arg(unstranded, len(bam_list))
    
    # Get gene_list from dataframe if gene_list is not a list
    df = None
    if type(gene_list) == dict:
        new_tx_dict = gene_list
        gene_list = gene_list.keys()
        
    elif type(gene_list) != list:
        df = gene_list
        gene_list = df.index
    
    for tx in gene_list:
        num_ax = len(bam_list)+1
        if plot_junctions is True:
            num_ax += len(bam_list)
        
        fig, ax = plt.subplots(num_ax, figsize=(10,num_ax), sharex=True)
        fig.subplots_adjust(hspace=0)
        
        # Get transcript info from transcript_dictionary
        if df is None:
            try:
                info = new_tx_dict[tx+tx_suffix]
            except KeyError:
                info = new_tx_dict[tx]
            chrom = info[3]
            start = info[0]
            end = info[1]
            strand = info[2]
        
        # If dataframe was passed, get plotting information from dataframe instead
        else:
            if isinstance(df.columns, pd.core.index.MultiIndex):
                new_columns = [x[1] for x in df.columns if x[0] == 'Peaks']
                df = df[[x for x in df.columns if x[0] == 'Peaks']]
                df.columns = new_columns
            strand = df.loc[tx,'strand']
            chrom = df.loc[tx,'chromosome']
            if strand == '+':
                start = df.loc[tx,'position']-100
                end = df.loc[tx,'position'] + df.loc[tx,'intron size']+100
            elif strand == '-':
                start = df.loc[tx,'position']-df.loc[tx,'intron size']-100
                end = df.loc[tx,'position']+100
            start = int(start)
            end = int(end)
            
            tx = df.loc[tx,'transcript']
        
        # Generate read series for each transcript
        max_y = 0
        junc_ymax = 0
        for n, bam in enumerate(bam_list):
            try:
                bam_iter = open_bams[bam].fetch(chrom, start, end)
            except ValueError:
                chrom = fix_info[chrom]
                bam_iter = open_bams[bam].fetch(chrom, start, end)
            if end_only[n] is not False:
                s = SP.generate_read_series_A(bam_iter, chrom, start, end, strand)
                linewidth = 2
            else:
                if PE is False:
                    s = SP.generate_read_series_B(bam_iter, chrom, start, end, strand)
                else:
                    s = SP.generate_read_series_PE(bam_iter, chrom, start, end, strand)
                linewidth = 1
            
            # Get reads from otherstrand if the library type is unstranded
            if unstranded[n] is True:
                bam_iter = open_bams[bam].fetch(chrom, start, end)
                if end_only[n] is not False:
                    s2 = SP.generate_read_series_A(bam_iter, chrom, start, end, fix_info[strand])
                    linewidth = 2
                else:
                    if PE is False:
                        s2 = SP.generate_read_series_B(bam_iter, chrom, start, end, fix_info[strand])
                    else:
                        s2 = SP.generate_read_series_PE(bam_iter, chrom, start, end, fix_info[strand])
                    linewidth = 1
                s = s.add(s2)
            
            # Normalize to rpm (will just divide by 1 if rpm is False)
            s = s.divide(total_list[n])
            if log_scale[n] is True:
                s = s.apply(np.log2)
            
            # Plot!
            ax[n].bar(s.index, s, linewidth=linewidth, color=colors[n], edgecolor=colors[n], zorder=2)
            ax[n].tick_params(axis='both', which='major', labelsize=14)
            
            max_y = max([max_y,max(s)])
            
            if plot_junctions is True:
                m = n+len(bam_list)
                intron_dict = get_junctions(open_bams[bam], chrom, start, end, strand)
                ax[m].plot((start, end),(0,0),'-',c='k')
                for coords, heights in intron_dict.iteritems():
                    ax[m].plot(coords, heights, '-', linewidth=2, color=colors[n])
                    ax[m].fill_between(coords, 0, heights, facecolor=colors[n], interpolate=True, alpha=0.5)
                if same_yaxis is True:
                    junc_ymax = max([junc_ymax, max(zip(*intron_dict.values())[1])])
            
        # Add diagram of gene below traces
        if tx in tx_dict:
            strand = gene_patches(tx, tx_dict, ax[-1])
            ax[-1].set_xlim(start, end)
        else:
            try:
                new_tx = tx.split(' ')[0]
                if new_tx[-2] == 'T' or new_tx[-2] == '.':
                    new_tx = new_tx[:-2]
                strand = gene_patches(new_tx, tx_dict, ax[-1])
                ax[-1].set_xlim(start, end)
            except KeyError:
                print "Transcript unknown"
                
        
        # Flip minus strand transcripts if indicated
        if transcript_direction is True:
            if strand == '-':
                ax[-1].invert_xaxis()

        # Set x and y limits
        for n in range(len(bam_list)):
            ax[n].set_xlim(start, end)
            if same_yaxis is True:
                ax[n].set_ylim(0,max_y+0.1*max_y)
                
                if plot_junctions is True:
                    ax[n+len(bam_list)].set_ylim(0,junc_ymax+0.1*junc_ymax)
            
            if strand == '-':
                ax[n].invert_xaxis()

        ax[0].set_ylabel('RPM', fontsize=16)
        ax[0].set_title(tx, fontsize=16)
        #ax[0].get_xaxis().set_ticks([])
        plt.show()
        
        # Save if indicated
        if save_dir is not None:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            fig.savefig(save_dir+tx+'.eps', format='eps')
            
        plt.clf()