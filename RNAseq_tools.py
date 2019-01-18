import sys
import os
import subprocess
import pandas as pd
from sklearn import cluster
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib_venn import venn2, venn2_circles
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.append('/home/jordan/RNA-is-awesome/')
import SPTools as SP
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/GeneTools/')
sys.path.append('/home/jordan/RNA-is-awesome/GeneTools/')
import Annotation_tools
import SeqTools
from math import factorial
from decimal import Decimal
import json
import matplotlib.patches as patches
from collections import OrderedDict
import seaborn as sns
from scipy import stats
import hdbscan

font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

def align_fastq(directory, threads=1, organism='crypto', PE=False):
    '''Automatically aligns all fastq.gz files in a drectory using TopHat'''
    if directory[-1] != '/':
        directory = directory+'/'
        
    if 'crypto' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/H99_bt2'
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    elif 'cerev' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/S288C/S288C_2'
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        if 'DBK' in organism:
            bowtie_ix = '/home/jordan/GENOMES/S288C/S288C_DBK'
    elif 'pombe' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/POMBE/Spombe-2'
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
    elif organism == 'candida':
        bowtie_ix = '/home/jordan/GENOMES/C_albicans2'
        gff3 = '/home/jordan/GENOMES/C_albicans_SC5314_version_A21-s02-m09-r10_features.gff'
        
    fastq_list = []
    for file in os.listdir(directory):
        if file.lower().endswith("fastq.gz") or file.lower().endswith("fastq"):
            fastq_list.append(directory+file)
            
    if PE is True:
        fastq_list = sorted(fastq_list)
        fastq_list2 = []
        for fastq in fastq_list:
            if '_R1_' in fastq:
                R2_name = fastq.split('_R1_')[0]+'_R2_'+fastq.split('_R1_')[1]
                if R2_name.split('/')[-1] not in os.listdir(directory):
                    print "No paired fastq file for "+fastq
                    print "Aborting!"
                    return None
                else:
                    fastq_list2.append(R2_name)
                    fastq_list.remove(R2_name)
        print len(fastq_list)
        print len(fastq_list2)
    
    ##Check to make sure tophat hasn't already completed on these files.
    for name in os.listdir(directory):
        if os.path.isdir(directory+name):
            for fastq in fastq_list:
                if directory+name == fastq.split('.fastq')[0]:
                    if 'accepted_hits.bam' in os.listdir(directory+name):
                        fastq_list.remove(fastq)
                        print "Tophat already completed on "+fastq
    
    if len(fastq_list) > 0:
        for n, fastq in enumerate(fastq_list):
            prefix = fastq.split('.fastq')[0]
            args = 'tophat2 --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 --min-anchor-length 8 --splice-mismatches 0 --min-intron-length 20 --max-intron-length 2000 --max-insertion-length 3 --max-deletion-length 3 --num-threads {0} --max-multihits 2 --library-type fr-firststrand --segment-mismatches 3 --no-coverage-search --segment-length 20 --min-coverage-intron 10 --max-coverage-intron 100000 --min-segment-intron 50 --max-segment-intron 500000 -G {1} -o {2} {3} {4}'.format(str(threads), gff3, prefix, bowtie_ix, fastq)
            
            if PE is True:
                args = args+' '+fastq_list2[n]
            
            print args
            p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
            p.wait()
        
            print '\n'
            print prefix

            if 'align_summary.txt' in os.listdir(prefix):
                with open(prefix+'/align_summary.txt','r') as f:
                    for line in f:
                        print line
                        print '\n'
            
    print 'Finished'
    
def align_fastq_STAR(directory, threads=1, organism='crypto', PE=False, no_introns=False):
    '''Automatically aligns all fastq.gz files in a drectory using STAR'''
    if directory[-1] != '/':
        directory = directory+'/'
        
    if 'crypto' in organism.lower():
        star_ix = '/home/jordan/GENOMES/STAR/CNA3/'
        gff3 = '/data/jade/H99/150901_Chr.GTF'
    elif 'cerev' in organism.lower():
        #gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        print "STAR not supported yet for S. cerevisiae"
        return None
    elif 'pombe' in organism.lower():
        #bowtie_ix = '/home/jordan/GENOMES/POMBE/Spombe-2'
        star_ix = '/home/jordan/GENOMES/STAR/POMBE'
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
    elif organism == 'candida':
        gff3 = '/home/jordan/GENOMES/C_albicans_SC5314_version_A21-s02-m09-r10_features.gff'
        print "STAR not supported yet for C. albicans"
        return None
        
    fastq_list = []
    for file in os.listdir(directory):
        if file.lower().endswith("fastq.gz") or file.lower().endswith("fastq"):
            bam1 = file.split('.fastq')[0]+'_Aligned.sortedByCoord.out.bam'
            bam2 = file.split('.fastq')[0]+'_sorted.bam'
            if bam1 not in os.listdir(directory) and bam2 not in os.listdir(directory):
                if PE and '_R2_' in file:
                    R1 = file.split('_R2_')[0]+'_R1_'+file.split('_R2_')[1]
                    bam1 = R1.split('.fastq')[0]+'_Aligned.sortedByCoord.out.bam'
                    bam2 = R1.split('.fastq')[0]+'_sorted.bam'
                    if bam1 not in os.listdir(directory) and bam2 not in os.listdir(directory):
                        fastq_list.append(directory+file)
                else:
                    fastq_list.append(directory+file)
            else:
                print "STAR already completed on "+file
            
    if PE is True:
        fastq_list = sorted(fastq_list)
        fastq_list2 = []
        for fastq in fastq_list:
            if '_R1_' in fastq:
                R2_name = fastq.split('_R1_')[0]+'_R2_'+fastq.split('_R1_')[1]
                if R2_name.split('/')[-1] not in os.listdir(directory):
                    print "No paired fastq file for "+fastq
                    print "Aborting!"
                    return None
                else:
                    fastq_list2.append(R2_name)
                    fastq_list.remove(R2_name)
    
    if len(fastq_list) > 0 and fastq_list[0].endswith('gz'):
        processes = []
        for n, fastq in enumerate(fastq_list):
            gzip_args = 'gunzip -k {0}'.format(fastq)
            print gzip_args
            processes.append(subprocess.Popen(gzip_args, shell=True, universal_newlines=True))
                
            if PE:
                fastq2 = fastq_list2[n]
                gzip_args = 'gunzip -k {0}'.format(fastq2)
                print gzip_args
                processes.append(subprocess.Popen(gzip_args, shell=True, universal_newlines=True))
            
            for p in processes:
                p.wait()
        

            prefix = fastq.split('.fastq')[0]+'_'

            if no_introns:
                args = 'STAR --alignIntronMax 1 --runThreadN {0} --genomeDir {1} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {2} --sjdbGTFfile {3} --readFilesIn {4}'.format(str(threads), star_ix, prefix, gff3, fastq.split('.gz')[0])
            else:
                args = 'STAR --runThreadN {0} --genomeDir {1} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {2} --sjdbGTFfile {3} --readFilesIn {4}'.format(str(threads), star_ix, prefix, gff3, fastq.split('.gz')[0])
            
            if PE:
                args = args+' '+fastq_list2[n].split('.gz')[0]
                
            print args
            p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
            p.wait()
        
            print '\n'
            print prefix

            prefix2 = prefix.split('/')[-1]
            if prefix2+'Log.final.out' in os.listdir(directory):
                with open(prefix+'Log.final.out','r') as f:
                    for line in f:
                        print line
            os.remove(fastq.split('.gz')[0])
            if PE:
                os.remove(fastq_list2[n].split('.gz')[0])
            
    print 'Finished'
    
def index_STAR_bam(directory):
    sorted_bams = []
    for name in os.listdir(directory):
        if name.endswith('Aligned.sortedByCoord.out.bam'):
            new_name = name.split('Aligned.sortedByCoord.out.bam')[0]+'sorted.bam'
            os.rename(directory+name, directory+new_name)
            sorted_bams.append(new_name)
    
    for bam in sorted_bams:
        args = 'samtools index {0}'.format(bam)
        print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
        p.wait()
        
def sort_index_bam_files(base_dir):
    '''Finds all TopHat output files and sorts and indexes BAM files'''
    try:
        check_samtools = subprocess.check_output(['samtools','--version'])
        version = check_samtools.split('\n')[0].split(' ')[1]
    except subprocess.CalledProcessError:
        version = 'old'
    print version
    
    
    if base_dir[-1] != '/':
        base_dir = base_dir+'/'
        
    tophat_out = []
    names = []
    for name in os.listdir(base_dir):
        if os.path.isdir(base_dir+name):
            if 'accepted_hits.bam' in os.listdir(base_dir+name):
                if name.split('/')[-1]+'_sorted.bam' not in os.listdir(base_dir):
                    tophat_out.append(base_dir+name+'/accepted_hits.bam')
                    names.append(base_dir+name)
        elif name.endswith('.bam') and 'sorted' not in name:
            if name.split('/')[-1].split('.bam')[0]+'_sorted.bam' not in os.listdir(base_dir):
                tophat_out.append(base_dir+name)
                names.append(base_dir+name)
    n=0
    for n in range(len(tophat_out)):
        name = names[n]
        args = 'samtools sort {0} -o {1}'.format(tophat_out[n], name+'_sorted.bam')
        if version == 'old':
            args = 'samtools sort {0} {1}'.format(tophat_out[n], name+'_sorted')
        print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
        p.wait()
    
    sorted_bams = []
    for file in os.listdir(base_dir):
        if file.endswith('_sorted.bam'):
            if file+'.bai' not in os.listdir(base_dir):
                sorted_bams.append(base_dir+file)
              
    for bam in sorted_bams:
        args = 'samtools index {0}'.format(bam)
        print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
        p.wait()
        
def count_with_HTseq(base_dir, organism='crypto'):
    '''Counts reads in every transcript with HTseq'''
    if base_dir[-1] != '/':
        base_dir = base_dir+'/'
    
    if 'crypto' in organism.lower():
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    elif 'cerev' in organism.lower():
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        if 'anti' in organism.lower():
            gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113_anti.gff3'
    elif 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
    
    
    bam_files = []
    for file in os.listdir(base_dir):
        if file.endswith('_sorted.bam') or file.endswith('sortedByCoord.out.bam'):
            if file.split('sorted')[0].rstrip('_')+'.htseq' not in os.listdir(base_dir):
                bam_files.append(file)
    
    for bam in bam_files:
        args = 'htseq-count -f bam -s reverse -t gene -i ID {0} {1}'.format(base_dir+bam, gff3)
        print args
        
        out = subprocess.check_output(args.split(' '))
        with open(base_dir+bam.split('sorted')[0].rstrip('_')+'.htseq','w') as fout:
            fout.write(out)

def count_for_QuantSeq(directory, organism=None, gff3=None, fa=None, threads=1, extend=200, insert_size=200, script_location='/home/jordan/CodeBase/RNA-is-awesome/QuantSeq_counting.py'):
    if organism is not None:
        if 'crypto' in organism.lower():
            gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
            fa = '/home/jordan/GENOMES/CNA3-gobs.fa'
        elif 'cerev' in organism.lower():
            gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
            fa = '/home/jordan/GENOMES/S288C/S288C_genome.fa'
        elif 'pombe' in organism.lower():
            gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
            fa = '/home/jordan/GENOMES/POMBE/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa'
    else:
        if fa is None or gff3 is None:
            print "Must provide fasta file and gff3 if not using built in organism!"
            return None
    
    processes = []
    for bam in [x for x in os.listdir('./') if x.endswith('bam')]:
        out_name = bam.split('_S')[0]
        if out_name+'.quant' not in os.listdir('./'):
            print out_name
            args = "python {0} --bam_file {1} --gff3 {2} --extend {3} --remove_pA_adj_reads --fasta {4} --insert_size {5} --name {4}".format(script_location, bam, gff3, extend, fa, insert_size, out_name)
            processes.append(subprocess.Popen(args.split(' ')))
        
        if len(processes) == threads:
            for p in processes:
                p.wait()
    
    for p in processes:
        p.wait()
            
def main():
    '''Runs the first part of this module to align, sort, index and count reads'''
    # Input will be RNAseq_tools.py directory number_of_threads_tophat --organism organism_name <--PE>
    # Current organism options: 'crypto', 'pombe', 'cerevisiae', 'candida'
    PE = False
    organism = None
    directory = './'
    threads = 1
    STAR = False
    no_introns = False
    quant_seq = False
    if '--help' in sys.argv or '-h' in sys.argv:
        print "Usage: python RNAseq_tools.py --directory directory --threads num_tophat_threads --organism organism_name <--PE> <--STAR> <--no_introns> <--QuantSeq --insert_size avg_library_insertion_size --extend extend_3UTR_n_bp>"
        return None
    
    for n, arg in enumerate(sys.argv):
        if arg == '--organism':
            if 'crypto' in sys.argv[n+1]:
                organism = 'crypto'
            elif 'pombe' in sys.argv[n+1]:
                organism = 'pombe'
            elif 'cerev' in sys.argv[n+1]:
                organism = 'cerevisiae'
            elif 'candida' in sys.argv[n+1]:
                organism = 'candida'
        elif arg == '--directory':
            directory = sys.argv[n+1]
        elif arg == '--threads':
            threads = int(sys.argv[n+1])
    
    if '--PE' in sys.argv:
        PE = True
        print "Data are paired end"
    
    if '--STAR' in sys.argv:
        STAR = True
        
    if '--no_introns' in sys.argv:
        no_introns = True
        
    if '--QuantSeq' in sys.argv:
        quant_seq = True
        if '--extend' in sys.argv:
            extend = sys.argv[sys.argv.index('--extend')+1]
        else:
            extend = 200
        
        if '--insert_size' in sys.argv:
            insert_size = sys.argv[sys.argv.index('--insert_size')+1]
        else:
            insert_size = 200
    
    if organism is None:
        organism = 'crypto'
        
    if not STAR:
        print "********Aligning reads with Tophat********\n"
        align_fastq(directory, threads=threads, organism=organism, PE=PE)
        print "********Sorting and indexing BAM files********\n"
        sort_index_bam_files(directory)
    else:
        print "********Aligning reads with STAR********\n"
        align_fastq_STAR(directory, threads=threads, organism=organism, PE=PE, no_introns=no_introns)
        index_STAR_bam(directory)
        
    
    if not quant_seq:
        print "********Counting reads in transcripts with HTseq********\n"
        count_with_HTseq(directory, organism=organism)
    else:
        print "********Counting reads in transcripts with QuantSeq_counting********\n"
        count_for_QuantSeq(directory, organism=organism, threads=threads, extend=extend, insert_size=insert_size)
    
if __name__ == "__main__":
    main()

    
def add_col_level(df, new_name):
    new_col_tuples = []
    for column in df:
        #print column
        new_col_tuples.append((new_name,column))
    new_cols = pd.MultiIndex.from_tuples(new_col_tuples)
    df.columns = new_cols
    return df
    
def load_DESeq2_results(csv_list):
    '''This function loads DESeq2 output csv files into a dataframe.
    
    Parameters
    -----------
    csv_list : list or str
             Either a list of the csv files - ["file1.csv","file2.csv"] or "all".
             If "all", will read all csv files in the current working directory.
    
    Returns
    -------
    pandas.core.frame.DataFrame : a Pandas dataframe with data from all samples'''
    
    if csv_list == 'all':
        csv_list = []
        cwd = os.getcwd()
        for file in os.listdir(cwd):
            if file.endswith('.csv'):
                csv_list.append(file)
        print csv_list
    
    n=0
    for n, file in enumerate(csv_list):
        if n == 0:
            df = pd.read_csv(file, index_col=0)
            df = add_col_level(df, file.split('/')[-1].split('.')[0])
        else:
            new_df = pd.read_csv(file, index_col=0)
            new_df = add_col_level(new_df, file.split('/')[-1].split('.')[0])
            df = df.merge(new_df, left_index=True, right_index=True)
    return df

def load_HTSeq_results(csv_list):
    '''This function loads HTseq output csv files into a dataframe.
    
    Parameters
    -----------
    csv_list : list or str
             Either a list of the csv files - ["file1.csv","file2.csv"] or "all".
             If "all", will read all csv files in the current working directory.
    
    Returns
    -------
    pandas.core.frame.DataFrame : a Pandas dataframe with data from all samples'''
    
    if csv_list == 'all':
        csv_list = []
        cwd = os.getcwd()
        for file in os.listdir(cwd):
            if file.endswith('.htseq'):
                csv_list.append(file)
        print csv_list
    
    n=0
    for n, file in enumerate(csv_list):
        if n == 0:
            df = pd.read_csv(file, index_col=0, names=[file.split('/')[-1].split('.htseq')[0]], sep='\t')
        else:
            new_df = pd.read_csv(file, index_col=0, names=[file.split('/')[-1].split('.htseq')[0]], sep='\t')
            df = df.merge(new_df, left_index=True, right_index=True)
            
    df = df[~df.index.str.contains('__')]
    
    for column in df.columns:
        df[column] = df[column].divide(float(sum(df[column])*1e-6))
        #df[column] = df[column].divide(float(sum(df[column])*1e-6)).apply(np.log2)
    df = df.replace([np.inf,np.inf*-1],np.NaN)
    df = df.dropna(how='any')
    
    return df

def RNAseq_clustered_heatmap(dataframe, sample_names=None, n_clusters=10, use_hdbscan=False, min_cluster_size=15, save_name=None):
    '''Generates a clustered heatmap using kmeans clustering and returns a dataframe with cluster indeces.
    
    Parameters
    ----------
    dataframe : pandas.core.frame.DataFrame
              Return from load_DESeq2_results
    sample_names : list, default ``None``
              Provide a list of sample names. Please use dataframe.columns to ensure the sample names are correctly formatted.
              If no sample names are provided, clusters across all samples in the dataframe.
    n_clusters : int, default 10 
              Number of clusters to create
    use_hdbscan : bool, default ``False``
              If True, will use the hdbscan clustering algorithm instead of K-means. This algorithm is actually much nicer because it does not require you to provide the number of clusters and it does not assume that the clusters are globular. It also does not necessarily place all data points into clusters. However, you do need to provide min_cluster_size, which is the minimum number of points for a group to be considered a cluster.
    min_cluster_size : int, default 15
              See hdbscan variable above for explanation.
    
    Returns
    -------
    clustered_results : pandas.core.frame.DataFrame
              The same dataframe with cluster indeces added as an additional column
              '''
    
    column_names = []
    if sample_names is not None:
        for name in sample_names:
            column_names.append((name,'log2FoldChange'))
    else:
        sample_names = []
        for column in dataframe.columns:
            if column[1] == 'log2FoldChange':
                column_names.append(column)
                sample_names.append(column[0])
        
    data = dataframe[column_names]
    data = data.dropna(how='any')

    # Convert DataFrame to matrix
    mat = data.as_matrix()

    # Cluster data
    if not use_hdbscan:
        km = cluster.KMeans(n_clusters=n_clusters)
        km.fit(mat)
        labels = km.labels_
    else:
        hdb = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
        hdb.fit(mat)
        labels = [x+1 for x in hdb.labels_]

    # Format results as a DataFrame
    data[('cluster','cluster')] = labels
    data = data.sort_values(('cluster','cluster'))
    cluster_sizes = []
    base=0
    for index in set(data[('cluster','cluster')].tolist()):
        cluster_sizes.append(len(data[data[('cluster','cluster')] == index])+base)
        base += len(data[data[('cluster','cluster')] == index])
        
        
    # Get sizes of clusters
    all_data = []
    for column in data[data.columns[:-1]].columns:
        all_data = all_data + data[column].tolist()

    # Plot data as heatmap
    data_min = min(all_data)
    data_max = max(all_data)
    both_max = max(data_min*-1, data_max)
    fig = plt.figure(figsize=(3,12))
    ax = fig.add_subplot(111)
    pcm = ax.pcolor(np.arange(len(data.columns)), np.arange(len(data.index)),
                    data[data.columns[:-1]].as_matrix(), cmap='PuOr_r',
                    norm=colors.Normalize(vmin=both_max*-1, vmax=both_max))
    ax.yaxis.set_ticks(cluster_sizes)
    ax.yaxis.set_ticklabels(range(len(cluster_sizes)))
    plt.xticks([x+0.5 for x in range(len(data.columns))], sample_names, rotation="vertical")
    for size in cluster_sizes:
        ax.plot([0,len(data.columns)-1], [size,size], '-', color='0.5')
        
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(pcm, cax=cax)
    if save_name is not None:
        fig.savefig(save_name+'.pdf', format='pdf')
    plt.show()
    
    new_df = dataframe.merge(pd.DataFrame(data[('cluster','cluster')]), right_index=True, left_index=True)
    
    return new_df

def RNAseq_clustered_heatmap2(dataframe, sample_names=None, n_clusters=10):
    '''Generates a clustered heatmap using kmeans clustering and returns a dataframe with cluster indeces.
    
    Parameters
    ----------
    dataframe : pandas.core.frame.DataFrame
              Return from load_DESeq2_results
    sample_names : list, default ``None``
              Provide a list of sample names. Please use dataframe.columns to ensure the sample names are correctly formatted.
              If no sample names are provided, clusters across all samples in the dataframe.
    n_clusters : int, default 10 
              Number of clusters to create
    
    Returns
    -------
    clustered_results : pandas.core.frame.DataFrame
              The same dataframe with cluster indeces added as an additional column
              '''
    
    if sample_names is None:
        sample_names = dataframe.columns
        
    data = dataframe[sample_names]
    data = data.dropna(how='any')

    # Convert DataFrame to matrix
    mat = data.as_matrix()

    # Using sklearn, cluster data
    km = cluster.KMeans(n_clusters=n_clusters)
    km.fit(mat)

    # Get cluster assignment labels
    labels = km.labels_

    # Format results as a DataFrame
    data['cluster'] = labels
    data = data.sort_values('cluster')
    cluster_sizes = []
    base=0
    for index in set(data['cluster'].tolist()):
        cluster_sizes.append(len(data[data['cluster'] == index])+base)
        base += len(data[data['cluster'] == index])
        
        
    # Get sizes of clusters
    all_data = []
    for column in data[data.columns[:-1]].columns:
        all_data = all_data + data[column].tolist()

    # Plot data as heatmap
    data_min = min(all_data)
    data_max = max(all_data)
    both_max = max(data_min*-1, data_max)
    fig = plt.figure(figsize=(3,12))
    ax = fig.add_subplot(111)
    pcm = ax.pcolor(range(len(data.columns)), range(len(data.index)), data[data.columns[:-1]], cmap='PuOr_r', norm=colors.Normalize(vmin=both_max*-1, vmax=both_max))
    ax.yaxis.set_ticks(cluster_sizes)
    ax.yaxis.set_ticklabels(range(len(cluster_sizes)))
    plt.xticks([x+0.5 for x in range(len(data.columns))], sample_names, rotation="vertical")
    for size in cluster_sizes:
        ax.plot([0,len(data.columns)-1], [size,size], '-', color='0.5')
        
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(pcm, cax=cax)
    plt.show()
    
    new_df = dataframe.merge(pd.DataFrame(data['cluster']), right_index=True, left_index=True)
    
    return new_df

def list_of_genes_in_cluster(data, cluster_index, name=None, annotate=False, organism=None):
    '''Extracts the indicatedcluster from RNAseq_clustered_heatmap and generates a new csv file with/out annotation.
    
    Parameters
    ---------
    data : pandas.core.frame.DataFrame
        Dataframe returned by RNAseq_clustered_heatmap
    cluster_index : int
        The number at the top of the cluster on the heatmap display
    name : str, default ``None``
        Name for csv file
    annotate : bool, default ``False``
        Choose whether or not to annotate the csv file
    organism : str, default ``None``
        Required if annotate is True. "pombe" or "crypto" to apply the correct annotation
    
    Output
    ------
    csv_file : csv formatted file with genes in cluster
        '''
    
    in_cluster = data[data[('cluster','cluster')] == cluster_index]
    
    if name is None:
        name = 'cluster'+str(cluster_index)
    
    multi_index = False
    try:
        if type(in_cluster.columns) == pd.core.indexes.multi.MultiIndex: 
            multi_index = True
    except AttributeError:
        if type(in_cluster.columns) == pd.indexes.multi.MultiIndex:
            multi_index = True
    
    in_cluster.to_csv(name+'_genes.csv')
            
    if annotate is True:
        if organism is None:
            print "Must include the organism!"
        elif 'crypto' in organism.lower():
            Annotation_tools.crypto_annotate(name+'_genes.csv', multi_index=multi_index)
        elif 'pombe' in organism.lower():
            Annotation_tools.pombe_annotate(name+'_genes.csv', multi_index=multi_index)
        else:
            print "Organism not recognized. Only 'pombe' and 'crypto' supported at this time"

def volcano_plot(df, sample_names=None, annotate=False, organism=None, change_thresh=1, padj_thresh=0.01):
    '''Generates volcano plots and csv files of significantly changed transcripts for each sample.
    
    Parameters
    ----------
    df : pandas.core.frame.DataFrame
         Output from load_DESeq2_results
    sample_names : list, default ``None``
         List of names if you don't want to analyze all the samples in the dataframe
    annotate : bool, default ``False``
        True if you want to generate a csv file with annotation
    organism : str, default ``None``
        Required if annotate is True. "pombe" or "crypto" to apply the correct annotation
        
    Output
    ------
    csv files : 3 comma separated value files
              1. All significantly changed transcripts
              2. Transcripts with increased expression
              3. Transcripts with decreased expression
              **Note if you use the annotate option, there will be 3 more files with annotation.
    volcano plots : PDF formated images of the volcano plots
        '''
    
    # Find all sample names
    if sample_names is None:
        sample_names = []
        for column in df.columns:
            if column[1] == 'log2FoldChange':
                sample_names.append(column[0])

    #df = df.dropna(how='any')
    for name in sample_names:
        ## Find significantly changed transcripts
        sig = df[df[(name,'padj')] <= padj_thresh][name]
        sig = sig[sig['log2FoldChange'] >= change_thresh]
        sig2 = df[df[(name,'padj')] <= padj_thresh][name]
        sig2 = sig2[sig2['log2FoldChange'] <= -1*change_thresh]
        print name
        print "Increased:"
        print len(sig)
        print "Decreased:"
        print len(sig2)
        
        #Build plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(df[(name,'log2FoldChange')], df[(name,'padj')].apply(np.log10).multiply(-1), color='0.5', s=15)
        ax.scatter(sig['log2FoldChange'], sig['padj'].apply(np.log10).multiply(-1), color='orange', s=15)
        ax.scatter(sig2['log2FoldChange'], sig2['padj'].apply(np.log10).multiply(-1), color='darkslateblue', s=15)
        xmin = (min(df[(name,'log2FoldChange')].replace([np.inf, np.inf*-1], np.NaN).dropna(how='any'))-
		0.1*min(df[(name,'log2FoldChange')].replace([np.inf, np.inf*-1], np.NaN).dropna(how='any')))
        xmax = (max(df[(name,'log2FoldChange')].replace([np.inf, np.inf*-1], np.NaN).dropna(how='any'))+
		0.1*max(df[(name,'log2FoldChange')].replace([np.inf, np.inf*-1], np.NaN).dropna(how='any')))
        x_all = max([xmin*-1,xmax])
        print x_all
        xmin= -1*x_all
        xmax= x_all
        
        ymed = np.percentile(df[(name,'padj')].apply(np.log10).multiply(-1).replace([np.inf, np.inf*-1], np.NaN).dropna(how='any'), 99)
	ymax = ymed + 1.5*ymed
        
        if ymax < 4:
            ymax = 20
        
        ax.set_ylim([-5,ymax])
        ax.set_xlim([xmin,xmax])
        ax.plot([xmin,xmax],[np.log10(padj_thresh)*-1,np.log10(padj_thresh)*-1], '--', color='0.7')
        ax.plot([change_thresh,change_thresh], [-5,ymax], '--', color='0.7')
        ax.plot([-1*change_thresh,-1*change_thresh], [-5,ymax], '--', color='0.7')
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 padj")
        ax.set_title(name)
        plt.show()

        fig.savefig(name+'_volcano.pdf',format='pdf')
        plt.clf()
        
        new_sig = sig.append(sig2)
        print len(new_sig)
        new_sig.to_csv(name+'_sig_changed.csv')
        split_DE_genes(name+'_sig_changed.csv')
        
        if annotate is True:
            if organism is None:
                print "Must include the organism!"
            elif 'crypto' in organism.lower():
                Annotation_tools.crypto_annotate(name+'_sig_changed.csv')
                Annotation_tools.crypto_annotate(name+'_increased.csv')
                Annotation_tools.crypto_annotate(name+'_decreased.csv')
            elif 'pombe' in organism.lower():
                Annotation_tools.pombe_annotate(name+'_sig_changed.csv')
                Annotation_tools.pombe_annotate(name+'_increased.csv')
                Annotation_tools.pombe_annotate(name+'_decreased.csv')
            else:
                print "Organism not recognized. Only 'pombe' and 'crypto' supported at this time"

def split_DE_genes(csv):
    csv_df = pd.read_csv(csv, index_col=0)
    
    #Split dataframe
    down_df = csv_df[csv_df['log2FoldChange'] < 0]
    up_df = csv_df[csv_df['log2FoldChange'] > 0]
    
    down_df.to_csv(csv.split('sig_changed.csv')[0]+'decreased.csv')
    up_df.to_csv(csv.split('sig_changed.csv')[0]+'increased.csv')
                             
def NchooseR(n,r):
    ''' Called by gene_venn'''
    c = factorial(n)/Decimal(factorial(r)*factorial(n-r))
    return c

def hypergeometric(N, n, K, J, k):
    #print "N: "+str(N)
    #print "n: "+str(n)
    #print "K: "+str(max(K,J))
    #print "k: "+str(k)
    
    '''Called by gene_venn'''
    K = max(K,J)
    if N-K < 0:
        print "Too many genes in your spreadsheet! Check your formatting."
        return None, None
    else:
        p=0
        for i in range(k,min([n,K])):
            p += (NchooseR(K,i)*NchooseR(N-K,n-i))/NchooseR(N,n)
    return p

def venn_2sample(n,K,k,J, name1, name2, colors, p):
    '''Called by gene_venn'''
    A_only = K-k
    B_only = J-k
    AB = k
    s = (A_only, B_only, AB)

    fig, ax = plt.subplots(1)
    v = venn2(subsets=s, set_labels=(name1, name2), ax=ax)

    # Subset colors
    v.get_patch_by_id('10').set_color(colors[0])
    v.get_patch_by_id('01').set_color(colors[1])
    if AB > 0:
        v.get_patch_by_id('11').set_color(colors[2])

    # Subset alphas
    v.get_patch_by_id('10').set_alpha(0.8)
    v.get_patch_by_id('01').set_alpha(0.8)
    if AB > 0:
        v.get_patch_by_id('11').set_alpha(0.9)

    # Border styles
    c = venn2_circles(subsets=s, linestyle='solid', ax=ax)
    for sub_c in c:
        sub_c.set_lw(1.0)       # Line width

    p = '%.1E' % p
    
    ax.set_title(name1+' vs. '+name2+'\n p-value: '+p, fontdict=font)
    fig.savefig(name1+'_'+name2+'_venn.pdf', format='pdf', bbox_inches='tight')
    plt.show()
    plt.clf()
    
def gene_venn(csv_files, organism):
    '''Finds overlap between 2 or 3 lists of genes.
    
    Parameters
    ----------
    csv_files : list
               2 or 3 csv files where the first column is the gene name (make sure the gene name format matches).
    organism : str
               Options are 'crypto', 'cerevisiae' or 'pombe'
    
    Output
    ------
    PDF files of venn diagrams (pairwise) and merged csv files containing the overlapping genes.'''
        
    if 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        organism = 'pombe'
    elif 'crypto' in organism.lower() or 'h99' in organism.lower():
        organism = None
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    elif 'cerev' in organism.lower():
        organism = None
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
    
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    transcripts = tx_dict.keys()
    genes = set([x[:-2] for x in transcripts])
    
    df_dict = {}
    names = []
    for csv in csv_files:
        name = csv.split('/')[-1].split('.')[0]
        names.append(name)
        line = next(open(csv))
        if len(line.split(',')) > 1:
            df_dict[name] = pd.read_csv(csv, index_col=0)
            df_dict[name] = add_col_level(df_dict[name], name)
        elif len(line.split('\t')) > 1:
            df_dict[name] = pd.read_csv(csv, index_col=0, sep='\t')
            df_dict[name] = add_col_level(df_dict[name], name)
        else:
            df_dict[name] = pd.read_csv(csv, index_col=0)
    
    # N = genome size
    # n = number of genes in analysis (so len(a)+len(b))
    # K = number of genes in group 1 (len(a))
    # k = overlap(len(a&b))
    
    N = len(genes)
    K = len(df_dict[names[0]])
    overlap = set(df_dict[names[0]].index).intersection(df_dict[names[1]].index)
    k = len(overlap)
    n = len(df_dict[names[0]]) + len(df_dict[names[1]]) - k
    J = len(df_dict[names[1]])
    
    p_ab= hypergeometric(N,n,K,J,k)
    if p_ab is not None:
        venn_2sample(n,K,k,J, names[0], names[1], ['crimson','deepskyblue','darkorchid'], p_ab)
        df_ab = df_dict[names[0]].merge(df_dict[names[1]], right_index=True, left_index=True)
        df_ab.to_csv('{0}_{1}_overlap.csv'.format(names[0], names[1]))

    if len(names) == 3:
        ## Compare sample 1 to sample 3
        overlap_ac = set(df_dict[names[0]].index).intersection(df_dict[names[2]].index)
        k_ac = len(overlap_ac)
        n_ac = len(df_dict[names[0]]) + len(df_dict[names[2]]) - k_ac
        J_ac = len(df_dict[names[2]])
        
        p_ac = hypergeometric(N,n_ac,K,J_ac,k_ac)
        if p_ac is not None:
            venn_2sample(n_ac,K,k_ac,J_ac, names[0], names[2], ['crimson','gold','darkorange'], p_ac)
            df_ac = df_dict[names[0]].merge(df_dict[names[2]], right_index=True, left_index=True)
            df_ac.to_csv('{0}_{1}_overlap.csv'.format(names[0], names[2]))
            
        ## Compare sample 2 to sample 3
        overlap_bc = set(df_dict[names[1]].index).intersection(df_dict[names[2]].index)
        k_bc = len(overlap_bc)
        n_bc = len(df_dict[names[1]]) + len(df_dict[names[2]]) - k_bc
        J_bc = len(df_dict[names[2]])
        K_bc = len(df_dict[names[1]])
        
        p_bc = hypergeometric(N,n_bc,K_bc,J_bc,k_bc)
        if p_bc is not None:
            venn_2sample(n_bc,K_bc,k_bc,J_bc, names[1], names[2], ['deepskyblue','gold','forestgreen'], p_bc)
            df_bc = df_dict[names[1]].merge(df_dict[names[2]], right_index=True, left_index=True)
            df_bc.to_csv('{0}_{1}_overlap.csv'.format(names[1], names[2]))
            
def map_to_chromosomes(csv, organism, fig_name="chrom_map"):
    '''Maps a list of transcripts onto chromosomes. If log2FoldChange column is in csv file, 
    will color code by increased or decreased.
    
    Parameters
    ----------
    csv : str
         File in csv format. Can also just be a list of genes in a text file.
    organism : str
         Options are 'crypto' or 'pombe'
    fig_name : str, default "chrom_map"
         Provide a name for the output file.
    
    Output
    ------
    PDF file of chromosome diagrams'''
    
    df = pd.read_csv(csv, index_col=0)
    if df.index[0].startswith('gene'):
        df.index = [x.split('gene:')[1] for x in df.index]
    if len(df.columns) <= 1:
	df = pd.read_csv(csv, sep='\t', index_col=0)
    
    if organism == 'crypto': 
        fa_json = '/home/jordan/GENOMES/H99_fa.json'
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
        cen_dict = '/home/jordan/GENOMES/H99_centromeres.json'
    elif organism == 'pombe': 
        fa_json = '/home/jordan/GENOMES/POMBE/Sp_fasta_dict.json'
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        cen_dict = '/home/jordan/GENOMES/POMBE/Sp_centromeres.json'
    #elif 'cerev' in organism.lower(): 
    #    fa_json = '/home/jordan/GENOMES/S288C/S288C_genome.fa'
    #    gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        
    with open(fa_json) as f:
        fa_dict = json.load(f)
        
    with open(cen_dict) as f:
        cen_dict = json.load(f)

    chrom_sizes = {k:len(v) for k, v in fa_dict.items() if k.startswith('chr')}
    chrom_sizes = OrderedDict(sorted(chrom_sizes.items(), key=lambda t: t[0]))
    
    if organism != 'pombe': 
        organism = None
        figsize = (10,12)
    else:
        figsize = (15,3)
    tx_dict = SeqTools.build_transcript_dict(gff3, organism=organism)
    
    tx_by_chrom = {k:set() for k, v in fa_dict.items()}
    for tx, info in tx_dict.iteritems():
        tx_by_chrom[info[3]].add(tx[:-2])
        
    divis = max(chrom_sizes.values())/0.85
    #print divis
    
    f, ax = plt.subplots(len(chrom_sizes), figsize=figsize, sharex=True, sharey=True)
    for n, chrom in enumerate(chrom_sizes.keys()):
        ax[n].add_patch(patches.Rectangle((0.1, 0.1), chrom_sizes[chrom]/divis, 0.6, fill=False, linewidth=1.5))
        ax[n].set_xticks([])
        ax[n].set_yticks([0.5])
        ax[n].axis('off')
        #ax[n].set_title(chrom, loc='left')
        
        ax[n].text(0.08, 0.4, chrom, horizontalalignment='right', verticalalignment='center', transform=ax[n].transAxes)
        
        chrom_df = df[df.index.isin(tx_by_chrom[chrom])]
        
        starts = []
        stops = []
        for gene in chrom_df.index:
            tx = gene + 'T0'
            if organism == 'pombe': 
                tx = gene+'.1'
            
            starts.append(tx_dict[tx][0])
            if tx_dict[tx][1] <= chrom_sizes[chrom]:
                stops.append(tx_dict[tx][1])
            else:
                stops.append(chrom_sizes[chrom])
        
        chrom_df.loc[:,'start'] = starts
        chrom_df.loc[:,'stop'] = stops
        chrom_df = chrom_df.sort_values('start')
        
        chrom_patches = []
        
        #Draw centromere
        cen_start = 0.1+cen_dict[chrom][0]/divis
        cen_width = (cen_dict[chrom][1]-cen_dict[chrom][0])/divis
        ax[n].add_patch(patches.Rectangle((cen_start, 0.1), cen_width, 0.6, color='0.1'))
        
        #Draw affected transcripts
        for gene, r in chrom_df.iterrows():
            if 'log2FoldChange' in chrom_df.columns:
                if r['log2FoldChange'] > 0:
                    color = 'orange'
                else:
                    color = 'slateblue'
            else:
                color = 'red'
            x = 0.1+r['start']/divis
            width = (r['stop']-r['start'])/divis
            chrom_patches.append((patches.Rectangle((x, 0.13), width, 0.57, color=color, linewidth=1, alpha=0.6)))
       
        for p in chrom_patches:
            ax[n].add_patch(p)
            
    f.savefig(fig_name+'.pdf', format='pdf')

    plt.show()
    plt.clf()
    
def RNAseq_log2foldchange_scatter(csv1, csv2, gene_list=None, sep=',', color='0.3'):
    '''Creates scatter of log2FoldChange from 2 DESeq2 output files. Provide a gene_list (can be a spreadsheet) if you want to filter which points are shown. Each point is a transcript.
    
    Parameters
    ----------
    csv1 : str
         DESeq2 output file in csv format.
    csv2 : str
         DESeq2 output file in csv format.
    gene_list : str or list, default ``None``
         Some kind of text file with the genes you want to plot as the first column. Can also be a Python list.
    sep : str, default ","
         Type of separator in the csv1 and csv2 files
    color : str, default "0.3", which is dark grey
         Color of points in scatter plot. Look up acceptable names for colors in matplotlib.
    
    Returns
    ------
    fig : matplotlib figure object
          fig can be saved using the command fig.savefig("somename.pdf", format="pdf", bbox_inches="tight")'''
    
    df1 = pd.read_csv(csv1, sep=sep, index_col=0)
    df2 = pd.read_csv(csv2, sep=sep, index_col=0)
    
    if type(gene_list) == str:
        gene_list_open = []
        with open(gene_list) as f:
            for line in f:
                gene_list_open.append(line.split(' ')[0].split(',')[0].split('\t')[0].strip())
        gene_list = gene_list_open
    
    df = df1.merge(df2, right_index=True, left_index=True, how='inner')
    if gene_list is not None:
        df = df[df.index.isin(gene_list)]
    
    fig, ax = plt.subplots(1, figsize=(5,5))
    
    cols = [x for x in df.columns if 'log2FoldChange' in x]
    
    ax.scatter(df[cols[0]],df[cols[1]], s=15, color=color, alpha=0.5)
    r, p = stats.pearsonr(df[cols[0]],df[cols[1]])
    print "Pearson R: " + "%.2f" % r
    print "p-value: " + "%.2E" % p
    
    min_xy = min(ax.get_xlim()[0], ax.get_ylim()[0])
    max_xy = max(ax.get_xlim()[1], ax.get_ylim()[1])
    
    ax.set_xlim(min_xy, max_xy)
    ax.set_ylim(min_xy, max_xy)
    ax.plot(ax.get_xlim(), ax.get_ylim(), '--', color='0.7', zorder=0)
    
    ax.set_xlabel(csv1.split('.csv')[0].split('/')[-1]+'\nlog2FoldChange', fontsize=14)
    ax.set_ylabel(csv2.split('.csv')[0].split('/')[-1]+'\nlog2FoldChange', fontsize=14)
    
    plt.show()
    return fig
