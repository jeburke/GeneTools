# Scripts
## RNAseq_tools.py
Starting with fastq.gz files, aligns to genome with Tophat (using Bowtie 1), generates sorted, indexed bam files and counts reads in annotated transcripts using htseq-count. Generates all files necessary to proceed with DESeq2 analysis.
### Usage:  
```Usage: python RNAseq_tools.py --directory directory --threads num_tophat_threads --organism organism_name <--PE> <--STAR>```
  
Required arguments:  
>--directory : directory containing fastq.gz files (e.g. ./ or /data/jordan/JEB001/)  
>--threads : processors to be used when running Tophat  
>--organism : yeast genus or species - currently working for crypto/pombe/cerevisiae for TopHat and crypto/pombe for STAR
  
Optional arguments:   
>--PE : Include if paired-end data
>--STAR : Align with STAR instead of Tophat/Bowtie2
  
## ChIP_tools.py
Starting with fastq.gz files, trims adaptor with Cutadapt, aligns to genome with Bowtie 1 and generates sorted, indexed bam files. See log files for Cutadapt and Bowtie summaries.
### Usage:
```python ChIP_tools.py --directory fastq_directory --threads num_threads --organism crypto/pombe/cerevisiae/candida <--adaptor GATCGGAAGA> <--gff3 gff3_file> <--index bowtie_index_prefix>```  
  
Required positional arguments:  
>--directory : directory with fastq.gz files (e.g. ./ or /data/jordan/JEB001/)  
>--threads : number of processors to use  
>--organism : yeast genus or species - currently can be crypto/pombe/cerevisiae/candida  

Optional arguments:  
>--adaptor : change only if library preparation differs from typical Madhani lab ChIP-seq protocol  
>--gff3 and --index : can provide custom gff3 file and bowtie index i.e. for a different organism or updated versions of the genome. Must provide both.  

## RNAi_tools.py
Counts reads in sorted, indexed bam files for all annotated transcripts on sense and antisense strand and generates a spreadsheet. Also compares two genotypes/conditions for changes in siRNA abundance and generates histograms of the RNA fragment sizes for each provided bam file. Designed for experiments run in duplicate only for now.
### Usage:
```python RNAiTools.py [bam_files] output_name```  
  
Required positional arguments:
>bam_files : Indexed, sorted bam file locations - must be 4 or 8 files with control replicates first (wild type or untreated).  
>output_name : Prefix for naming output spreadsheet and figures (must be last argument).  
  
## Generate_bedgraphs.py
Creates bedgraph files for every indexed, sorted bam file in a directory. Bedgraphs are scaled based on the number of aligned reads in the bam file (final bedgraph scale is in reads per million (RPM)). Can also normalize to a whole cell extract or untagged sample and smooth by applying a rolling mean.
### Usage: 
```python Generate_bedgraphs.py directory organism <--start_only> <--stranded> <--normalize untagged_sample_name> <--smooth window>```  
  
Required positional arguments:  
>directory : directory containing bam files  
>organism : crypto, pombe, candida or cerevisiae  
  
Optional arguments:  
>--start-only : the bedgraph will include only the read start positions (good for NET-seq and Spliceoosome 3' end profiling)  
>--stranded : generate separate bedgraphs for Watson and Crick strands  
>--normalize : generate an additional bedgraph that is divided by a background sample, like whole cell extract. Provide the name of the background sample  
>--smooth : generate smoothed bedgraph files over the indicated window in base-pairs (note: this will probably be half what you normally use as a window for C. Homer's script)  

# Functions
Access these function by importing the GeneTools module. Must set up a display environment for matplotlib (works best with %matplotlib inline in a Jupyter notebook) 
  
## Functions for RNA-seq analysis
#### load_DESeq2_results(csv_list)  
    '''This function loads DESeq2 output csv files into a dataframe.
    
    Parameters
    -----------
    csv_list : list or str
             Either a list of the csv files - ["file1.csv","file2.csv"] or "all".
             If "all", will read all csv files in the current working directory.
    
    Returns
    -------
    pandas.core.frame.DataFrame : a Pandas dataframe with data from all samples'''
#### load_HTSeq_results(csv_list)
    '''This function loads HTseq output csv files into a dataframe.
    
    Parameters
    -----------
    csv_list : list or str
             Either a list of the csv files - ["file1.csv","file2.csv"] or "all".
             If "all", will read all csv files in the current working directory.
    
    Returns
    -------
    pandas.core.frame.DataFrame : a Pandas dataframe with data from all samples'''
#### RNAseq_clustered_heatmap(dataframe, sample_names=None, n_clusters=10)
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
              The same dataframe with cluster indeces added as an additional column'''
#### list_of_genes_in_cluster(data, cluster_index, name=None, annotate=False, organism=None)
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
    csv_file : csv formatted file with genes in cluster'''
 #### volcano_plot(df, sample_names=None, annotate=False, organism=None, change_thresh=1, padj_thresh=0.01)
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
    volcano plots : PDF formated images of the volcano plots'''
#### gene_venn(csv_files, organism)
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
#### map_to_chromosomes(csv, organism, fig_name="chrom_map")
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
#### RNAseq_log2foldchange_scatter(csv1, csv2, gene_list=None, sep=',', color='0.3'):
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
  
## Functions for ChIP-seq analysis  
#### ChIP_rpkm_scatter(WCE_bam, WT1_bam, WT2_bam, Mut1_bam, Mut2_bam, gff3, plot_name, Z_change=False, cen_tel=False):
    '''Plots RPKM as scatter plots from replicates of two different genotypes/conditions - can do promoters or just centromeres and telomeres
    
    Parameters
    ----------
    WCE_bam : str
            Whole cell extract bam file
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
#### MACS_peak_RPKM_scatters(xls_pair1, xls_pair2, untagged_xls, bam_list, WCE_bam_list, name, organism='crypto', min_overlap=0.5, exclude_chrom=None)
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
            
    Outputs
    ------
    Comparison csv files : csv files containing the reproducible, reduced set of MACS peaks
    Genes with peaks : files with a list of genes that had at least one reproducible peak
    Scatter spreadsheet : The data the corresponds to the scatter plots - RPKM of each peak normalized to whole cell extract
    Scatter plot : pdf file with scatter plot'''
#### get_peak_sequence(input_file, fasta_file, gff3_file, window=1000):
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
      
## Other functions
#### igv_plots_general(bam_list, gene_list, organism, colors=None, names=None, save_dir=None, unstranded=False, end_only=False, same_yaxis=False, specific_range=None, transcript_direction=True, log_scale=False, rpm=True, PE=False, plot_junctions=False)
    '''Generates figures of IGV-like plots using pyplot of reads from indexed, sorted bam files.
    
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
            directory to save eps files (will be created if it does not yet exist). If None, does not save files
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
            
    Outputs
    -------
    figures : EPS format files in indicated save_dir, vector graphic that can be edited with Illustrator'''
#### crypto_annotate(text_file, sep=',', multi_index=False)
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
    csv file : File (comma separated) with annotation appended to each row.'''
#### pombe_annotate(text_file, sep=',', multi_index=False)
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
    csv file : File (comma separated) with annotation appended to each row.'''
#### count_reads_in_window(bam, chrom, start, end, strand)
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
#### count_aligned_reads(bam_file)
    '''Counts aligned reads in bam file reads using Samtools
    
    Parameters
    ----------
    bam_file : str
            bam file from Bowtie or STAR
    
    Returns
    ------
    total : float
         Million aligned reads'''
#### PE_fragment_size(bam_file)
    '''Calculates average and standard deviation of insert fragment sizes from paired end data. Necessary for GEO deposition
    
    Parameters
    ----------
    bam_file : str
            bam file from Bowtie or STAR from paired end data
    
    Output
    ------
    Prints the average and standard deviation of the library fragment size'''
#### generate_read_series(bam_iterator, chrom, start, end, strand, baseline=0)
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
