import os
import subprocess
import argparse

def run_cutadapt(directory, adaptor="TGGAATTCTC", threads=10):
    print('Trimming reads...\n')
    processes = []
    for fq in [x for x in os.listdir(directory) if x.endswith("fastq.gz")]:
        prefix = fq.split('/')[-1].split('_S')[0]
        print prefix
        
        with open(directory+prefix+'_cutadapt.log','w') as logfile:
            cutadapt_args = "cutadapt -a {0} -m 16 -o {1}_trim.fastq {2}".format(adaptor, directory+prefix, directory+fq)
            if prefix+'_trim.fastq' not in os.listdir(directory):
                processes.append(subprocess.Popen(cutadapt_args, shell=True, universal_newlines=True, stdout=logfile, stderr=logfile))
                
            if len(processes) == threads:
                processes[0].wait()
                    
    for p in processes:
        p.wait()
        
def run(cmd, logfile):
    '''Function to open subprocess, wait until it finishes and write all output to the logfile'''
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile, stderr=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

def run_bowtie(directory, threads=10):
    print('Aligning with Bowtie1...\n')
    with open('bowtie.log','w') as logfile:
        for fq in [x for x in os.listdir(directory) if x.endswith("trim.fastq")]:
            prefix = fq.split('/')[-1].split('_trim')[0]
            print prefix
            bowtie_args = "bowtie -p{0} -v2 -M1 --best --max {1}_multi --un {1}_un.fastq /home/jordan/GENOMES/Crypto_for_gobs -q {2} --sam {1}.sam".format(threads, directory+prefix, directory+fq)
            if prefix+'.sam' not in os.listdir(directory):
                logfile.write("\n***"+prefix+"***\n")
                ret_code = run(bowtie_args, logfile)
            else:
                ret_code = None
    return ret_code

def sort_index(directory, threads=10):
    print('Converting to bam format...\n')
    processes = []
    for sam in [x for x in os.listdir(directory) if x .endswith('.sam')]:
        name = directory+sam.split('.sam')[0]
        args = "samtools view -Sbo {0}.bam {0}.sam".format(name)
        processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
        
        if len(processes) == threads:
            processes[0].wait()
                
    for p in processes: p.wait()

    print('Sorting...\n')
    processes = []
    for bam in [x for x in os.listdir(directory) if x .endswith('.bam')]:
        name = directory+bam.split('.bam')[0]
        args = "samtools sort {0}.bam -o {0}_sorted.bam".format(name)
        processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
        
        if len(processes) == threads:
            processes[0].wait()
            
    for p in processes: p.wait()

    print('Indexing...\n')
    processes = []
    for bam in [x for x in os.listdir(directory) if x .endswith('_sorted.bam')]:
        args = "samtools index {0}".format(directory+bam)
        processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
          
        if len(processes) == threads:
            processes[0].wait()
          
    for p in processes: p.wait()
          
def count_antisense_reads(directory, gff3, threads=10):
    print('Counting antisense reads...\n')
    processes = []
    for bam in [x for x in os.listdir(directory) if x.endswith('_sorted.bam')]:
        out_name = bam.split('_sorted.bam')[0].split('/')[-1]
        if out_name+'.quant' not in os.listdir(directory):
            args = "python /home/jordan/CodeBase/RNA-is-awesome/QuantSeq_counting.py --bam_file {0} --gff3 {1} --name {2}".format(directory+bam,gff3,directory+out_name)
            processes.append(subprocess.Popen(args.split(' ')))
        
        if len(processes) == threads:
            processes[0].wait()
          
    for p in processes:
        p.wait()
        
    for file in os.listdir(directory):
        if file.endswith('.sam'): os.remove(directory+file)
        elif file.endswith('.bam') and 'sorted' not in file: os.remove(directory+file)
        elif file.endswith('trim.fastq'): os.remove(directory+file)
            
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--directory", default='./', help="Working directory containing fastq files")
    parser.add_argument("--adaptor", default="TGGAATTCTC", help="Adaptor sequence to trim from 3' end of reads")
    parser.add_argument("--gff3", default="/home/jordan/GENOMES/siRNA_all.gff3", help="GFF3 file for read counting")
    parser.add_argument("--threads", default=10, help="Number of processors to use for analysis")
    parser.add_argument("--count_only", action='store_true', help="Only count reads, do not run alignment and formatting")
    args = parser.parse_args()

    if not args.directory.endswith('/'):
        args.directory = args.directory+'/'
        
    if not args.count_only:
        run_cutadapt(args.directory, adaptor=args.adaptor, threads=args.threads)
        code = run_bowtie(args.directory, threads=args.threads)
        sort_index(args.directory, threads=args.threads)
    count_antisense_reads(args.directory, args.gff3, threads=args.threads)
    
if __name__ == "__main__":
    main()
