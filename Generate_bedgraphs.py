'''Usage: python Generate_bedgraphs.py directory organism <--start_only> <--stranded> <--normalize untagged_sample_name> <--smooth window>
Arguments in <> are optional
Organism can be crypto, pombe or cerevisiae
Include the start_only argument to map only the 5' ends of reads'''

import sys
import os
import warnings; warnings.simplefilter('ignore')
import argparse
script_path = os.path.dirname(os.path.realpath(__file__)).split('GeneTools')[0]
sys.path.append(script_path)
import GeneTools as GT
from multiprocessing import Pool

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument("directory", default='./', help="Working directory containing fastq files")    
    parser.add_argument("organism", default="crypto", help="Organisms available: crypto, pombe, cerevisiae, candida")    
    parser.add_argument("--threads", default=1, type=int, help="Number of processors")   
    parser.add_argument("--start_only", action='store_true', help="Include only the start of the read in bedgraph")    
    parser.add_argument("--stranded", action="store_true", help="Create separate bedgraphs for W and C strands")
    parser.add_argument("--subtract", action="store_true", help="Subtract background levels using a rolling window")
    parser.add_argument("--normalize", default=None, help='Normalize to provided bam file name')
    parser.add_argument("--smooth", type=int, default=0, help='Normalize to provided bam file name')
    parser.add_argument("--bam_names", default=None, nargs='+', help='Specify bam names')
    args = parser.parse_args()
                        
    directory = args.directory
    if not directory.endswith('/'):
        directory = directory+'/'
        
    
    normalize = False
    if args.normalize is not None:
        normalize = True
        untagged = args.normalize
    
    if not untagged.endswith('.bam'):
        try:
            untagged_bam = [directory+x for x in os.listdir(directory) if untagged in x and x.endswith('.bam')][0]
        except IndexError:
            "Can't find background file for normalization... aborting."
            return None
        if args.bam_names is not None:
            args.bam_names.append(untagged_bam)
    
    smooth = False
    if args.smooth != 0:
        smooth = True
        window = args.smooth

    if 'crypto' not in args.organism.lower() and 'pombe' not in args.organism.lower() and 'candida' not in args.organism.lower() and 'cerev' not in args.organism.lower():
        try:
            with open(organism) as f:
                for line in f:
                    continue
        except IOError:
            print "Unrecognized organism"
            return None 

    expand = False
    if normalize  or smooth:
        expand = True
    
    print "Generating scaled bedgraphs..."
    GT.generate_scaled_bedgraphs2(directory, untagged, organism=args.organism, start_only=args.start_only, stranded=args.stranded, threads=args.threads, expand=expand, bam_list=args.bam_names)


    if normalize is True:
        print "\nNormalizing to untagged..."
        if untagged.endswith('.bam'):
            untagged = untagged.split('/')[-1].split('.bam')[0]

        bg_list = [directory+x for x in os.listdir(directory) if x.endswith('.bedgraph')]
        bg_list = [x for x in bg_list if 'norm' not in x and 'smooth' not in x]
        try:
            untagged_bg = [x for x in bg_list if untagged in x][0]
            bg_list.remove(untagged_bg)
        except IndexError:
            untagged_bg = [directory+x for x in os.listdir(directory) if untagged in x and x.endswith('.bedgraph')]
            if len(untagged_bg) > 1:
                untagged_bg = [x for x in untagged_bg if 'smooth' not in x and 'norm' not in x]
            elif len(untagged_bg) == 1: 
                untagged_bg = untagged_bg[0]
            elif len(untagged_bg) == 0:
                untagged_bg = untagged.split('/')[-1].split('.bam')[0]+'.bedgraph'
        
        last = False
        for n, bg in enumerate(bg_list):
            if n == len(bg_list)-1:
                last = True
            GT.normalize_bedgraph(bg, untagged_bg, smooth=smooth, last=last)

    if smooth is True:
        print "\nSmoothing with {0} bp window...".format(str(window))
        bg_list = [directory+x for x in os.listdir(directory) if x.endswith('.bedgraph')]
        bg_list = [x for x in bg_list if 'smooth' not in x]
        GT.smooth_bedgraphs(bg_list, window)
        
    if args.subtract:
        print "\nSubtracting background..."
        bg_list = [directory+x for x in os.listdir(directory) if x.endswith('.bedgraph')]
        p = Pool(threads/2)
        p.map(GT.background_subtraction, bg_list)
        
if __name__ == "__main__":
    main()
