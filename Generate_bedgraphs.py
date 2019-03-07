'''Usage: python Generate_bedgraphs.py directory organism <--start_only> <--stranded> <--normalize untagged_sample_name> <--smooth window>
Arguments in <> are optional
Organism can be crypto, pombe or cerevisiae
Include the start_only argument to map only the 5' ends of reads'''

import sys
script_path = os.path.dirname(os.path.realpath(__file__)).split('GeneTools')[0]
sys.path.append(script_path)
import GeneTools as GT
import os
from multiprocessing import Pool

def main():
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print '\nUsage:\n python Generate_bedgraphs.py directory organism <--threads n> <--start_only> <--stranded> <--normalize untagged_sample_name> <--smooth window>'
        print 'Arguments in <> are optional'
        print "directory : ./ for current directory or specify. Sorted bam files must be in specified directory. All files wil be output to specified directory."
        print 'organism : crypto, pombe, cerevisiae or candida'
        print "--threads : number of processors to use, default=1"
        print "--start_only : Map only the 5' ends of reads"
        print "--stranded : Create separate bedgraphs for the Watson and Crick strands.\n   Disclaimer - this also produces a combined bedgraph where reads on the Crick strand are negative. Please note that areas with overlapping transcription will cancel out in this representation."
        print "--normalize : Normalize to an untagged or whole cell extract sample (untagged_sample_name)"
        print "--smooth : Smooth data using rolling mean. Provide window in base pairs (e.g. 100)\n"
        return None
    
    directory = sys.argv[1]
    file_provided = False
    if directory.endswith('.bam'):
        file_provided = True
        print "File provided"
    elif not directory.endswith('/'):
        directory = directory+'/'
    
    organism = sys.argv[2]
    if 'crypto' not in organism.lower() and 'pombe' not in organism.lower() and 'candida' not in organism.lower() and 'cerev' not in organism.lower():
        try:
            with open(organism) as f:
                for line in f:
                    continue
        except IOError:
            print "Unrecognized organism"
            return None 

    threads=1
    start_only = False
    stranded = False
    normalize = False
    smooth = False
    untagged = None
    window = None
    sub = False
    for n, arg in enumerate(sys.argv):
        if arg == "--threads":
            try:
                threads = int(sys.argv[n+1])
            except IndexError:
                "Must provide number of threads"
                return None
            except ValueError:
                "Threads must be an integer (or number of threads not provided)"
                return None
        elif arg == "--start_only":
            start_only = True
        elif arg == "--stranded":
            stranded = True
        elif arg == "--normalize":
            normalize = True
            try:
                untagged = sys.argv[n+1]
            except IndexError:
                print "Must provide untagged sample name"
                return None
            
            if untagged.startswith('--'):
                print "Must provide untagged sample name"
                return None
        elif arg == "--smooth":
            smooth = True
            try:
                window = int(sys.argv[n+1])
            except IndexError:
                print "Must provide window size"
                return None
            except ValueError:
                print "Must provide window size"
                return None
        elif arg == "--subtract_background":
            sub = True

    expand=False
    if normalize is True or smooth is True or sub is True:
        expand = True
    
    print "Generating scaled bedgraphs..."
    GT.generate_scaled_bedgraphs2(directory, untagged, organism=organism, start_only=start_only, stranded=stranded, threads=threads, file_provided=file_provided, expand=expand)
    
    base_dir = directory.split('/')[:-1]
    base_dir = '/'.join(base_dir)+'/'

    if normalize is True:
        print "\nNormalizing to untagged..."
        if untagged.endswith('.bam'):
            untagged = untagged.split('/')[-1].split('.bam')[0]

        bg_list = [base_dir+x for x in os.listdir(base_dir) if x.endswith('.bedgraph')]
        untagged_bg = [x for x in bg_list if untagged in x][0]
        bg_list.remove(untagged_bg)
        
        last = False
        for n, bg in enumerate(bg_list):
            if n == len(bg_list)-1:
                last = True
            GT.normalize_bedgraph(bg, untagged_bg, smooth=smooth, last=last)

    if smooth is True:
        print "\nSmoothing with {0} bp window...".format(str(window))
        bg_list = [base_dir+x for x in os.listdir(base_dir) if x.endswith('.bedgraph')]
        #if file_provided:
        #    name = directory.split('/')[-1].split('.bam')[0]
        #    bg_list = [x for x in bg_list if name in bg_list]
        GT.smooth_bedgraphs(bg_list, window)
        
    if sub is True:
        print "\nSubtracting background..."
        bg_list = [base_dir+x for x in os.listdir(base_dir) if x.endswith('.bedgraph')]
        p = Pool(threads/2)
        p.map(GT.background_subtraction, bg_list)
        
if __name__ == "__main__":
    main()
