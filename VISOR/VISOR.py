#!/usr/bin/python3 env

import os
import sys
import argparse
from argparse import HelpFormatter, RawTextHelpFormatter
from datetime import datetime

from VISOR import __version__

#v1.0->v1.1 added logo

def barcodefile(path):
	return os.path.abspath(path)

def main():

	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This program was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI). Extensive documentation is available at: https://davidebolo1993.github.io/visordoc/''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='HACk,SHORtS,LASeR,XENIA')

	## HACk ##

	#v1.0->v1.1 general code improvements
	#v1.0->v1.1 -bed changed to -b (-bed is not accepted anymore)
	#v1.0->v1.1 removed logging module dependency. Now printing warnings/errors/messages. Easier debugging
	#v1.0->v1.1 final FASTA is 60-char per line FASTA. Indexed using pyfaidx. Post-processing by samtools is not required
	#v1.1->v1.2_fix Fix index error when parsing initial .bed file

	parser_hack = subparsers.add_parser('HACk', help='HAplotype Creator. Generates one or more haplotypes in FASTA format containing SVs specified in BED-like format')

	required = parser_hack.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome in FASTA format', metavar='FASTA', required=True)
	required.add_argument('-b', '--bedfile', help='one or more variant file in BED-like format (one for each haplotype), as specified in https://davidebolo1993.github.io/visordoc/usage/usage.html#hack-bed', metavar='BED', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', help='output folder with FASTA haplotypes', metavar='FOLDER', required=True)

	#this won't be used yet. This is just for future reference and implementation (Store a VCF with variants inserted)


	optional = parser_hack.add_argument_group('Additional parameters')

	optional.add_argument('--vcf', help=argparse.SUPPRESS, action='store_true')

	parser_hack.set_defaults(func=run_subtool)

	## SHORtS ##
	
	#v1.0->v1.1 general code improvements
	#v1.0->v1.1 -bed changed to -b (-bed is not accepted anymore)
	#v1.0->v1.1 removed logging module dependency. Now printing warnings/errors/messages. Easier debugging
	#v1.0->v1.1 mapping performed using minimap2 (short-read preset with paired-end reads) and directly piped into samtools sort (sorted BAM files are required for final merging)
	#v1.0->v1.1 not tagging by default
	#v1.0->v1.1 Removed several calls to samtools index which were not required. A unique index is built in the end, after merging
	#v1.0->v1.1 Short-read simulations are performed using pywgsim (no need to call wgsim from bash)
	#v1.0->v1.1 SHORtS now integrates a previous script that performed simulations of multiple cells by calling SHORtS in strand-seq mode multiple times. This is now handled internally by SHORtS.
	#v1.0->v1.1 Removed calls to external bash scripts

	parser_shorts = subparsers.add_parser('SHORtS', help='SHOrt Reads Simulator. Simulate piared-end short-read BAM from FASTA haplotypes using regions specified in BED-like format')

	required = parser_shorts.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome in FASTA format', metavar='FASTA', required=True)
	required.add_argument('-s','--sample', help='one or more folders containing FASTA haplotypes with variants implanted by HACk. If multiple folders are given, each sample is considered a subclone and a percentage must be specified accordingly', metavar='FOLDER',  nargs='+', action='append', required=True)
	required.add_argument('-b','--bedfile', help='one or more regions to simulate in BED-like format, as specified in https://davidebolo1993.github.io/visordoc/usage/usage.html#visor-shorts-and-laser', metavar='BED', required=True)
	required.add_argument('-o','--output', help='output folder with simulated clone/haplotype-tagged BAM', metavar='FOLDER', required=True)

	wgsim = parser_shorts.add_argument_group('Paired-end FASTQ simulation using pywgsim (https://github.com/ialbert/pywgsim)') #wgsim core inside

	wgsim.add_argument('--coverage', help='mean coverage for simulated regions [30.0]', metavar='', default=30.0, type=float) #required number of read calculated internally
	wgsim.add_argument('--error', help='base error rate [0.02]', metavar='', default=0.02, type=float)
	wgsim.add_argument('--distance', help='outer distance between the two ends [500]', metavar='', default=500, type=int)
	wgsim.add_argument('--stdev', help='standard deviation of the distance between the two ends [50]',metavar='', default=50, type=int)
	wgsim.add_argument('--length', help='length of reads [150]', metavar='', default=150, type=int)
	wgsim.add_argument('--mutation', help='mutation rate [0.001]', metavar='', default=0.001, type=float)
	wgsim.add_argument('--indels', help='indels rate [0.15]', metavar='', default=0.15, type=float)
	wgsim.add_argument('--extindels', help='indels extension rate [0.25]', metavar='', default=0.25, type=float)

	bulk = parser_shorts.add_argument_group('Simulation of bulk data [default behaviour]')

	bulk.add_argument('--clonefraction', help='ordered percentages (floats) for each input clone (if multiple) [None]', metavar='', nargs='+', action='append', default=None)

	strandseq = parser_shorts.add_argument_group('Simulation of strand-seq data. Coverage must be adjusted coherently')

	strandseq.add_argument('--strandseq', help='perform simulation of strand-seq data', action='store_true')
	strandseq.add_argument('--cells', help='number of cells to simulate [1]', metavar='', type=int, default=1)
	strandseq.add_argument('--refcells', help='percentage of reference cells to simulate. The others are simulated from the input sample. Reference is assumed diploid [0.0]', metavar='', type=float, default=0)
	strandseq.add_argument('--sce', help='one or more regions containing sister chromatid exchange (SCE) events to simulate in BED-like format, as aspecified in https://davidebolo1993.github.io/visordoc/usecases/usecases.html#visor-shorts-and-laser', metavar='', default=None)
	strandseq.add_argument('--noise', help='percentage of noise (reads with incorrect orientation) to add to strand-seq BAM files [0.0]', type=float, metavar='', default=0.0)

	optional = parser_shorts.add_argument_group('Additional parameters')
		
	optional.add_argument('--threads', help='number of cores to use for mapping (minimap2 with preset for short reds) [1]', metavar='', type=int, default=1)
	optional.add_argument('--tag', help='tag simulated BAM by clone (CL-tag) and haplotype (HP-tag). Does not apply to strand-seq data, where cells and haplotypes are separated', action='store_true')
	optional.add_argument('--fastq', help='store synthetic read pairs in FASTQ format in the output folder. Does not work for strand-seq data', action='store_true')
	optional.add_argument('--compress', help='gzip compress output FASTQ', action='store_true')

	parser_shorts.set_defaults(func=run_subtool)

	## LASeR ##

	#v1.0->v1.1 general code improvements
	#v1.0->v1.1 -bed changed to -b (-bed is not accepted anymore)
	#v1.0->v1.1 removed logging module dependency. Now printing warnings/errors/messages. Easier debugging
	#v1.0->v1.1 mapping performed using minimap2 (long-read preset map-pb/map-ont) and directly piped into samtools sort (sorted BAM files are required for final merging)
	#v1.0->v1.1 not tagging by default
	#v1.0->v1.1 Removed several calls to samtools index which were not required. A unique index is built in the end, after merging
	#v1.0->v1.1 Long-read simulations are performed using BadRead. This is slower than Pbsim, but can be accelerated using multiple cores. BadRead gives users control over a wide variety of parameters, which were not supported by Pbsim
	#v1.0->v1.1 Removed calls to external bash scripts
	#v1.1->v1.2_fix Change badread to badread v4.0, allow floats, change models

	parser_long = subparsers.add_parser('LASeR', help='Long reAds SimulatoR. Simulate long-read BAM from FASTA files using regions specified in BED-like format')

	required = parser_long.add_argument_group('Required I/O arguments')
	
	required.add_argument('-g','--genome', help='reference genome in FASTA format', metavar='FASTA', required=True)
	required.add_argument('-s','--sample', help='one or more folders containing FASTA haplotypes with variants implanted by HACk. If multiple folders are given, each sample is considered a subclone and a percentage must be specified accordingly', metavar='FOLDER',  nargs='+', action='append', required=True)
	required.add_argument('-b','--bedfile', help='one or more regions to simulate in BED-like format, as specified in https://davidebolo1993.github.io/visordoc/usage/usage.html#visor-shorts-and-laser', metavar='BED', required=True)
	required.add_argument('-o','--output', help='output folder with simulated clone/haplotype-tagged BAM', metavar='FOLDER', required=True)

	badread= parser_long.add_argument_group('FASTQ simulation using Badread (https://github.com/rrwick/Badread)')

	badread.add_argument('--coverage', help='mean coverage for simulated regions [30.0]', metavar='', default=30.0, type=float) #badread handles this internally
	badread.add_argument('--length_mean', help='mean length of simulated reads [15000]', metavar='', default=15000, type=int)
	badread.add_argument('--length_stdev', help='length stdev of simulated reads [13000]', metavar='', default=13000, type=int)
	badread.add_argument('--identity_min', help='minimum sequencing identity [90.0]', metavar='', default=90.0, type=float)
	badread.add_argument('--identity_max', help='maximum sequencing identity [95.0]', metavar='', default=95.0, type=float)
	badread.add_argument('--identity_stdev', help='stdev of sequencing identity [5.0]', metavar='', default=5.0, type=float)
	badread.add_argument('--error_model', help='error model. Can be "nanopore2023", "nanopore2020", "nanopore2018", "pacbio2016", "random" or a model provided by the user using instructions at https://github.com/rrwick/Badread/wiki/Generating-error-and-qscore-models [nanopore2023]', metavar='', default='nanopore2023', type=str)
	badread.add_argument('--qscore_model', help='quality score model. Can be "nanopore2023", "nanopore2020", "nanopore2018", "pacbio2016", "random", "ideal" or a model provided by the user using instructions at https://github.com/rrwick/Badread/wiki/Generating-error-and-qscore-models [nanopore2023]', metavar='', default='nanopore2020', type=str)
	badread.add_argument('--junk_reads', help='percentage of low-complexity reads [0.5]', metavar='', default=0.5, type=float)
	badread.add_argument('--random_reads', help='percentage of randomly generated reads [0.0]', metavar='', default=0.0, type=float)
	badread.add_argument('--chimera_reads', help='percentage of different reads that are joined together [0.0]', metavar='', default=0.0, type=float)
	badread.add_argument('--glitches_rate', help='number of bases between portions of a read where the sequence is messed up (glitch) [5000]', metavar='', default=5000, type=int)
	badread.add_argument('--glitches_size', help='length of glitches [50]', metavar='', default=10, type=int)
	badread.add_argument('--glitches_skip', help='number of bases lost in glitches [50]', metavar='', default=10, type=int)

	bulk = parser_long.add_argument_group('Simulation of bulk data')

	bulk.add_argument('--clonefraction', help='ordered percentages (floats) for each input clone (if multiple) [None]', metavar='', nargs='+', action='append', default=None)

	optional = parser_long.add_argument_group('Additional parameters')
		
	optional.add_argument('--read_type', help='type of reads, used for mapping. Can be either "nanopore" or "pacbio" [nanopore]', default='nanopore', choices=['nanopore', 'pacbio'], metavar='')
	optional.add_argument('--threads', help='number of cores to use for simulation (BadRead is rather slow with a single core) and mapping (minimap2 with preset for nanopore or pacbio reads) [1]', metavar='', type=int, default=1)
	optional.add_argument('--tag', help='tag simulated BAM by clone (CL-tag) and haplotype (HP-tag)', action='store_true')
	optional.add_argument('--fastq', help='store synthetic reads in FASTQ format in the output folder', action='store_true')
	optional.add_argument('--compress', help='Compress output FASTQ using gzip', action='store_true')

	parser_long.set_defaults(func=run_subtool)

	## XENIA ## [BETA. Tested and it works, but it needs to be used by others to be sure this fits]

	#v1.0->v1.1 general code improvements
	#v1.0->v1.1 -bed changed to -b (-bed is not accepted anymore)
	#v1.0->v1.1 removed logging module dependency. Now printing warnings/errors/messages. Easier debugging
	#v1.0->v1.1 Short-read simulations are performed using pywgsim (no need to call wgsim from bash)
	#v1.0->v1.1 Removed calls to external bash scripts
	xeniatext = """
Simulate linked-read FASTQ using regions specified in standard BED format. XENIA can optionally take a file of nucleotide barcodes, 1 per line. \
You can specify the linked-read barcode type to use (--barcode-type) as well as the output format of FASTQ files (default is the same as barcode type). For example, you \
can provide a file of 96 barcodes (commonly haplotagging style), select --barcode-type stlfr (combinatorial 3-barcode on R2 read), and have \
--output-format tellseq (@seqid:barcode header format).

Barcode Types
10x/tellseq   single barcode on R1
stlfr         combinatorial 3-barcode on R2
haplotagging  R1 and R2 each have different combinatorial 2-barcodes

Output Types
10x           barcode in sequence at start of R1
haplotagging  barcode appended to sequence header via 'BX:Z' tag
              encoded in ACBD combinatorial format
standard      barcode appended to sequence header via `BX:Z` tag
              but no specific barcode format (uses input style)
stlfr         barcode appended to sequence ID via '#'
              encoded in 1_2_3 combinatorial format
tellseq       barcode appended to sequence ID via ':'
	"""
	parser_linkedreads = subparsers.add_parser('XENIA', help = '10X gENomics sImulAtor. Simulate many kinds of linked-read data types.', description = xeniatext, formatter_class=RawTextHelpFormatter)

	required = parser_linkedreads.add_argument_group('I/O arguments')

	required.add_argument('bedfile', help='one or more regions to simulate in BED format', metavar='BEDFILE')
	required.add_argument('fasta', help='folder containing FASTA haplotypes', metavar='FASTA_FOLDER')
	required.add_argument('-o','--output', help='output directory (required)', metavar='FOLDER', required=True)
	required.add_argument('-p','--prefix', help='prefix for FASTQ files [SIM]', metavar='PREFIX')
	required.add_argument('-O','--output-format', help='output format for FASTQ files [--barcode-type]\n{10x, haplotagging, standard, stlfr, tellseq}', metavar='CHOICE', choices = ["10x", "stlfr", "standard", "haplotagging", "tellseq"])

	wgsim = parser_linkedreads.add_argument_group('Paired-end FASTQ simulation using pywgsim')

	wgsim.add_argument('--coverage', help='mean coverage for simulated regions [30.0]', metavar='', default=30.0, type=float) #required number of read calculated internally
	wgsim.add_argument('--error', help='base error rate [0.02]', metavar='', default=0.02, type=float)
	wgsim.add_argument('--distance', help='outer distance between the two ends [500]', metavar='', default=500, type=int)
	wgsim.add_argument('--stdev', help='standard deviation of --distance [50]',metavar='', default=50, type=int)
	wgsim.add_argument('--length', help='length of reads [150]', metavar='', default=150, type=int)
	wgsim.add_argument('--mutation', help='mutation rate [0.001]', metavar='', default=0.001, type=float)
	wgsim.add_argument('--indels', help='indels rate [0.15]', metavar='', default=0.15, type=float)
	wgsim.add_argument('--extindels', help='indels extension rate [0.25]', metavar='', default=0.25, type=float)

	molecules=parser_linkedreads.add_argument_group('Linked-read simulation')
	bc_default = os.path.abspath(os.path.dirname(__file__) + '/XENIA/4M-with-alts-february-2016.txt.gz')
	molecules.add_argument('-b', '--barcodes', help='file of barcodes to simulate, one per line [10X barcodes]', metavar='FILE', default = bc_default, type=barcodefile)
	molecules.add_argument('-B', '--barcode-type', help='Type of linked-read experiment [10x]\n{10x, haplotagging, stlfr, tellseq}', metavar='CHOICE', default = "10x", choices = ["10x", "stlfr", "haplotagging", "tellseq"])
	molecules.add_argument('--molecule-length', help='mean length of molecules [80000]', metavar='', default=80000, type=int)
	molecules.add_argument('--molecule-number', help='number of fragments per barcode on average [10]', metavar='', default=10, type=int)
	molecules.add_argument('--molecule-coverage', help='if <1, mean coverage per molecule [0.2]\nif >=1, mean number of reads per molecule', metavar='', default=0.2, type=float)

	optional = parser_linkedreads.add_argument_group('Additional parameters')

	optional.add_argument('--threads', help='number of cores to use for simulation [1]', metavar='', type=int, default=1)

	parser_linkedreads.set_defaults(func=run_subtool)


	#print help if no subcommand nor --help provided

	print(fr"""

	 ___      ___ ___  ________  ________  ________     
	|\  \    /  /|\  \|\   ____\|\   __  \|\   __  \    
	\ \  \  /  / | \  \ \  \___|\ \  \|\  \ \  \|\  \   
	 \ \  \/  / / \ \  \ \_____  \ \  \\\  \ \   _  _\  
	  \ \    / /   \ \  \|____|\  \ \  \\\  \ \  \\  \| 
	   \ \__/ /     \ \__\____\_\  \ \_______\ \__\\ _\ 
	    \|__|/       \|__|\_________\|_______|\|__|\|__| v{__version__}
	                     \|_________|                   
                                                    	                                                                        
	""", file = sys.stderr)
	
	
	if len(sys.argv)==1:
    	
		parser.print_help(sys.stderr)
		sys.exit(1)

	#case-insensitive submodules
	
	if sys.argv[1].lower() == 'hack':

		sys.argv[1] = 'HACk'

	elif sys.argv[1].lower() == 'shorts':

		sys.argv[1] = 'SHORtS'

	elif sys.argv[1].lower() == 'laser':

		sys.argv[1] = 'LASeR'

	elif sys.argv[1].lower() == 'xenia':

		sys.argv[1] = 'XENIA'

	args = parser.parse_args()

	now = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print(f'[{now}] Running VISOR with the following arguments:\n', file = sys.stderr)

	maxl = max(len(arg) for arg in vars(args))
	for arg in vars(args):
		if arg not in ['func']:
			arg_value = getattr(args, arg)
			print('{:<{width}}: {}'.format(arg, arg_value, width=maxl), file = sys.stderr)

	print("", file = sys.stderr)
	args.func(parser, args)



class CustomFormat(HelpFormatter):

	'''
	Customize how help is diplayed
	'''

	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()



def run_subtool(parser, args):


	if args.command == 'HACk': 

		from .HACk import HACk as submodule
	
	elif args.command == 'SHORtS': 

		from .SHORtS import SHORtS as submodule

	elif args.command == 'LASeR':

		from .LASeR import LASeR as submodule

	elif args.command == 'XENIA':

		from .XENIA import XENIA as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':

	main()

