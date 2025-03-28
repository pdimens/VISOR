#!/usr/bin/python3 env

#python 3 standard library

import os
import sys
import glob
import re
import math
#import resource
import gzip
import multiprocessing
from datetime import datetime
from shutil import which

#additional modules

import pybedtools
import pyfaidx
from pywgsim import wgsim
import numpy as np

from VISOR import __version__

#barcodepath=os.path.abspath(os.path.dirname(__file__) + '/4M-with-alts-february-2016.txt.gz')

class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	OUT = ''
	BED = ''
	SAMPLE = ''

	#pywgsim

	coverage=0
	regioncoverage=0
	error=0
	distance=0
	stdev=0
	length=0
	mutation=0
	indels=0
	extindels=0

	#bulk

	ffiles=None
	fperc=0.0
	ffile=None
	hapnumber=0
	threads=0

	#molecules
	barcodepath=""
	molnum=0
	mollen=0
	molcov=0
	barcodebp = 0
	barcodes=set()


class Molecule(object):

	'''
	Molecule instance
	'''
		
	def __init__(self,length,start,end,index):
				
		self.seqidx=index
		self.length=length
		self.index_droplet=0
		self.barcode=None
		self.start=start
		self.end=end

	def __str__(self):
		
		return str(self.__class__) + ": " + str(self.__dict__)


def redirect_stdout():

	'''
	Suppress c/c++/fortran stdout, keep python print calls
	'''

	sys.stdout.flush() # <--- important when redirecting to files
	newstdout = os.dup(1)
	devnull = os.open(os.devnull, os.O_WRONLY)
	os.dup2(devnull, 1)
	os.close(devnull)
	sys.stdout = os.fdopen(newstdout, 'w')


def atoi(text):

	'''
	Convert text to integers
	'''

	return int(text) if text.isdigit() else text


def natural_keys(text):

	'''
	Natural sort
	'''
	
	return [ atoi(c) for c in re.split(r'(\d+)', text)]


def Chunks(l,n):

	'''
	Split list in chunks based on number of threads
	'''

	return [l[i:i+n] for i in range(0, len(l), n)]


def readfq(fp): # this is a fast generator function

	'''
	Yield FASTQ record
	'''

	read=[]

	for line in fp:

		read.append(line.rstrip())

		if len(read) == 4:

			yield read[0][1:], read[1], read[3]

			read=[]


def randomlong(Par,seq_,EXPM):

	'''
	Length of molecules is randomly distributed according to an exponential distribution
	'''

	index=0	
	lensingle=len(seq_)
	molecules=[]
				
	for i in range(EXPM):
						
		start=int(np.random.uniform(low=0,high=lensingle))
		length=int(np.random.exponential(scale=c.mollen))

		if length!=0:
					
			end=start+length-1

			if end>lensingle:
							 
				Molseq=seq_[start:lensingle]
				lengthnew=lensingle-start
				NewMol=Molecule(lengthnew,start,lensingle,index)
				molecules.append(NewMol)

			else:
							 
				Molseq=seq_[start:end]
				NewMol=Molecule(length-1,start,end,index)
				molecules.append(NewMol)
		
			index+=1

	return molecules


def deternumdroplet(molecules,molnum):

	'''
	Determine the number of droplets
	'''

	large_droplet=4000000
	assign_drop=[]

	frag_drop = np.random.poisson(molnum,large_droplet)
	totalfrag=0
	
	for i in range(large_droplet):
				
		totalfrag=totalfrag+frag_drop[i]
		
		if totalfrag<=len(molecules):
		
			assign_drop.append(frag_drop[i])
		
		else:
			
			last=len(molecules)-(totalfrag-frag_drop[i])
			assign_drop.append(last)
			break

	return assign_drop


def selectbarcode(drop,molecules,c):

	'''
	Select barcode
	'''

	permutnum=np.random.permutation(len(molecules))
	N_droplet=len(drop)
	assigned_barcodes=set()
	droplet_container=[]

	start=0
	
	for i in range(N_droplet):
				
		num_molecule_per_partition=drop[i]
		index_molecule=permutnum[start:start+num_molecule_per_partition]
		#print(index_molecule)
		totalseqlen=0
		temp=[]
		start=start+num_molecule_per_partition
		bc = c.barcodes.pop()
		for j in range(num_molecule_per_partition):
						
			index=index_molecule[j]
			temp.append(index)
			molecules[index].index_droplet=i
			molecules[index].barcode=bc
			totalseqlen=totalseqlen+molecules[index].length
		
		assigned_barcodes.add(bc)
		droplet_container.append(temp)

	return droplet_container, assigned_barcodes


def Gzipper(sli,):

	'''
	Compressing here is quite slow. Multi-processing saves some time
	'''

	for s in sli:
		with open(s, 'rb') as textin, gzip.open(s+'.gz', 'wb', compresslevel=6) as gzout:
			gzout.writelines(textin)
		os.remove(s)


def MolSim(processor,molecule,hfa,w,c):

	'''
	Parallelize linked reads simulation
	'''

	for mol in molecule:

		moleculenumber=str(mol.seqidx+1)
		moleculedroplet=str(mol.index_droplet+1)
		barcodestring=str(mol.barcode)
		chromstart=str(w.start+mol.start)
		chromend=str(w.start+mol.end)

		header='MOL:' + moleculenumber + '_GEM:' + moleculedroplet + '_BAR:' + barcodestring + '_CHROM:' + w.chrom + '_START:' + chromstart + '_END:' + chromend
		seq__=hfa[w.chrom][w.start+mol.start-1:w.start+mol.end].seq

		truedim=mol.length-seq__.count('N')
		N=int(truedim*c.molcov)/(c.length*2)

		R1A=os.path.abspath(c.OUT + '/SIM_S1_L' + str(c.hapnumber).zfill(3) + '_R1_001.fastq')
		R2A=os.path.abspath(c.OUT + '/SIM_S1_L' + str(c.hapnumber).zfill(3) + '_R2_001.fastq')
		
		if N != 0:
			molfa=os.path.abspath(c.OUT+'/' + processor + '_' + moleculenumber + '.fa')

			with open(molfa, 'w') as faout:
				faout.write('>' + header + '\n' + '\n'.join(re.findall('.{1,60}', seq__)) + '\n')

			R1tmp=os.path.abspath(c.OUT+'/' + processor + '.R1.tmp.fq')
			R2=os.path.abspath(c.OUT+'/' + processor + '.R2.fq')

			wgsim.core(
				r1 = R1tmp,
				r2 =R2,
				ref = molfa,
				err_rate = c.error,
				mut_rate = c.mutation,
				indel_frac = c.indels,
				indel_ext = c.extindels,
				N = N,
				dist = c.distance,
				stdev = c.stdev,
				size_l = c.length - (c.barcodebp), # (c.barcodebp + 6) if using random 6mer
				size_r = c.length,
				max_n = 0.05,
				is_hap = 0,
				is_fixed = 0,
				seed = 0
			)

			os.remove(molfa)
			#RANDOM6MER=''.join(np.random.choice(['A','T','G','C','N'],6, replace=True))

			if os.stat(R1tmp).st_size == 0:
				os.remove(R1tmp)
				os.remove(R2)

			else:

				with open(R1tmp,'r') as infile, open(R1A,'a') as outfile:
					for name,seq,qual in readfq(infile):
						#read=[f'@{name}', f"{barcodestring}{RANDOM6MER}{seq}", '+', str(qual[0])*(c.barcodebp+6)+qual]
						read=[f'@{name}', f"{barcodestring}{seq}", '+', f'{qual[0] * c.barcodebp}{qual}']
						outfile.write('\n'.join(read) + '\n')
				os.remove(R1tmp)

				with open(R2,'r') as infile, open(R2A,'a') as outfile:
					for name,seq,qual in readfq(infile):
						read=[f'@{name}',seq,'+',qual]
						outfile.write('\n'.join(read) + '\n')
				os.remove(R2)


def LinkedSim(w,c):

	'''
	Perform linked-reads simulation

	'''

	hfa=pyfaidx.Fasta(c.ffile)

	if w.chrom not in hfa.keys():
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Warning] Chromosome {w.chrom} not found in {c.ffile}. Skipped simulation', file = sys.stderr)
	else:
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Message] Preparing simulation from {c.ffile}. Haplotype {c.hapnumber}', file = sys.stderr)

		chr_= hfa[w.chrom]
		seq_ = chr_[w.start-1:w.end].seq
		#tmpfa=os.path.abspath(c.OUT + '/' + 'htmp.fa')
		region=w.chrom+'_'+str(w.start)+'_'+str(w.end)

		#with open(tmpfa, 'w') as tmpfout: #write temporary fa for sampling reads

			#tmpfout.write('>' + region + '\n' + '\n'.join(re.findall('.{1,60}', seq_)) + '\n')

		Ns=seq_.count('N') #normalize coverage on Ns

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Message] Number of available barcodes: {len(c.barcodes)}', file = sys.stderr)

		MRPM=(c.molcov*c.mollen)/(c.length*2)
		TOTALR=round(((c.regioncoverage*(len(seq_)-Ns))/c.length)/2)
		EXPM=round(TOTALR/MRPM)

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Message] Average number of paired reads per molecule: {MRPM}', file = sys.stderr)
		print(f'[{now}][Message] Number of reads required to get the expected coverage: {TOTALR}', file = sys.stderr)
		print(f'[{now}][Message] Expected number of molecules: {EXPM}', file = sys.stderr)

		molecules=randomlong(c,seq_,EXPM) #MolSet

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Message] Molecules generated: {len(molecules)}', file = sys.stderr)

		drop=deternumdroplet(molecules,c.molnum)

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Message] Assigned molecules to: {len(drop)} partitions', file = sys.stderr)

		droplet_container,assigned_barcodes=selectbarcode(drop,molecules,c)

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Message] Assigned a unique barcode to each molecule', file = sys.stderr)
		# remove the barcodes that were used
		#[c.barcodes.remove(x) for x in assigned_barcodes]
		#c.barcodes=[x for x in c.barcodes if x not in assigned_barcodes]

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Message] {len(c.barcodes)} barcodes left', file = sys.stderr)
		print(f'[{now}][Message] Simulating', file = sys.stderr)

		chunk_size=len(molecules)/c.threads
		slices=Chunks(molecules,math.ceil(chunk_size))

		processes=[]

		for i,molecule in enumerate(slices):

			processor=f'p{i+1}'
			p=multiprocessing.Process(target=MolSim, args=(processor,molecule,hfa,w,c))
			p.start()
			processes.append(p)
		
		for p in processes:
		
			p.join()


def run(parser,args):

	'''
	Check arguments, run functions
	'''

	redirect_stdout()# block pywgsim stdout

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print(f'[{now}][BETA] VISOR XENIA v{__version__}', file = sys.stderr)

	#fill container

	c.OUT=os.path.abspath(args.output)
	c.BED=os.path.abspath(args.bedfile)
	c.SAMPLE=os.path.abspath(args.sample)
	c.threads=args.threads

	if not os.path.exists(c.OUT):

		try:

			os.makedirs(c.OUT)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print(f'[{now}][Error] Cannot create the output folder', file = sys.stderr)
			sys.exit(1)

	else:

		if not os.access(os.path.abspath(c.OUT),os.W_OK):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print(f'[{now}][Error] Missing write permissions on the output folder', file = sys.stderr)
			sys.exit(1)
			
		elif os.listdir(os.path.abspath(c.OUT)):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print(f'[{now}][Error] The output folder is not empty: specify another output folder or clean the current one ({c.OUT})', file = sys.stderr)
			sys.exit(1)


	if which('bedtools') is None:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Error] bedtools must be in PATH', file = sys.stderr)
		sys.exit(1)

	try:

		bedfile=pybedtools.BedTool(c.BED)
		bedsrtd=bedfile.sort()

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Error] BED {c.BED} does not exist, is not readable, or is not a valid BED', file = sys.stderr)
		sys.exit(1)

	#fill c with wgsim parameters for simulations 

	c.coverage=args.coverage
	c.error=args.error
	c.distance=args.distance
	c.stdev=args.stdev
	c.length=args.length
	c.mutation=args.mutation
	c.indels=args.indels
	c.extindels=args.extindels
	c.molnum=args.molecule_number
	c.mollen=args.molecule_length
	c.molcov=args.molecule_coverage
	c.barcodepath = args.barcodes
	#if args.barcodes:
	#	c.barcodepath = os.path.abspath(args.barcodes)
	fasta_files = [
    	f for f in glob.glob(f'{os.path.abspath(c.SAMPLE)}/*') 
    	if re.search(r'\.(fa|fasta)$', f, re.IGNORECASE)
	]
	if len(fasta_files) == 0:
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Error] No FASTA files detected in {c.SAMPLE}. If your FASTA files are gzipped, please decompress them.', file = sys.stderr)
		sys.exit(1)
	c.ffiles=sorted(fasta_files, key=natural_keys) #list all FASTA in folder
	c.regioncoverage=c.coverage/len(c.ffiles)

	try:

		with gzip.open(c.barcodepath, 'rt') as filein:
			c.barcodes = filein.read().splitlines()
	except gzip.BadGzipFile:
		with open(c.barcodepath, 'r') as filein:
			c.barcodes = filein.read().splitlines()
	except:
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Error] Cannot open {c.barcodepath} for reading', file = sys.stderr)
		sys.exit(1)
	# validate barcodes are all the same length
	bc_lens = set()
	for i in c.barcodes:
		bc_lens.add(len(i))
		if len(bc_lens) > 1:
			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print(f'[{now}][Error] Barcodes in {c.barcodepath} must all the same length.', file = sys.stderr)
			sys.exit(1)
	else:
		c.barcodebp = next(iter(bc_lens))
	# validate barcodes are only ATCGU nucleotides
	for bc in c.barcodes:
		if not bool(re.fullmatch(r'^[ATCGU]+$', bc, flags = re.IGNORECASE)):
			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print(f'[{now}][Error] Barcodes can only contain nucleotides A,T,C,G,U, but invalid barcode(s) provided: {bc}. This was first invalid barcode identified, but it may not be the only one.', file = sys.stderr)
			sys.exit(1)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print(f'[{now}][Message] Preparing for bulk simulations with a single clone', file = sys.stderr) #maybe we want to add more here in the future. 

	for k,s in enumerate(c.ffiles):
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print(f'[{now}][Message] Processing haplotype {k+1}', file = sys.stderr)
		c.hapnumber=str(k+1)
		c.ffile=c.ffiles[k]

		for w in bedsrtd: #do not use multi-processing on this as minimap2 may require too much memory

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print(f'[{now}][Message] Simulating from region {w.chrom}:{w.start}-{w.end}', file = sys.stderr)
			LinkedSim(w,c)

	allfastq=glob.glob(os.path.abspath(c.OUT) + '/*.fastq')

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print(f'[{now}] Compressing FASTQ', file = sys.stderr)

	#gzip multiprocessing

	chunk_size=len(allfastq)/c.threads
	slices=Chunks(allfastq,math.ceil(chunk_size))
	processes=[]

	for i,sli in enumerate(slices):

		p=multiprocessing.Process(target=Gzipper, args=(sli,))
		p.start()
		processes.append(p)
		
	for p in processes:
		
		p.join()

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print(f'[{now}] Done', file = sys.stderr)
