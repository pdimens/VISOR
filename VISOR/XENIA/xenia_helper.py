#! /usr/bin/python3 env

import os
import sys
import glob
import re
import math
import gzip
import pysam
import multiprocessing
from datetime import datetime
from itertools import product
from shutil import which

#additional modules

import pybedtools
import pyfaidx
from pywgsim import wgsim
import numpy as np


class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	OUT = ''
	BED = ''
	FASTADIR = ''

	#pywgsim

	coverage=0
	regioncoverage=0
	error=0
	distance=0
	stdev=0
	length=0
	length_r1=0
	length_r2=0
	mutation=0
	indels=0
	extindels=0

	#bulk

	ffiles=None
	fperc=0.0
	ffile=None
	outformat=None
	hapnumber=0
	threads=0

	#molecules
	barcodepath=None
	molnum=0
	mollen=0
	molcov=0
	barcodetype=None
	barcodebp=0
	barcodes=None	# will be an iterable to use with next()
	bc_generator = None
	used_bc={}
	totalbarcodes=0
	remainingbarcodes=0
	whitelist=None

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

def get_now():
	return datetime.now().strftime('%d/%m/%Y %H:%M:%S')

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

def validate_barcodes(bc_list):
	'''
	Takes a list of barcodes read from input barcode file and validates them according to the linked-read type.
	Validations include: ATGCU nucleotides, barcodes same length (or same length within column).
	'''
	# check first row for multiple columns, if there are multiple, it's haplotagging
	if len(bc_list[0].strip().split()) != 1:
		print(f'[{get_now()}][Error] Barcode file is expected to only have one barcode per line)', file = sys.stderr)
		sys.exit(1)
	else:
		bc_lens = set()
		for i in bc_list:
			bc_lens.add(len(i))
			if len(bc_lens) > 1:
				print(f'[{get_now()}][Error] Barcodes provided must all be the same length. If these are haplotagging barcodes, then barcodes within a column must be the same length (but the length can vary across rows)', file = sys.stderr)
				sys.exit(1)
		# validate barcodes are only ATCGU nucleotides
		for bc in bc_list:
			if not bool(re.fullmatch(r'^[ATCGU]+$', bc, flags = re.IGNORECASE)):
				print(f'[{get_now()}][Error] Barcodes can only contain nucleotides A,T,C,G,U, but invalid barcode(s) provided: {bc}. This was first invalid barcode identified, but it may not be the only one.', file = sys.stderr)
				sys.exit(1)

def interpret_barcodes(infile, lr_type):
	"""
	Takes an open file connection and reads it line by line. Performs barcode validations and returns:
	- either an iter() or generator of barcodes (to use with next())
	- the total barcode	length (int)
	- the total number of barcodes [or combinations] (int)
	"""
	bc = list(set(i.strip() for i in infile.read().splitlines()))
	validate_barcodes(bc)
	bc_len = len(bc[0]) 
	if lr_type == "haplotagging":
		# 2 barcodes per
		return product(bc,bc,bc,bc), 2 * bc_len, len(bc)**4
	if lr_type == "stlfr":
		return product(bc,bc,bc), 3 * bc_len, len(bc)**3
	else:
		return iter(bc), bc_len, len(bc)


def format_linkedread(name, bc, seq, qual):
	'''
	Given a linked-read output type, will format the read accordingly and return it
	'''
	if c.outformat == "10x":
		read = [f'@{name}', f"{bc}{seq}", '+', f'{qual[0] * c.barcodebp}{qual}']
		if bc not in c.used_bc:
			c.used_bc[bc] = True
			c.whitelist.write(f"{bc}\n")
	elif c.outformat == "tellseq":
		read = [f'@{name}:{bc}', seq, '+', qual]
		if bc not in c.used_bc:
			c.whitelist.write(f"{bc}\n")
			c.used_bc[bc] = True
	elif c.outformat == "haplotagging":
		if bc not in c.used_bc:
			acbd = "".join(next(c.bc_generator))
			c.used_bc[bc] = acbd
			c.whitelist.write(f"{bc}\n")
		else:
			acbd = c.used_bc[bc]
		read = [f'@{name}\tOX:Z:{bc}\tBX:Z:{acbd}', seq, '+', qual]
	elif c.outformat == "stlfr":
		if bc not in c.used_bc:
			stlfr_bc = "_".join([str(i) for i in next(c.bc_generator)])
			c.used_bc[bc] = stlfr_bc
			c.whitelist.write(f"{bc}\n")
		else:
			stlfr_bc = c.used_bc[bc]
		read = [f'@{name}#{stlfr_bc}', seq, '+', qual]
	return read

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

	large_droplet=c.totalbarcodes
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
	Select barcode to use for droplet/partition
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
		try:
			bc = next(c.barcodes) if c.barcodetype in ["10x", "tellseq"] else "".join(next(c.barcodes))
		except StopIteration:
			print(f'[{get_now()}][Error] No more barcodes left for simulation. The requested parameters require more barcodes.', file = sys.stderr)
			sys.exit(1)
		c.remainingbarcodes -= 1
		for j in range(num_molecule_per_partition):
						
			index=index_molecule[j]
			temp.append(index)
			molecules[index].index_droplet=i
			molecules[index].barcode=bc
			totalseqlen=totalseqlen+molecules[index].length
		
		assigned_barcodes.add(bc)
		droplet_container.append(temp)

	return droplet_container, assigned_barcodes


def BGzipper(sli,):

	'''
	Use pysam/htslib BGzip and multi-processing to save some time
	'''
	for s in sli:
		pysam.tabix_compress(s, f'{s}.gz', force=True)
		os.remove(s)


def MolSim(processor,molecule,hfa,w,c):

	'''
	Parallelize linked reads simulation
	'''

	for mol in molecule:

		moleculenumber = mol.seqidx + 1
		moleculedroplet = mol.index_droplet + 1
		barcodestring = mol.barcode
		chromstart= w.start + mol.start
		chromend = w.start + mol.end

		header=f'MOL:{moleculenumber}_GEM:{moleculedroplet}_BAR:{barcodestring}_CHROM:{w.chrom}_START:{chromstart}_END:{chromend}'
		seq__=hfa[w.chrom][w.start+mol.start-1:w.start+mol.end].seq

		truedim=mol.length-seq__.count('N')
		N=int(truedim*c.molcov)/(c.length*2)

		R1A=os.path.abspath(c.OUT + '/SIM_S1_L' + str(c.hapnumber).zfill(3) + '_R1_001.fastq')
		R2A=os.path.abspath(c.OUT + '/SIM_S1_L' + str(c.hapnumber).zfill(3) + '_R2_001.fastq')

		if N != 0:
			molfa=os.path.abspath(f'{c.OUT}/{processor}_{moleculenumber}.fa')

			with open(molfa, 'w') as faout:
				faout.write(f'>{header}\n' + '\n'.join(re.findall('.{1,60}', seq__)) + '\n')

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
				size_l = c.len_r1,
				size_r = c.len_r2,
				max_n = 0.05,
				is_hap = 0,
				is_fixed = 0,
				seed = 0
			)

			os.remove(molfa)

			if os.stat(R1tmp).st_size == 0:
				os.remove(R1tmp)
				os.remove(R2)

			else:

				with open(R1tmp,'r') as infile, open(R1A,'a') as outfile:
					for name,seq,qual in readfq(infile):
						read = format_linkedread(name, barcodestring, seq, qual)
						outfile.write('\n'.join(read) + '\n')
				os.remove(R1tmp)

				with open(R2,'r') as infile, open(R2A,'a') as outfile:
					for name,seq,qual in readfq(infile):
						if c.outformat == "10x":
							read = [f'@{name}',seq,'+',qual]
						else:
							read = format_linkedread(name, barcodestring, seq, qual)
						outfile.write('\n'.join(read) + '\n')
				os.remove(R2)


def LinkedSim(w,c):

	'''
	Perform linked-reads simulation
	'''

	hfa=pyfaidx.Fasta(c.ffile)

	if w.chrom not in hfa.keys():
		print(f'[{get_now()}][Warning] Chromosome {w.chrom} not found in {c.ffile}. Skipped simulation', file = sys.stderr)
	else:
		print(f'[{get_now()}] Preparing simulation from {c.ffile}. Haplotype {c.hapnumber}', file = sys.stderr)

		chr_= hfa[w.chrom]
		seq_ = chr_[w.start-1:w.end].seq
		#tmpfa=os.path.abspath(c.OUT + '/' + 'htmp.fa')
		region=w.chrom+'_'+str(w.start)+'_'+str(w.end)

		#with open(tmpfa, 'w') as tmpfout: #write temporary fa for sampling reads

			#tmpfout.write('>' + region + '\n' + '\n'.join(re.findall('.{1,60}', seq_)) + '\n')

		Ns=seq_.count('N') #normalize coverage on Ns
		
		print(f'[{get_now()}] Number of available barcodes: {c.remainingbarcodes}', file = sys.stderr)

		MRPM=(c.molcov*c.mollen)/(c.length*2)
		TOTALR=round(((c.regioncoverage*(len(seq_)-Ns))/c.length)/2)
		EXPM=round(TOTALR/MRPM)

		print(f'[{get_now()}] Average number of paired reads per molecule: {MRPM}', file = sys.stderr)
		print(f'[{get_now()}] Number of reads required to get the expected coverage: {TOTALR}', file = sys.stderr)
		print(f'[{get_now()}] Expected number of molecules: {EXPM}', file = sys.stderr)

		molecules=randomlong(c,seq_,EXPM) #MolSet

		print(f'[{get_now()}] Molecules generated: {len(molecules)}', file = sys.stderr)

		drop=deternumdroplet(molecules,c.molnum)
		
		print(f'[{get_now()}] Assigned molecules to: {len(drop)} partitions', file = sys.stderr)

		droplet_container,assigned_barcodes=selectbarcode(drop,molecules,c)

		print(f'[{get_now()}] Assigned a unique barcode to each molecule', file = sys.stderr)

		print(f'[{get_now()}] {c.remainingbarcodes} barcodes left', file = sys.stderr)

		print(f'[{get_now()}] Simulating', file = sys.stderr)

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