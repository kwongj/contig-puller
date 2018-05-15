#!/usr/bin/env python3
# Script by JK
# Extracts targeted contigs from a multi-FASTA file eg. draft genome assembly

import subprocess
import os
import sys
import shutil
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import UnknownSeq
import csv

# Usage
import argparse
parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='Extracts contigs containing a target sequence (eg. gene) from several multi-FASTA assemblies \n and aligns the contigs at the target sequence',
		usage='\n  %(prog)s --db DBFILE --out OUTFILE --anno ANNOFILE [OPTIONS] FASTA-1 FASTA-2 ... FASTA-N')
parser.add_argument('--db', metavar='DBFILE', required=True, help='target sequence (FASTA) to search for (REQUIRED to create BLAST database)')
parser.add_argument('--id', metavar='ID', default='100', help='percentage identity cutoff (default=100)')
parser.add_argument('--out', metavar='OUTFILE', default='contigs.gbk', help='output file (default=contigs.gbk)')
parser.add_argument('--anno', metavar='ANNOFILE', help='reference proteins.faa file to annotate from (optional | requires Prokka to be installed)')
parser.add_argument('--cpus', metavar='CPUS', default='8', help='number of cpus to use (default=8)')
#parser.add_argument('--log', metavar='LOGFILE', default='tempdb/temp.log', help='Log file')
parser.add_argument('assembly', metavar='FASTA', nargs='+', help='FASTA assemblies to search')
parser.add_argument('--version', action='version', version='v2.0')
args = parser.parse_args()


# Functions
def progexit():
	shutil.rmtree('./tempdb')
	banner()
	print('Done. Output written to {}.'.format(outfile))
	banner()
	sys.exit(0)

def banner():
	print('--------------------------')
	return()

# Check output file
outfile = args.out
tmpout = 'tempdb/tmpout'
if (os.path.isfile(outfile)):
	print('ERROR: {} already exists. Please specify another output file.'.format(outfile))
	sys.exit(1)

# Create local BLAST database
print('Creating BLAST database ...')
if os.path.isdir('./tempdb'):
	print('ERROR: ./tempdb already exists. Please remove/rename this directory first.')
	sys.exit(1)
subprocess.call(['mkdir', 'tempdb'])
subprocess.call(['makeblastdb', '-in', args.db, '-parse_seqids', '-dbtype', 'nucl', '-out', 'tempdb/db'])
banner() 

# Run BLAST and find lengths of contigs and orientation/position of target gene
seqlist = []
seqplus = []
seqminus = []
contigs = []
genestart = []
geneend = []
for f in args.assembly:
	seqlist.append(f)
	blastn_for = NcbiblastnCommandline(query=f, db='tempdb/db', word_size=32, strand='plus', dust='no', perc_identity=args.id, evalue=1E-99, outfmt=6, out='tempdb/bplus.txt')
	blastn_for()
	blastn_rev = NcbiblastnCommandline(query=f, db='tempdb/db', word_size=32, strand='minus', dust='no', perc_identity=args.id, evalue=1E-99, outfmt=6, out='tempdb/bminus.txt')
	blastn_rev()
	if os.stat('tempdb/bplus.txt').st_size > 0:
		print('Searching {} ...'.format(f))
		bplus = list(csv.reader(open('tempdb/bplus.txt', 'r'), delimiter='\t'))
		bplus.insert(0, f)
		seqplus.append(bplus)
		cont = str(bplus[1][0])
		contigs.append(cont)
		start = int(bplus[1][6])
		end = int(bplus[1][7])
		for s in SeqIO.parse(f, 'fasta'):
			if s.id == cont:
				print("Found '{}' (+) ...".format(cont))
				genestart.append(start - 1)
				geneend.append(len(s) - end)
	if os.stat('tempdb/bminus.txt').st_size > 0:
		print('Searching {} ...'.format(f))
		bminus = list(csv.reader(open('tempdb/bminus.txt', 'r'), delimiter='\t'))
		bminus.insert(0, f)
		seqminus.append(bminus)
		cont = str(bminus[1][0])
		contigs.append(cont)
		start = int(bminus[1][6])
		end = int(bminus[1][7])
		for s in SeqIO.parse(f, 'fasta'):
			if s.id == cont:
				print("Found '{}' (-) ...".format(cont))
				genestart.append(len(s) - end)
				geneend.append(start - 1)
maxstart = max(genestart)
maxend = max(geneend)
banner()


# Write target contigs to FASTA file
print('Writing forward sequence matches ... ')
for row in seqplus:
	seqname = row[0]
	cont = row[1][0]
	begin = int(row[1][6])
	finish = int(row[1][7])
	for s in SeqIO.parse(seqname, 'fasta'):
		if s.id == cont:
			b1len = maxstart - (begin - 1)							# Determine start buffer
			b2len = maxend - (len(s) - finish)						# Determine end buffer
			buff1 = UnknownSeq(b1len, character = 'N')
			buff2 = UnknownSeq(b2len, character = 'N')
			newseq = str(buff1 + s.seq + buff2)						# Buffer sandwich
			record = SeqRecord(Seq(newseq), id=seqname, description='')
			with open(tmpout, 'a') as f:
				SeqIO.write(record, f, 'fasta')						# Write new sequence to file

print('Writing reverse sequence matches (reverse complement) ...')
for row in seqminus:
	seqname = row[0]
	cont = row[1][0]
	begin = int(row[1][6])
	finish = int(row[1][7])
	for s in SeqIO.parse(seqname, 'fasta'):
		if s.id == cont:
			b1len = maxstart - (len(s) - finish)					# Determine start buffer
			b2len = maxend - (begin - 1)							# Determine end buffer	
			buff1 = UnknownSeq(b1len, character = 'N')
			buff2 = UnknownSeq(b2len, character = 'N')
			rc = s.reverse_complement()
			newseq = str(buff1 + rc.seq + buff2)					# Reverse complement
			record = SeqRecord(Seq(newseq), id=seqname, description='')
			with open(tmpout, 'a') as f:
				SeqIO.write(record, f, 'fasta')						# Write new sequence to file


# Optional annotation
if args.anno == None:
	shutil.move(tmpout, outfile)
	progexit()
else:
	for seq in SeqIO.parse(tmpout, 'fasta'):
		seq.id = os.path.splitext(os.path.basename(seq.id))[0]
		strname = str(seq.id)
		splitfile = 'tempdb/strname'
		banner()
		print('Annotating {} ...'.format(strname))
		banner()
		SeqIO.write(seq, splitfile, 'fasta')
		subprocess.call(['prokka', '--outdir', 'tempdb/split', '--force', '--prefix', strname, '--locustag', strname, '--compliant', '--proteins', args.anno, '--cpus', args.cpus, splitfile])
	import glob
	read_files = glob.glob('tempdb/split/*.gbk')
	with open(outfile, 'w') as write_file:
		for f in read_files:
			with open(f, 'r') as infile:
				write_file.write(infile.read())
	progexit()
