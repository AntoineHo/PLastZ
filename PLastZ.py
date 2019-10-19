#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import errno
import sys
import shutil
import subprocess
import datetime
import argparse
from multiprocessing import Pool, TimeoutError
import subprocess

from Bio import SeqIO # Need BIOPYTHON SEQ/IO

class FH :
	def __init__(self, query, target, outdir) :
		self.outdir = outdir
		self.query = None
		self.target = None
		self.tmpdir = None
		if os.path.isfile(query) :
			self.query = os.path.abspath(query)
		else :
			raise Exception("ERROR: Query .fa file does not exist!")

		if os.path.isfile(target) :
			self.target = os.path.abspath(target)
		else :
			raise Exception("ERROR: Target .fa file does not exist!")

		self.make_outdir()

	def make_outdir(self) :
		try :
			os.makedirs(self.outdir)
			self.outdir = os.path.abspath(self.outdir)
		except OSError as e :
			if e.errno != errno.EEXIST :
				raise Exception("ERROR: Cannot create output directory!")
			else :
				self.outdir = os.path.abspath(self.outdir)
				print("WARNING: Directory already exists!")

		try :
			self.tmpdir = os.path.join(self.outdir, "tmp")
			os.makedirs(self.tmpdir)
		except OSError as e :
			if e.errno != errno.EEXIST :
				raise Exception("ERROR: Cannot create output directory!")

	def remove_tmp(self) :
		shutil.rmtree(self.tmpdir)

	def __str__(self) :
		return "Query: {}\nTarget: {}\nOut directory: {}".format(self.query, self.target, self.outdir)

def str_to_bool(v) :
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

def parseArgs() :
	parser = argparse.ArgumentParser(description='Parallel LastZ of a query against a target (can handle self).')
	parser.add_argument('Query',nargs=1,type=str,help="A fasta file containing the target sequences")
	parser.add_argument('Target',nargs=1,type=str,help="A fasta file containing the query sequences")
	parser.add_argument('Output',nargs=1,type=str,help="A directory path")
	parser.add_argument('--lastz-options','-lo',nargs=1,type=str,default=[None],required=False,help="A list of arguments e.g.: \"--strand=plus --noxtrim\"")
	parser.add_argument('--processes','-p',nargs=1,type=int,default=[4],required=False,help="Maximum threads to use. Default: 4.")
	parser.add_argument('--keep-temp', '-kp',type=str_to_bool, nargs='?', const=True, default=False, help="Keep temporary files (single alignments and single fasta). Unset by default.")
	return parser.parse_args()

def create_jobs(query, target, tmpdir, lastz_options) :
	extract_job_list = []
	align_job_list = []
	extracted_sequences = []
	pairs_done = []

	qCtgs = []
	for qry in SeqIO.parse(query, "fasta") :
		qCtgs.append(qry.id)
	tCtgs = []
	for tgt in SeqIO.parse(target, "fasta") :
		tCtgs.append(tgt.id)

	for qry in qCtgs :
		if qry not in extracted_sequences :
			ex = "samtools faidx {} {} > {}".format(query, qry, os.path.join(tmpdir, qry + ".fa") )
			extract_job_list.append(ex)
			extracted_sequences.append(qry)

		for tgt in tCtgs :
			if [qry, tgt] in pairs_done :
				continue
			pairs_done.append([qry, tgt])
			pairs_done.append([tgt, qry])

			if tgt not in extracted_sequences :
				ex = "samtools faidx {} {} > {}".format(target, tgt, os.path.join(tmpdir, tgt + ".fa") )
				extract_job_list.append(ex)
				extracted_sequences.append(tgt)

			tmpfile = os.path.join(tmpdir, "{}_V_{}.TMP".format(qry, tgt))
			qrypath = os.path.join(tmpdir, qry + ".fa")
			tgtpath = os.path.join(tmpdir, tgt + ".fa")
			cmd = ""
			if lastz_options != None :
				cmd = "lastz {} {} {} > {}".format(qrypath, tgtpath, lastz_options, tmpfile)
			else :
				cmd = "lastz {} {} > {}".format(qrypath, tgtpath, tmpfile)
			align_job_list.append(cmd)

	return extract_job_list, align_job_list

def run(cmd) :
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	proc.communicate()

def job_runner(job_list, proc) :
	p = Pool(processes=proc)
	p.map(run, job_list)
	return 1

def main() :
	args = parseArgs()
	query = args.Query[0]
	target = args.Target[0]
	out = args.Output[0]
	lastz_options = args.lastz_options[0] # LastZ options
	proc = args.processes[0]
	keep_temp = args.keep_temp
	iFH = FH(query, target, out) # File Handler
	print("Prepping files...")
	extract_job_list, align_job_list = create_jobs(iFH.query, iFH.target, iFH.tmpdir, lastz_options)
	print("Extracting queries and targets sequences...")
	job_runner(extract_job_list, proc)
	print("Aligning sequences...")
	job_runner(align_job_list, proc)
	print("Concatenate alignments...")
	run("cat {}/*.TMP > {}/plastz.out".format(iFH.tmpdir, iFH.outdir))
	if not keep_temp :
		iFH.remove_tmp()

	print("Done")
	sys.exit(0)

if __name__ == '__main__':
	main()
