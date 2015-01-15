#!/usr/bin/python

########################################################################
# 9 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import argparse
import os, re, sys
import subprocess
import shutil
import pyngspipe
from pyngspipe import tools, downloader, sqlite_scripts
import pkg_resources
import time

def rnacleanup(gse, gsm):
	command = "mv {0}/{1}/accepted_hits.bam {0}/{1}/{1}.bam".format(gse, gsm)
	subprocess.call(command.split())
	if os.path.isdir("{}/{}/tmp".format(gse, gsm)):
		shutil.rmtree("{}/{}/tmp".format(gse, gsm))
	if os.path.isdir("{}/{}/logs".format(gse, gsm)):
		shutil.rmtree("{}/{}/logs".format(gse, gsm))
	useless_files = ["prep_reads.info", "tophat_report.txt", "unmapped.bam"]
	for ufile in useless_files:
		if os.path.isfile("{}/{}/{}".format(gse, gsm, ufile)):
			os.remove("{}/{}/{}".format(gse, gsm, ufile))

def get_paths(path1, genome):
	if genome == "hg19":
		bowtie_ref = "/home/patrick/Reference_Genomes/pyrnapipe_references/hg19/hg19"
		gtf = path1 + "hg19.gtf"
		refbed = path1 + "hg19_Ensembl.bed"
	elif genome == "mm10":
		bowtie_ref = "/home/patrick/Reference_Genomes/pyrnapipe_references/mm10/mm10"
		gtf = path1 + "mm10.gtf"
		refbed = path1 + "mm10_Ensembl.bed"
	return bowtie_ref, gtf, refbed

def rnaseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, threads):
	if paired:
		fastq1 = "{0}/{1}/{1}_1.fastq".format(gse, gsm)
		fastq2 = "{0}/{1}/{1}_2.fastq".format(gse, gsm)
		reverse, insert = tools.infer_experiment(fastq1, fastq2, bowtie_ref, refbed)
		tools.paired_rnaseq_process(fastq1, fastq2, gse, gsm, bowtie_ref, gtf, reverse, insert, threads)
		rnacleanup(gse, gsm)
		gzip(fastq1)
		gzip(fastq2)
	else:
		fastq = "{0}/{1}/{1}.fastq".format(gse, gsm)
		tools.single_rnaseq_process(fastq, gse, gsm, bowtie_ref, gtf)
		rnacleanup(gse, gsm)
		gzip(fastq)

def chipseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, threads):
	if paired:
		fastq1 = "{0}/{1}/{1}_1.fastq".format(gse, gsm)
		fastq2 = "{0}/{1}/{1}_2.fastq".format(gse, gsm)
		#reverse, insert = tools.infer_experiment(fastq1, fastq2, bowtie_ref, refbed)
		tools.paired_chipseq_process(fastq, gse, gsm, bowtie_ref, genome, threads)
		gzip(fastq1)
		gzip(fastq2)
	else:
		fastq = "{0}/{1}/{1}.fastq".format(gse, gsm)
		tools.single_chipseq_process(fastq, gse, gsm, bowtie_ref, genome, threads)
		gzip(fastq)

def read_tophat_report(gse, gsm):
	align_file = "{}/{}/align_summary.txt".format(gse, gsm)
	data = None
	if os.path.isfile(align_file):
		align = open(align_file, "r")
		lines = align.readlines()
		if lines:
			data = lines
	return data

def read_bowtie_report(gse, gsm):
	align_file = "{0}/{1}/{1}_report.txt".format(gse, gsm)
	data = None
	if os.path.isfile(align_file):
		align = open(align_file, "r")
		lines = align.readlines()
		if lines:
			data = lines
	return data

def gzip(ifile):
	if os.path.isfile(ifile):
		command = "gzip {}".format(ifile)
		subprocess.call(command.split())

def create_gsm_dict(gsm_dict, GSE, GSM, details, srx, genome, aligner, exp_type, submitter):
	gsm_dict[GSM] = {}
	date = time.strftime("%d/%m/%Y")
	if exp_type == "rnaseq":
		alignment_report = read_tophat_report(GSE, GSM)
	elif exp_type == "chipseq":
		alignment_report = read_bowtie_report(GSE, GSM)
	gsm_dict[GSM]["gse"] = GSE
	gsm_dict[GSM]["details"] = details
	gsm_dict[GSM]["srx"] = srx
	gsm_dict[GSM]["genome"] = genome
	gsm_dict[GSM]["aligner"] = aligner
	gsm_dict[GSM]["exp_type"] = exp_type
	gsm_dict[GSM]["completion_date"] = date
	gsm_dict[GSM]["alignment_report"] = alignment_report
	gsm_dict[GSM]["submitter"] = submitter
	gsm_dict[GSM]["paired"] = paired
	return gsm_dict

def main():
	parser = argparse.ArgumentParser(description='Pyrnapipe is a RNA-seq pipeline. \n')
	parser.add_argument('-g', '--GSE', help='GSE accession for processing. Will try all samples in the accession', required=False)
	parser.add_argument('-m', '--GSM', help='Individual GSM samples for processing', required=False)
	parser.add_argument('-d', '--db', help='Optional Sqlite database file which will be updated with sample information.', required=False)
	parser.add_argument('-t', '--threads', help='Number of threads to use for alignment, default=1', default=1, required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	
	gsm_dict = {}
	
	path1 = pkg_resources.resource_filename('pyngspipe', 'data/')

	if args["GSE"]:
		gsms = downloader.download_gse(args["GSE"])
		for gsm in sorted(gsms):
			gse, genome, paired, details, sra, exp_type = downloader.download_gsm(gsm)
			bowtie_ref, gtf, refbed = get_paths(path1, genome)
			if exp_type == "rnaseq":
				rnaseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, args["threads"])
			elif exp_type == "chipseq":
				chipseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, args["threads"])
			
			if args["db"]:
				create_gsm_dict(gsm_dict, gse, args["GSM"], details, sra, genome, "tophat2", exp_type, "PATRICK")
				sqlite_scripts.insert_data(args["db"], gsm_dict) #Updating dictionary

	elif args["GSM"]:
		gse, genome, paired, details, sra, exp_type = downloader.download_gsm(args["GSM"]) 
		bowtie_ref, gtf, refbed = get_paths(path1, genome)
		if exp_type == "rnaseq":
			rnaseq_process_gsm(paired, gse, args["GSM"], gtf, bowtie_ref, refbed, args["threads"])
		elif exp_type == "chipseq":
			chipseq_process_gsm(paired, gse, args["GSM"], gtf, bowtie_ref, refbed, args["threads"])

		if args["db"]:
			create_gsm_dict(gsm_dict, gse, args["GSM"], details, sra, genome, "tophat2", exp_type, "PATRICK")
			sqlite_scripts.insert_data(args["db"], gsm_dict)
		

