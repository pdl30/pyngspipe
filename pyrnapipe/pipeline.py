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
import pyrnapipe
from pyrnapipe import tools, downloader, sqlite_scripts
import pkg_resources
import time

def cleanup(gse, gsm):
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

def process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed):
	if paired:
		fastq1 = "{0}/{1}/{1}_1.fastq".format(gse, gsm)
		fastq2 = "{0}/{1}/{1}_2.fastq".format(gse, gsm)
		reverse, insert = tools.infer_experiment(fastq1, fastq2, bowtie_ref, refbed)
		tools.paired_process(fastq1, fastq2, gse, gsm, bowtie_ref, gtf, reverse, insert)
		cleanup(gse, gsm)
	else:
		fastq = "{0}/{1}/{1}.fastq".format(gse, gsm)
		tools.single_process(fastq, gse, gsm, bowtie_ref, gtf)
		cleanup(gse, gsm)

def read_alignment_report(gse, gsm):
	align_file = "{}/{}/align_summary.txt".format(gse, gsm)
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

def create_gsm_dict(gsm_dict, GSE, GSM, details, srx, genome, aligner, submitter):
	gsm_dict[gsm] = {}
	date = time.strftime("%d/%m/%Y")
	alignment_report = read_alignment_report(gse, gsm)
	gsm_dict[gsm]["gse"] = GSE
	gsm_dict[gsm]["details"] = details
	gsm_dict[gsm]["srx"] = srx
	gsm_dict[gsm]["genome"] = genome
	gsm_dict[gsm]["aligner"] = aligner
	gsm_dict[gsm]["completion_date"] = date
	gsm_dict[gsm]["alignment_report"] = alignment_report
	gsm_dict[gsm]["submitter"] = submitter
	gsm_dict[gsm]["paired"] = paired
	return gsm_dict

def main():
	parser = argparse.ArgumentParser(description='Pyrnapipe is a RNA-seq pipeline. \n')
	parser.add_argument('-g', '--GSE', help='GSE accession for processing. Will try all samples in the accession', required=False)
	parser.add_argument('-m', '--GSM', help='Individual GSM samples for processing', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	
	path1 = pkg_resources.resource_filename('pyrnapipe', 'data/')
	sqlite_database = "/home/patrick/Scripts/pyrnapipe/database/"
	gsm_dict = {}

	if args["GSE"]:
		gsms = downloader.download_gse(args["GSE"])
		for gsm in gsms:
			gse, genome, paired, details, sra = downloader.download_gsm(gsm)
			bowtie_ref, gtf, refbed = get_paths(path1, genome)
			process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed)
			create_gsm_dict(gsm_dict, gse, args["GSM"], details, sra, genome, "tophat2","PATRICK")
			insert_data(sqlite_database, gsm_dict) #Updating dictionary

	elif args["GSM"]:
		gse, genome, paired, details, sra = downloader.download_gsm(args["GSM"]) 
		bowtie_ref, gtf, refbed = get_paths(path1, genome)
		process_gsm(paired, gse, args["GSM"], gtf, bowtie_ref, refbed)
		create_gsm_dict(gsm_dict, gse, args["GSM"], details, sra, genome, "tophat2","PATRICK")
		insert_data(sqlite_database, gsm_dict)
