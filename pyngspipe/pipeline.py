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
	useless_files = ["prep_reads.info", "tophat_report.txt", "unmapped.bam", "trimmed_trun1.fq", "trimmed_trun2.fq", "trimmed_trun_sort.bam", "trimmed_trun_sort.bam.bai"]
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
		tools.single_rnaseq_process(fastq, gse, gsm, bowtie_ref, gtf, threads)
		rnacleanup(gse, gsm)
		gzip(fastq)

def chipseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, threads, genome):
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
		command = "gzip -f {}".format(ifile)
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

def check_dir(idir):
	necessary_ends = {"bam":False, "txt":False, "count":False, "fastqc.zip":False, "fastq.gz":False}
	ifiles = [f for f in os.listdir(idir)]
	for end in necessary_ends:
		for ifile in ifiles:
			if ifile.endswith(end):
				necessary_ends[end] = True
	for ends in necessary_ends:
		if necessary_ends[ends] == False:
			raise Exception("Not all files present")

def upload_folder(idir, gse, gsm):
	command = "rsync -ave ssh {} super:/home/pdl30/RNAseq_compendium/samples/{}/".format(idir, gse)
	subprocess.call(command.split())

def main():
	parser = argparse.ArgumentParser(description='Pyrnapipe is a pipeline for RNA-seq and ChIP-seq samples from the GEO database\n')
	parser.add_argument('-g', '--GSE', help='GSE accession for processing. Will try all samples in the accession', required=False)
	parser.add_argument('-m', '--GSM', help='Individual GSM samples for processing', required=False)
	parser.add_argument('-d', '--db', help='Optional Sqlite database file which will be updated with sample information.', required=False)
	parser.add_argument('-t', '--threads', help='Number of threads to use for alignment, default=1', default=1, required=False)
	parser.add_argument('-u', action='store_true', help='Will upload files to remote server', required=False)
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
			directory = "{}/{}".format(gse, gsm)
			if exp_type == "rnaseq":
				rnaseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, args["threads"])
			elif exp_type == "chipseq":
				chipseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, args["threads"], genome)
			
			if args["db"]:
				create_gsm_dict(gsm_dict, gse, args["GSM"], details, sra, genome, "tophat2", exp_type, "PATRICK")
				sqlite_scripts.insert_data(args["db"], gsm_dict) #Updating dictionary
			if args["u"]:
				check_dir(directory)
				upload_folder(directory, gse, gsm)
				shutil.rmtree(directory)

	elif args["GSM"]:
		gse, genome, paired, details, sra, exp_type = downloader.download_gsm(args["GSM"]) 
		bowtie_ref, gtf, refbed = get_paths(path1, genome)
		directory = "{}/{}".format(gse, args["GSM"])
		if exp_type == "rnaseq":
			rnaseq_process_gsm(paired, gse, args["GSM"], gtf, bowtie_ref, refbed, args["threads"])
		elif exp_type == "chipseq":
			chipseq_process_gsm(paired, gse, args["GSM"], gtf, bowtie_ref, refbed, args["threads"], genome)

		if args["db"]:
			create_gsm_dict(gsm_dict, gse, args["GSM"], details, sra, genome, "tophat2", exp_type, "PATRICK")
			sqlite_scripts.insert_data(args["db"], gsm_dict)
		if args["u"]:
			check_dir(directory)
			upload_folder(directory, gse, args["GSM"])
			shutil.rmtree(directory)
