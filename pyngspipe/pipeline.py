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
import ConfigParser

def rnacleanup(gse, gsm):
	command = "mv {0}/{1}/accepted_hits.bam {0}/{1}/{1}.bam".format(gse, gsm)
	subprocess.call(command.split())
	if os.path.isdir("{}/{}/tmp".format(gse, gsm)):
		shutil.rmtree("{}/{}/tmp".format(gse, gsm))
	if os.path.isdir("{}/{}/logs".format(gse, gsm)):
		shutil.rmtree("{}/{}/logs".format(gse, gsm))
	useless_files = ["prep_reads.info", "unmapped.bam", "trimmed_trun1.fq", "trimmed_trun2.fq", "trimmed_trun_sort.bam", "trimmed_trun_sort.bam.bai"]
	for ufile in useless_files:
		if os.path.isfile("{}/{}/{}".format(gse, gsm, ufile)):
			os.remove("{}/{}/{}".format(gse, gsm, ufile))

def get_paths(config, genome):
	if genome == "hg19":
		bowtie_ref = config["hg19_bowtie_index"]
		gtf = config["hg19_gtf"]
		refbed = config["hg19_ref"]
		anno_gtf = config["hg19_anno_gtf"]
	elif genome == "mm10":
		bowtie_ref = config["mm10_bowtie_index"]
		gtf = config["mm10_gtf"]
		refbed = config["mm10_ref"]
		anno_gtf = config["mm10_anno_gtf"]
	return bowtie_ref, gtf, refbed, anno_gtf

def ConfigSectionMap(section, Config):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

def rnaseq_process_gsm(paired, gse, gsm, gtf, anno_gtf, bowtie_ref, refbed, threads):
	if paired:
		fastq1 = "{0}/{1}/{1}_1.fastq".format(gse, gsm)
		fastq2 = "{0}/{1}/{1}_2.fastq".format(gse, gsm)
		reverse, insert = tools.infer_experiment(fastq1, fastq2, bowtie_ref, refbed)
		align_command, htseq_count = tools.paired_rnaseq_process(fastq1, fastq2, gse, gsm, bowtie_ref, gtf, anno_gtf, reverse, insert, threads)
		rnacleanup(gse, gsm)
		gzip(fastq1)
		gzip(fastq2)
		p = [reverse, insert[0], insert[1]]
		write_pygns_report(gse, gsm, align_command, htseq=htseq_count, paired=p)
	else:
		fastq = "{0}/{1}/{1}.fastq".format(gse, gsm)
		align_command, htseq_count = tools.single_rnaseq_process(fastq, gse, gsm, bowtie_ref, gtf, anno_gtf, threads)
		rnacleanup(gse, gsm)
		gzip(fastq)
		write_pygns_report(gse, gsm, align_command, htseq=htseq_count)

def chipseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, threads, genome):
	if paired:
		fastq1 = "{0}/{1}/{1}_1.fastq".format(gse, gsm)
		fastq2 = "{0}/{1}/{1}_2.fastq".format(gse, gsm)
		align_command, toucsc = tools.paired_chipseq_process(fastq1,fastq2, gse, gsm, bowtie_ref, genome, threads)
		gzip(fastq1)
		gzip(fastq2)
	else:
		fastq = "{0}/{1}/{1}.fastq".format(gse, gsm)
		align_command, toucsc = tools.single_chipseq_process(fastq, gse, gsm, bowtie_ref, genome, threads)
		gzip(fastq)
	write_pygns_report(gse, gsm, align_command, toucsc=toucsc)

def read_tophat_report(gse, gsm):
	align_file = "{}/{}/align_summary.txt".format(gse, gsm)
	data = None
	if os.path.isfile(align_file):
		align = open(align_file, "r")
		lines = align.readlines()
		data = []
		for line in lines:
			line = line.rstrip()
			line = line.lstrip()
			data.append(line)
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

def write_pygns_report(gse, gsm, align_command, htseq=None, toucsc=None, paired=None):
	output = open("{}/{}/pyngs_summary.txt".format(gse, gsm), "w")
	output.write("{}\n".format(align_command)),
	if htseq:
		output.write("{}\n".format(htseq))
	if toucsc:
		output.write("{}\n".format(toucsc))
	if paired:
		output.write("Experiment orientation: {}\nInsert size: {}\nSD: {}\n".format(paired[0], paired[1], paired[2])),
	output.close()

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
	#parser.add_argument('-d', '--db', help='Optional Sqlite database file which will be updated with sample information.', required=False)
	parser.add_argument('-t', '--threads', help='Number of threads to use for alignment, default=1', default=1, required=False)
	#parser.add_argument('-u', action='store_true', help='Will upload files to remote server', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	gsm_dict = {}


	if args["GSE"]:
		gsms = downloader.download_gse(args["GSE"])
		for gsm in sorted(gsms):
			print gsm
			gse, genome, paired, details, sra, exp_type, name = downloader.download_gsm(gsm)
			bowtie_ref, gtf, refbed, anno_gtf = get_paths(config, genome)
			directory = "{}/{}".format(gse, gsm)
			if exp_type == "rnaseq":
				rnaseq_process_gsm(paired, gse, gsm, gtf, anno_gtf, bowtie_ref, refbed, args["threads"])
			elif exp_type == "chipseq":
				chipseq_process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed, args["threads"], genome)
			
		#	if args["db"]:
		#		create_gsm_dict(gsm_dict, gse, args["GSM"], details, sra, genome, "tophat2", exp_type, "PATRICK")
		#		sqlite_scripts.insert_data(args["db"], gsm_dict) #Updating dictionary
			#if args["u"]:
			#	check_dir(directory)
			#	upload_folder(directory, gse, gsm)
			#	shutil.rmtree(directory)

	elif args["GSM"]:
		gse, genome, paired, details, sra, exp_type, name = downloader.download_gsm(args["GSM"]) 
		print gse, genome, paired, details, sra, exp_type, name 
		bowtie_ref, gtf, refbed, anno_gtf = get_paths(config, genome)
		directory = "{}/{}".format(gse, args["GSM"])
		if exp_type == "rnaseq":
			rnaseq_process_gsm(paired, gse, args["GSM"], gtf, anno_gtf, bowtie_ref, refbed, args["threads"])
		elif exp_type == "chipseq":
			chipseq_process_gsm(paired, gse, args["GSM"], gtf, bowtie_ref, refbed, args["threads"], genome)

	#	if args["db"]:
	#		create_gsm_dict(gsm_dict, gse, args["GSM"], details, sra, genome, "tophat2", exp_type, "PATRICK")
	#		sqlite_scripts.insert_data(args["db"], gsm_dict)
		#if args["u"]:
		#	check_dir(directory)
		#	upload_folder(directory, gse, args["GSM"])
		#	shutil.rmtree(directory)