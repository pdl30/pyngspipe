#!/usr/bin/python

########################################################################
# 28 July 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse
import ConfigParser
import time 
import subprocess

def ConfigSectionMap(Config, section):
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

def read_current_spreadsheet(spread):
	#id, GSE, GSM, name, details, srx, genome, aligner, exp_type, alignment_report, date, comments, submitter
	current_samples = {}
	name = spread.split("_")
	version = name[2].strip(".tsv$")
	print version
	with open(spread) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 10:
				current_samples[int(word[0])] = word
	return current_samples, version

def read_idir(idir):
	ifiles = [f for f in os.listdir(idir)]

def find_gsms(conditions):
	gsms = {}
	for key in conditions:
		gsm = os.path.basename(key)
		gsms[key] = gsm
	return gsms

def read_tophat_report(idir):
	align_file = "{}/align_summary.txt".format(idir)
	data = []
	if os.path.isfile(align_file):
		align = open(align_file, "r")
		lines = align.readlines()
		if lines:
			for line in lines:
				line = line.rstrip()
				line = line.lstrip()
				data.append(line)
	return data

def find_details(gsms, conditions):
	gsm_dict = {}
	for path in gsms:
		GSM = gsms[path]
		gsm_dict[GSM] = {}
		date = time.strftime("%d/%m/%Y")
		alignment_report = read_tophat_report(path)
		gse, genome, details, sra, exp_type, name = download_gsm_info(GSM)
		gsm_dict[GSM]["gse"] = gse
		gsm_dict[GSM]["details"] = details
		gsm_dict[GSM]["srx"] = sra
		gsm_dict[GSM]["genome"] = genome
		gsm_dict[GSM]["aligner"] = "tophat2"
		gsm_dict[GSM]["exp_type"] = exp_type
		gsm_dict[GSM]["completion_date"] = date
		gsm_dict[GSM]["alignment_report"] = alignment_report
		gsm_dict[GSM]["submitter"] = "PATRICK"
		gsm_dict[GSM]["name"] = name
	return gsm_dict

def write_new_spreadsheet(current_samples, spread, version, gsm_dict):
	for gsm in gsm_dict:
		if gsm_dict["genome"] == "hg19":
			
	for key in sorted(current_samples):
		output.write("\t".join(current_samples[key])),
	for gsm in gsm_dict:
		if gsm in current_samples:
			print "Conflict with {}\n".format(gsm)
		else:
			output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}")

def download_gsm_info(gsm):
	old_path = os.getcwd()
	download2 = "wget -c -nv -q 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='{0}'&targ=self&view=quick&form=text&token=_blank' -O {0}.soft".format(gsm)
	subprocess.call(download2, shell=True)
	f = open("{}.soft".format(gsm), "r")
	lines = f.readlines()
	extract_info = ""
	growth_info = ""
	os.remove("{}.soft".format(gsm))
	for line in lines:
		line = line.rstrip()
		if line.startswith("!Sample_title = "):
			line = line.lstrip("!Sample_title = ")
			name = line
		if line.startswith("!Sample_series_id = "):
			line = line.lstrip("!Sample_series_id = ")
			gse = line
		if line.startswith("!Sample_organism_ch1 = "):
			line = line.lstrip("!Sample_organism_ch1 = ")
			if line == "Homo sapiens":
				genome = "hg19"
			elif line == "Mus musculus":
				genome = "mm10"
		if line.startswith("!Sample_library_strategy = "):
			line = line.lstrip("!Sample_library_strategy = ")
			if line == "RNA-Seq":
				exp_type = "rnaseq"
			elif line == "ChIP-Seq":
				exp_type = "chipseq"
		
		if line.startswith("!Sample_extract_protocol_ch1"):
			extract_info = line.lstrip("!Sample_extract_protocol_ch1 = ")
			extract_info = extract_info.rstrip()
		elif line.startswith("!Sample_growth_protocol_ch1"):
			growth_info = line.lstrip("!Sample_growth_protocol_ch1 = ")
			growth_info = growth_info.rstrip()
	details = extract_info + growth_info
	sra = None
	gsm_sra = []
	for line in lines:
		line = line.rstrip()
		if line.startswith("!Sample_supplementary_file"):
			sra_path = line.split("ByExp/sra/")
			if len(sra_path) > 1:
				sra = sra_path[1]
	return gse, genome, details, sra, exp_type, name

def main():
	parser = argparse.ArgumentParser(description='Updates current spreadsheet with processed samples.')
	parser.add_argument('-c', '--config', help='Config pointing to finished samples', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap(Config, "Conditions")
	gsms = find_gsms(conditions)
	gsm_dict = find_details(gsms, conditions)
	print gsm_dict

main()