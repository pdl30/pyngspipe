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

def download_gse(gse):
	if not os.path.isdir(gse):
		os.mkdir(gse)
	gsm_samples = {}
	download2 = "wget -c -nv -q 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='{0}'&targ=self&view=quick&form=text&token=_blank' -O {0}.soft".format(gse)
	subprocess.call(download2, shell=True)
	with open("{}.soft".format(gse)) as f:
		for line in f:
			line  = line.rstrip()
			if line.startswith("!Series_sample_id = "):
				gsm = re.sub("!Series_sample_id = ", "", line)
				gsm_samples[gsm] = 1
	os.remove("{}.soft".format(gse))
	return gsm_samples

def download_gsm(gsm):
	old_path = os.getcwd()
	download2 = "wget -c -nv -q 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='{0}'&targ=self&view=quick&form=text&token=_blank' -O {0}.soft".format(gsm)
	subprocess.call(download2, shell=True)
	f = open("{}.soft".format(gsm), "r")
	lines = f.readlines()
	os.remove("{}.soft".format(gsm))
	for line in lines:
		line = line.rstrip()
		if line.startswith("!Sample_series_id = "):
			line = line.lstrip("!Sample_series_id = ")
			gse = line
		if line.startswith("!Sample_organism_ch1 = "):
			line = line.lstrip("!Sample_organism_ch1 = ")
			if line == "Homo sapiens":
				genome = "hg19"
			elif line == "Mus musculus":
				genome = "mm10"
		if line.startswith("!Sample_extract_protocol_ch1 = "):
			extract_info = line.lstrip("!Sample_organism_ch1 = ")
		elif line.startswith("!Sample_growth_protocol_ch1 = "):
			growth_info = line.lstrip("!Sample_growth_protocol_ch1 = ")
		details = extract_info + ":::" + growth_info
	if not os.path.isdir(gse):
		os.mkdir(gse)
	gsm_path = "{}/{}".format(gse, gsm)
	if not os.path.isdir(gsm_path):
		os.mkdir(gsm_path)
	os.chdir(gsm_path)
	sra = None
	gsm_sra = []
	for line in lines:
		line = line.rstrip()
		if line.startswith("!Sample_supplementary_file"):
			sra_path = line.split("ByExp/sra/")
			if len(sra_path) > 1:
				sra = sra_path[1]
				download3 = "wget -r -c --no-verbose -N -nd ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/{}".format(sra)
				subprocess.call(download3, shell=True)
	paired = combine_convert(old_path, gsm)
	return gse, genome, paired, details, sra

def combine_convert(old_path, gsm):
	new_path = os.getcwd()
	sras = [f for f in os.listdir(new_path) if f.endswith(".sra")]
	for sra in sras:
		command0 = "fastq-dump --split-3 {}".format(sra)
		subprocess.call(command0, shell=True)
		os.remove(sra)
	#Catting everything together
	fq1s = [f for f in os.listdir(new_path) if f.startswith("SRR") and f.endswith("_1.fastq")]
	if fq1s:
		paired = True
		command = "cat"
		for fq in fq1s:
			command += " {}".format(fq)
		command += " > {}_1.fastq".format(gsm)
		subprocess.call(command, shell=True)
		fq2s = [f for f in os.listdir(new_path) if f.startswith("SRR") and f.endswith("_2.fastq")]
		command = "cat"
		for fq in fq2s:
			command += " {}".format(fq)
		command += " > {}_2.fastq".format(gsm)
		subprocess.call(command, shell=True)
		for fq in fq1s:
			os.remove(fq)
		for fq in fq2s:
			os.remove(fq)
	else:
		paired = False
		fqs = [f for f in os.listdir(new_path) if f.startswith("SRR") and f.endswith(".fastq")]
		command = "cat"
		for fq in fqs:
			command += " {}".format(fq)
		command += " > {}.fastq".format(gsm)
		subprocess.call(command, shell=True)
		for fq in fqs:
			os.remove(fq)
	#Remove old SRR fastqs
	os.chdir(old_path)
	return paired