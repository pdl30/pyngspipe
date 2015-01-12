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
import pkg_resources

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
	return gse, genome, paired

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

def cleanup(gse, gsm):
	command = "mv {0}/{1}/accepted_hits.bam {0}/{1}/{1}.bam".format(gse, gsm)
	subprocess.call(command.split())
	if os.path.isdir("{}/{}/tmp".format(gse, gsm)):
		shutil.rmtree("{}/{}/tmp".format(gse, gsm))
	if os.path.isdir("{}/{}/logs".format(gse, gsm)):
		shutil.rmtree("{}/{}/logs".format(gse, gsm))

def head_file(ifile, outfile, num):
	output = open(outfile, "w")
	c = 0
	with open(ifile) as f:
		for line in f:
			if c < num:
				output.write(line),
			else:
				break
			c += 1

def infer_experiment(fastq1, fastq2, bowtie_ref, refbed):
	if not os.path.isdir("tmp"):
		os.mkdir("tmp")
	head_file(fastq1, "tmp/infer_test_1.fq", 1000000)
	head_file(fastq2, "tmp/infer_test_2.fq", 1000000)
	align_command = " bowtie2 -x {} -1 tmp/infer_test_1.fq -2 tmp/infer_test_2.fq -S tmp/tmp.sam".format(bowtie_ref)
	subprocess.call(align_command.split())
	infercommand = "infer_experiment.py -i tmp/tmp.sam -r {} > tmp/infer_res.txt".format(refbed)
	subprocess.call(infercommand, shell=True)
	reverse = read_infer()
	shutil.rmtree("tmp")
	return reverse

def read_infer():
	with open("tmp/infer_res.txt") as f:
		for line in f:
			line = line.rstrip()
			if line.startswith("Fraction of reads explained by \"1++,1--,2+-,2-+\": "):
				per1 = line.lstrip("Fraction of reads explained by \"1++,1--,2+-,2-+\": ")
			elif line.startswith("Fraction of reads explained by \"1+-,1-+,2++,2--\": "):
				per2 = line.lstrip("Fraction of reads explained by \"1+-,1-+,2++,2--\": ")
	if float(per1) > float(per2):
		reverse = False
	else:
		reverse = True
	return reverse

def paired_process(fastq1, fastq2, gse, gsm, bowtie_ref, gtf, reverse):
	#Need to look at insert size as well! Add that to infer_experiment?
	align_command = "pyrna_align.py tophat -p {} {} -i {} -g {} -t 2 -o {}/{}".format(fastq1, fastq2, bowtie_ref, gtf, gse, gsm)
	subprocess.call(align_command.split())
	#toucsc = "pyrna_ucsc.py -i {}/{}/accepted_hits.bam -g hg19 -ens".format(gse, gsm)
	#subprocess.call(toucsc.split())
	if reverse:
		htseq_count = "pyrna_count.py htseq -i {0}/{1}/accepted_hits.bam -g {2} -o {0}/{1}/{1}.count -s reverse".format(gse, gsm, gtf)
	else:
		htseq_count = "pyrna_count.py htseq -i {0}/{1}/accepted_hits.bam -g {2} -o {0}/{1}/{1}.count -s yes".format(gse, gsm, gtf)
	subprocess.call(htseq_count.split())

def single_process(fastq, gse, gsm, bowtie_ref, gtf):
	align_command = "pyrna_align.py tophat -f {} -i {} -g {} -t 2 -o {}/{}".format(fastq, bowtie_ref, gtf, gse, gsm)
	subprocess.call(align_command.split())
	#toucsc = "pyrna_ucsc.py -i {}/{}/accepted_hits.bam -g hg19 -ens".format(gse, gsm)
	#subprocess.call(toucsc.split())
	htseq_count = "pyrna_count.py htseq -i {0}/{1}/accepted_hits.bam -g {2} -o {0}/{1}/{1}.count".format(gse, gsm, gtf)
	subprocess.call(htseq_count.split())	

def get_paths(path1, genome):
	if genome == "hg19":
		bowtie_ref = "/home/patrick/Reference_Genomes/pyrnapipe_references/hg19/hg19"
		gtf = path2 + "hg19.gtf"
		refbed = path2 + "hg19_Ensembl.bed"
	elif genome == "mm10":
		bowtie_ref = "/home/patrick/Reference_Genomes/pyrnapipe_references/mm10/mm10"
		gtf = path2 + "mm10.gtf"
		refbed = path2 + "mm10_Ensembl.bed"
	return bowtie_ref, gtf, refbed

def process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed):
	if paired:
		fastq1 = "{0}/{1}/{1}_1.fastq".format(gse, gsm)
		fastq2 = "{0}/{1}/{1}_2.fastq".format(gse, gsm)
		infer_experiment(fastq1, fastq2, bowtie_ref, refbed)
		paired_process(fastq1, fastq2, gse, gsm, bowtie_ref, gtf)
	else:
		fastq = "{0}/{1}/{1}.fastq".format(gse, gsm)
		single_process(fastq, gse, gsm, bowtie_ref, gtf)
	cleanup(gse, gsm)

def main():
	parser = argparse.ArgumentParser(description='Pyrnapipe is a RNA-seq pipeline. \n')
	parser.add_argument('-g', '--GSE', help='GSE accession for processing. Will try all samples in the accession', required=False)
	parser.add_argument('-m', '--GSM', help='Individual GSM samples for processing', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	path1 = pkg_resources.resource_filename('pyrnapipe', 'data/')
	if args["GSE"]:
		gsms = download_gse(args["GSE"])
		for gsm in gsms:
			gse, genome, paired = download_gsm(gsm)
			bowtie_ref, gtf, refbed = get_paths(path1, path2, genome)
			process_gsm(paired, gse, gsm, gtf, bowtie_ref, refbed)

	elif args["GSM"]:
		gse, genome, paired = download_gsm(args["GSM"])
		bowtie_ref, gtf, refbed = get_paths(path1, path2, genome)
		process_gsm(paired, gse, args["GSM"], gtf, bowtie_ref, refbed)
