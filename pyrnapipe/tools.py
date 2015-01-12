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