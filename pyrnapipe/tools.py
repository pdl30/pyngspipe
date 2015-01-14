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
import math

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
	align_command = "bowtie2 -x {} -1 tmp/infer_test_1.fq -2 tmp/infer_test_2.fq -S tmp/tmp.sam".format(bowtie_ref)
	subprocess.call(align_command.split())
	infercommand = "infer_experiment.py -i tmp/tmp.sam -r {} > tmp/infer_res.txt".format(refbed)
	insert = get_insert("tmp/tmp.sam")
	subprocess.call(infercommand, shell=True)
	reverse = read_infer()
	shutil.rmtree("tmp")
	return reverse, insert

def read_insert():
	insert = []
	with open("tmp/insert_res.txt") as f:
		c = 0
		for line in f:
			line = line.rstrip()
			if c == 1:
				insert.append(float(line))
			elif c == 5:
				insert.append(float(line))
			c += 1
	return insert

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

def paired_process(fastq1, fastq2, gse, gsm, bowtie_ref, gtf, reverse, insert):
	#Need to look at insert size as well! Add that to infer_experiment?
	align_command = "pyrna_align.py tophat -p {} {} -i {} -g {} -t 2 -o {}/{} -a {} -b {}".format(fastq1, fastq2, bowtie_ref, gtf, gse, gsm, insert[0], insert[1])
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

def getmeanval(dic,maxbound=-1):
	nsum=0;  n=0;	
	for (k,v) in dic.items():
		if maxbound!=-1 and k>maxbound:
			continue
		nsum=nsum+k*v;
		n=n+v;
	meanv=nsum*1.0/n;
	nsum=0; n=0;
	for (k,v) in dic.items():
		if maxbound!=-1 and k>maxbound:
			continue;
		nsum=nsum+(k-meanv)*(k-meanv)*v;
		n=n+v;
		varv=math.sqrt(nsum*1.0/(n-1));
	return [meanv,varv];
 
def get_insert(samfile):
	plrdlen={};
	plrdspan={};
	objmrl=re.compile('([0-9]+)M$');
	objmtj=re.compile('NH:i:(\d+)');
	nline=0
	sam = open(samfile, "r")
	for lines in sam:
		field=lines.strip().split();
		nline=nline+1;
		if len(field)<12:
			continue;
		try:
			mrl=objmrl.match(field[5]);
			if mrl==None: # ignore non-perfect reads
				continue;
			readlen=int(mrl.group(1));
			if readlen in plrdlen.keys():
				plrdlen[readlen]=plrdlen[readlen]+1;
			else:
				plrdlen[readlen]=1;
			if field[6]!='=':
				continue;
			dist=int(field[8]);
			if dist<=0: # ignore neg dist
				continue;
			mtj=objmtj.search(lines);

			if dist in plrdspan.keys():
				plrdspan[dist]=plrdspan[dist]+1;
			else:
				plrdspan[dist]=1;
		except ValueError:
			continue;

	if len(plrdspan)==0:
		print('No qualified paired-end reads found. Are they single-end reads?');
	else:
		maxv=max(plrdspan,key=plrdspan.get);
		spanval=getmeanval(plrdspan,maxbound=maxv*3);
		return spanval
