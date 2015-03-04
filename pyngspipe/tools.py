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
	dev = open('/dev/null', 'w')
	head_file(fastq1, "tmp/infer_test_1.fq", 1000000)
	head_file(fastq2, "tmp/infer_test_2.fq", 1000000)
	align_command = "bowtie2 -x {} -1 tmp/infer_test_1.fq -2 tmp/infer_test_2.fq -S tmp/tmp.sam".format(bowtie_ref)
	subprocess.call(align_command.split(), stdout=dev)
	infercommand = "infer_experiment.py -i tmp/tmp.sam -r {} > tmp/infer_res.txt".format(refbed)
	insert = get_insert("tmp/tmp.sam")
	subprocess.call(infercommand, shell=True, stdout=dev)
	reverse = read_infer()
	shutil.rmtree("tmp")
	return reverse, insert

def read_infer():
	with open("tmp/infer_res.txt") as f:
		for line in f:
			line = line.rstrip()
			if line.startswith("Fraction of reads explained by \"1++,1--,2+-,2-+\": "):
				per1 = line.lstrip("Fraction of reads explained by \"1++,1--,2+-,2-+\": ")
			elif line.startswith("Fraction of reads explained by \"1+-,1-+,2++,2--\": "):
				per2 = line.lstrip("Fraction of reads explained by \"1+-,1-+,2++,2--\": ")
	if float(per1) > 0.8:
		infer = "yes"
	elif float(per2) > 0.8:
		infer = "reverse"
	else:
		infer = "no"
	return infer

def paired_rnaseq_process(fastq1, fastq2, gse, gsm, bowtie_ref, gtf, anno_gtf, reverse, insert, threads):
	#Need to look at insert size as well! Add that to infer_experiment?
	print "==> Running Tophat...\n"
	dev = open('/dev/null', 'w')
	align_command = "pyrna_align.py tophat -p {0} {1} -i {2} -g {3} -t {4} -o {5}/{6} -a {7} -b {8}".format(fastq1, fastq2, bowtie_ref, gtf, threads, gse, gsm, 
		int(round(insert[0])), int(round(insert[1])))
	subprocess.call(align_command.split(), stdout=dev)
	print "==> Running HTSeq-count...\n"
	if reverse == "reverse":
		htseq_count = "pyrna_count.py htseq -i {0}/{1}/{1}.bam -g {2} -o {0}/{1}/{1}.count -s reverse".format(gse, gsm, anno_gtf)
	elif reverse == "yes":
		htseq_count = "pyrna_count.py htseq -i {0}/{1}/{1}.bam -g {2} -o {0}/{1}/{1}.count -s yes".format(gse, gsm, anno_gtf)
	elif reverse == "no":
		htseq_count = "pyrna_count.py htseq -i {0}/{1}/{1}.bam -g {2} -o {0}/{1}/{1}.count -s no".format(gse, gsm, anno_gtf)
	subprocess.call(htseq_count.split(), stdout=dev)
	return align_command, htseq_count

def paired_chipseq_process(fastq1, fastq2, gse, gsm, bowtie_ref, genome, threads):
	#For human chipseq, use v1, mouse use v2
	dev = open('/dev/null', 'w')
	if genome == "mm10":
		v = 2
	elif genome == "hg19":
		v = 1
	print "==> Running Bowtie...\n"
	align_command = "pychip_align.py -p {0} {1} -i {2} -v {3} -n {4} -o {5}/{4} -t {6}".format(fastq1, fastq2, bowtie_ref, v, gsm, gse, threads)
	subprocess.call(align_command.split(), stdout=dev)
	print "==> Converting to BigWig...\n"
	toucsc = "pychip_ucsc.py -i {0}/{1}/{1}.sam -g {2} -p".format(gse, gsm, genome)
	subprocess.call(toucsc.split(), stdout=dev)
	return align_command, toucsc

def single_rnaseq_process(fastq, gse, gsm, bowtie_ref, gtf, anno_gtf, threads):
	print "==> Running Tophat...\n"
	dev = open('/dev/null', 'w')
	align_command = "pyrna_align.py tophat -f {0} -i {1} -g {2} -t {3} -o {4}/{5}".format(fastq, bowtie_ref, gtf, threads, gse, gsm)
	subprocess.call(align_command.split(), stdout=dev)
	print "==> Running HTSeq-count...\n"
	htseq_count = "pyrna_count.py htseq -s no -i {0}/{1}/{1}.bam -g {2} -o {0}/{1}/{1}.count".format(gse, gsm, anno_gtf)
	subprocess.call(htseq_count.split(), stdout=dev)
	return align_command, htseq_count

def single_chipseq_process(fastq, gse, gsm, bowtie_ref, genome, threads):
	#For human chipseq, use v1, mouse use v2
	dev = open('/dev/null', 'w')
	if genome == "mm10":
		v = 2
	elif genome == "hg19":
		v = 1
	print "==> Running Bowtie...\n"
	align_command = "pychip_align.py -f {0} -i {1} -v {2} -n {3} -o {4}/{3} -t {5}".format(fastq, bowtie_ref, v, gsm, gse, threads)
	subprocess.call(align_command.split(), stdout=dev)
	print "==> Converting to BigWig...\n"
	toucsc = "pychip_ucsc.py -i {0}/{1}/{1}.sam -g {2}".format(gse, gsm, genome)
	subprocess.call(toucsc.split(), stdout=dev)
	return align_command, toucsc

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
