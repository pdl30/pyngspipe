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
import tarfile
import shutil

upload_folder = "/home/pdl30/RNAseq_compendium/samples"
upload_server = "super:"

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

def make_tarfile(output_filename, source_dir):
	ifiles = [f for f in os.listdir(source_dir)]
	with tarfile.open(output_filename, "w:gz") as tar:
		tar.add(source_dir, arcname=os.path.basename(source_dir))

def upload_folder(idir):
	command = "scp -r {} {}{}".format(idir, upload_server, upload_folder)
	subprocess.call(command.split())

#You need to zip and tar the entire directory and then upload it to tobias
#def main():
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Pyrnapipe is a RNA-seq pipeline. \n')
	parser.add_argument('-i', '--idir', help='Input directory', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	check_dir(args["idir"])
	#name = os.path.normpath(args["idir"])
	#name = name+".tar.gz"
	#make_tarfile(name, args["idir"])
	upload_folder(args["idir"])
	shutil.rmtree(args["idir"])