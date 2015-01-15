#!/usr/bin/python

########################################################################
# 12 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import re, os, sys
import argparse
import sqlite3 as lite
import pyrnapipe
import pkg_resources
import time

def create_database(sqlite_db):
	con = None
	try:
	#	if os.path.isfile(pkg_path + "/GEO_samples.db")
		con = lite.connect(sqlite_db)
		with con:
			cur = con.cursor()   
			cur.execute("DROP TABLE IF EXISTS Samples") 
			cur.execute("CREATE TABLE Samples(Id INTEGER PRIMARY KEY AUTOINCREMENT, GSE TEXT, GSM TEXT, details TEXT, paired TEXT, srx TEXT, genome TEXT, aligner TEXT, exp_type TEXT completion_date TEXT, alignment_report TEXT, submitter TEXT)")
	except lite.Error, e:
		if con:
			con.rollback()
			print "Error %s:" % e.args[0]
			sys.exit(1)
	finally:
		if con:
			con.close() 

def insert_data(sqlite_db, gsm_dict):
	con = None
	try:
		if os.path.isfile(sqlite_db):
			con = lite.connect(sqlite_db)
		else:
			create_database(sqlite_db)
			con = lite.connect(sqlite_db)
		with con:
			cur = con.cursor()
			for gsm in gsm_dict:
				cur.execute("INSERT INTO Samples(GSE, GSM, details, paired, srx, genome, aligner, exp_type, completion_date, alignment_report, submitter) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", 
					(gsm_dict[gsm]["gse"], gsm, gsm_dict[gsm]["details"], gsm_dict[gsm]["paired"], gsm_dict[gsm]["srx"], gsm_dict[gsm]["genome"], gsm_dict[gsm]["aligner"], gsm_dict[gsm]["exp_type"], 
						gsm_dict[gsm]["date"], gsm_dict[gsm]["report"],gsm_dict[gsm]["submitter"]))
			con.commit()
	except lite.Error, e:
		if con:
			con.rollback()
			print "Error %s:" % e.args[0]
			sys.exit(1)
	finally:
		if con:
			con.close() 

def print_database(sqlite_db):
	con = None
	try:
	#	if os.path.isfile(pkg_path + "/GEO_samples.db"):
		con = lite.connect(sqlite_db)
	#	else:
	#		create_database(pkg_path)
		with con:
			cur = con.cursor()
			cur.execute("SELECT * FROM Samples")
			rows = cur.fetchall()
			for row in rows:
				print row
	except lite.Error, e:
		if con:
			con.rollback()
			print "Error %s:" % e.args[0]
			sys.exit(1)
	finally:
		if con:
			con.close() 

def print_schmena(sqlite_db):
	con = None
	try:
		if os.path.isfile(sqlite_db):
			con = lite.connect(sqlite_db)
		else:
			create_database(sqlite_db)
			con = lite.connect(sqlite_db)
		with con:
			cur = con.cursor()
			cur.execute('PRAGMA table_info(Samples)')
			data = cur.fetchall()
			print data
	except lite.Error, e:
		if con:
			con.rollback()
			print "Error %s:" % e.args[0]
			sys.exit(1)
	finally:
		if con:
			con.close() 

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Sqlite scripts for use with pyrnapipe\n')
	parser.add_argument('-d', '--db', help='Sqlite_database file to be used', required=False)
	parser.add_argument('-c', help='Create a database', action='store_true')
	parser.add_argument('-t', help='Test insert', action='store_true')
	parser.add_argument('-p', help='Print database', action='store_true')
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	if args["c"]:
		create_database(args["db"])
	elif args["t"]:
		date = time.strftime("%d/%m/%Y")
		gsm_dict = {"GSM11111": {"details": "not important", "srx": "SRX/SRX000", "paired": True, "genome": "mm10", "gse":"GSE1000", 
			"aligner": "tophat", "date": date, "report": "useless", "submitter": "ME"}}
		insert_data(args["db"], gsm_dict)
	elif args["p"]:
		print_database(args["db"])
	#	print_schmena(pkg_path)
