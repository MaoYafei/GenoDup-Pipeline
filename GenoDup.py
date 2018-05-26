
###Checking dependece####
try:
	import sys, getopt
	import os
	import subprocess
	import Bio
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
except Exception:
	print 'You have to install related python modules'
	sys.exit()

if os.system("mkdir test_test_test_test")!=0 or os.system("rm -r test_test_test_test")!=0:
        print 'You need to have authority to make/delete directory'
        sys.exit()

try:
	print 'Mafft:'
	subprocess.check_call(['which', 'mafft'])
except OSError:
	print 'You have to install MAFFT'
	sys.exit()

try:
	print 'Translatorx:'
        subprocess.check_call(['ls', 'translatorx_vLocal.pl'])
except OSError:
        print 'You have to install TRANSLATORX'
        sys.exit()

try:
	print 'Codeml:'
    	subprocess.check_call(['which','codeml'])
except OSError:
        print 'You have to install PAML'
        sys.exit()


####Help Page######

def Help():

	print "This python script will help you to caculate dS values of gene_family pairs or pairwise ohnologs \
in order to detect Whole (large-scale))-genome Duplication.\nMore information in README.file \
               \nParameters:\n\
        -g: Group file name, it contains the gene family cluster and could be produced by OrthoMCL\n \
	-s: Nuclear sequence file name, it contains all the nc squences in your analysis (fasta format)\n \
	-p: Protein sequence file name,	it contains all	the pro squences in your analysis (fasta format)\n \
	-c: Pairwise Ohnologs file name, it contains two Ohnologs in two colums, it could be produced by MCScanX or OrthoMCL or i-ADHoRe\n \
	-n: Maximum number for gene family which you want to analysis, only use with -g \n\
	-c and -g can not be choosed at same time"
	return;


###Reading parameters####
opts, args = getopt.getopt(sys.argv[1:], "hg:n:p:s:c:")
group_file=""
max_genefam=""
pro_file=""
seq_file=""
pair_file=""
for op, value in opts:
	if op == "-g":
    		group_file = value
  	elif op == "-n":
    		max_genefam = float(value)
	elif op == "-p":
		pro_file = value
	elif op == "-s":
		seq_file = value
        elif op == "-c":
                pair_file = value
  	elif op == "-h":
		Help()
    		sys.exit()

if group_file!='' and pair_file!='':
	print '-g and -c can not use together!'
	sys.exit()
if seq_file=='' or pro_file=='':
	print '-s or -p is mandatory parameter(s)'
	sys.exit()
if pair_file=='':
	if max_genefam=='':	
		print 'You have to choose the maximum number of gene famliy you wanna analysis'
		sys.exit()
	else:
		from genefam_dS_Cal import *
		Parameters=[group_file,pro_file,seq_file,max_genefam]
		genefam_dS_Cal(Parameters)
		os.system("grep 'dN/dS= ' PAML_result/led*cdmrlt |cut -f 7 -d '=' > dS_value.txt")
		sys.exit()
if group_file=='':
	from pairwise_dS_Cal import *
	Parameters=[pair_file,pro_file,seq_file]
	pairwise_dS_Cal(Parameters)
	os.system("grep 'dN/dS= ' PAML_result/led*cdmrlt |cut -f 7 -d '=' > dS_value.txt")
	sys.exit()

