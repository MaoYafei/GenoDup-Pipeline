
import sys, getopt
import os

import Bio

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def genefam_dS_Cal(Parameters):
	group_file=Parameters[0]
	pro_file=Parameters[1]
	seq_file=Parameters[2]
	max_genefam=Parameters[3]

###select gene family###

	group_file=open(group_file,'r')
	pro_dic={}
	seq_dic={}

	def famfasta(famlist):
		LED_name=famlist[0][:-1]
		outfile_aa=open(LED_name+'.aa','w')
		outfile_nc=open(LED_name+'.nc','w')
		for i in famlist[1:]:
        		outfile_aa.write('>'+i+'\n')
			outfile_nc.write('>'+i+'\n')
        		outfile_aa.write(pro_dic[i]+'\n')
			outfile_nc.write(seq_dic[i]+'\n')
		outfile_aa.close()
        	outfile_nc.close()
		dic={}

####Filter bad sequence###
        	outfile_del=open(LED_name+'.del','w')
        	for seq_record in SeqIO.parse(LED_name+'.aa','fasta'):
                	seq_seq=str(seq_record.seq)
                	if seq_seq.find('X')!=-1:
                        	seq_seq=seq_seq.replace('X','*')
                        	errfile.write( 'stop codon in string:'+' '+i+'\n')
                	dic[seq_record.id]=str(seq_seq)
        	for seq_record in SeqIO.parse(LED_name+'.nc','fasta'):
                	seq_id=str(seq_record.id)
                	seq_seq=str(seq_record.seq)
                	if seq_seq.find('n')!=-1 or seq_seq.find('K')!=-1 or seq_seq.find('W')!=-1 or seq_seq.find('N')!=-1 or seq_seq.find('M')!=-1:
                        	outfile_del.close()
                        	errfile.write('nnnn: found in'+LED_name+'\n')
                        	break
                	outfile_del.write('>'+seq_id+'\n')
                	Num=len(seq_seq)
                	for start in range(Num):
                        	if dic[seq_id]+'*'==str(Seq(seq_seq,IUPAC.unambiguous_dna).translate()) or dic[seq_id]==str(Seq(seq_seq,IUPAC.unambiguous_dna).translate()):
                                	break
                        	else:
                             		seq_seq=seq_seq[1:]
                	mu=len(seq_seq) % 3
                	if mu !=0:
                        	seq_seq=seq_seq[:-mu]
                	if str(seq_seq)[-3:] in ['TAA','TAG','TGA','taa','tag','tga']:
                        	outfile_del.write(str(seq_seq)[:-3]+'\n')
                	else:
                        	outfile_del.write(str(seq_seq)+'\n')
        	outfile_del.close()
        	return;

	for seq_record in SeqIO.parse(pro_file,'fasta'):
        	pro_dic[seq_record.id]=str(seq_record.seq)

	for seq_record in SeqIO.parse(seq_file,'fasta'):
        	seq_dic[seq_record.id]=str(seq_record.seq)
	print 'Sequence Load Finished'


	errfile=open('errfile.txt','a')
	os.system("rm -r  aln_genefam")
	os.system("mkdir aln_genefam")
	for line in group_file:
    		famlist=line.strip('\n').split(' ')
    		if len(famlist[1:])>=max_genefam:
			continue
		file_name=famlist[0][:-1]
    		famfasta(famlist)
        	print "algiment %s" %(file_name)

###alignment sequence###

        	os.system("mafft --localpair --maxiterate 1000 %s.aa > %s.aa.aln" % (file_name,file_name))
		os.system("perl translatorx_vLocal.pl -i %s.del -a %s.aa.aln -o %s" %(file_name,file_name,file_name))
		os.system('mv %s.nt_* aln_genefam' % (file_name))
        	os.system("rm %s*" %(file_name))
	group_file.close()

	L=os.listdir('aln_genefam')

####creat pairwise sequence for caculate dS###
	os.system('rm -r  pairwise')
	os.system('mkdir pairwise')

	s='pairwise/'


	def get_seq(dic):
        	s=''
        	for key in dic:
                	s=s+'>'+key+'\n'+dic[key]+'\n'
        	return s
	for i in L:
        	print i
        	Ltem=[]
        	for seq_record in SeqIO.parse('aln_genefam/'+i,'fasta'):
                	dic={}
                	dic[str(seq_record.id)]=str(seq_record.seq)
                	Ltem.append(dic)
        	n=0
        	print len(Ltem)
        	for num1 in range(len(Ltem)):
                	for num2 in range(num1+1,len(Ltem)):
                        	outfile=open(s+i+str(n),'w')
                        	n=n+1
                        	outfile.write(get_seq(Ltem[num1])+get_seq(Ltem[num2]))
                        	outfile.close()

###creat script to PAML####
	os.system("rm -r PAML_result")
	os.system("mkdir PAML_result")
	L=os.listdir('pairwise/')
	for i in L:
		file_name=i+'.ctl_cdm'
        	outfile=open(i+'.ctl_cdm','w')
        	outfile.write('seqfile = pairwise/'+i+'\noutfile = '+i+'.cdmrlt\nnoisy = 9\nverbose = 1\nrunmode = -2\nseqtype = 1\nCodonFreq = 2\nmodel = 0\nNSsites = 0\nicode = 0\nfix_kappa = 0\nkappa = 1\nfix_omega = 0\nomega = 0.5')
        	outfile.close()
		os.system("codeml %s" % (file_name))
		os.system("mv *cdmrlt PAML_result")
		os.system("rm %s*" % (i))

