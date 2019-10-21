'''
Author: Yarrow Madrona
A script to automate the process of designing primers for protein enjineering libraries using
the NEB Basechanger website(https://nebasechanger.neb.com/).
Dependancies: Python 3 including pandas and Selenium modules
This script uses Selenium, usually used for website testing as a work-around to automate high-throughput primer design
from the NEB basechanger website. This was necessary as the website is written with Java and includes many pop-up windows
The downside to this is that it is slower (though still faster than doing by hand) and includes some absolute xpaths.
Not having time to recode everything in relative xpaths, the code will break if there is any change to xpaths used here.
If this is the case, please let me know and I will fix it. Of course it would be better to reimplement the algorithm used
by NEB, they were not interested in giving me this information (big surprise). And researchers tend to want to stick to
what works, thus the implementation of Selenium on a Java driven website.

Inputs: 1. DNA fasta file that encodes the protein to be mutated.
		2. List of Amino acid changes using three letter codes (Sorry, I havn't added dictionary for one letter codes)
		The list should be a two column CSV with no headers as follows:
		18, Gln
		221, Tyr
		331, Arg
		
The first column contains the amino acid position to be mutated. The second position contains the three letter code you
want to mutate to. Again, sorry you can't use one letter codes.
The nucleotides to be changed to are codon optimized for expression in E.coli. I can add this as an alternative input if
there is any interest. The amino acids to change each position to are also hard coded as follows:
Leu', 'Val', 'Ile', 'Met', 'Phe', 'Trp', 'Ser', 'Thr', 'Tyr', 'Asn', 'Gln', 'Asp', 'Glu', 'Lys', 'Arg', 'His'
This can also be an optional input if anyone is interested.
The output contains temporary files stored in a newly created "Temp" directory. This is to ensure that we get an output in
case of getting booted off the website, though this has not been a problem so far. Final output is a list of all primers as well
as two seperate lists for forward and reverse primers. The primers are listed by the amino acid position that is being changed followed
by the amino acid it is changed to.
'''

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.keys import Keys
import time
import sys
import csv
import pandas as pd
import os
import argparse


def get_data(aa_csv, fasta, offset):
	os.mkdir('temp', exist_ok=True)
	aa_pos_list = []
	aa_list = []
	with open(aa_csv, 'r') as csvfile:
		csv_reader = csv.reader(csvfile, delimiter=',')
		for row in csv_reader:
			aa_pos_list.append(int(row[0])+offset)  #Use 92 for POL6
			aa_list.append(row[1])
	print ("PDB positions + offset of %s for construct positions:")%(offset)
	for pos in  aa_pos_list:
		print (pos)
	with open(fasta) as f:
		fasta = f.read()
	nt_start = [int(i)*3-2 for i in aa_pos_list]
	nt_stop = [int(i)*3 for i in aa_pos_list]
	nt_start_stop = zip(nt_start, nt_stop)
	print ('nucleotide positions:')
	for start, stop in (zip(nt_start, nt_stop)):
		print (start, stop)
	library_aa = ['Leu', 'Val', 'Ile', 'Met', 'Phe', 'Trp', 'Ser', 'Thr', 'Tyr', 'Asn', 'Gln', 'Asp', 'Glu',\
				  'Lys', 'Arg', 'His']
	codon_table = {'Gly':'GGC', 'Ala':'GCG', 'Leu':'CTG', 'Glu':'GAA', 'Gln':'CAG', 'Lys':'AAA', 'Phe':'TTT',\
				   'Tyr':'TAT', 'Arg':'CGT', 'Trp':'TGG', 'Thr':'ACC', 'Met':'ATG', 'Pro':'CCG', 'Val':'GTG',\
				   'Ile':'ATC', 'Ser':'AGC', 'His':'CAC', 'Asp':'GAT', 'Asn':'AAC', 'Cys':'TGC'}
	return nt_start_stop, fasta, library_aa, codon_table
	
#Get primers from NEB Basechanger website using Selenium
def get_primers(nt_start_stop, fasta, library_aa, codon_table):
	driver = webdriver.Chrome()
	driver.get('http://nebasechanger.neb.com')
	main_page = driver.current_window_handle
	time.sleep(2)
	add_seq_button = driver.find_element_by_id('addseqbutton')
	add_seq_button.click()
	time.sleep(2)
	add_seq_window = driver.current_window_handle
	driver.switch_to.window(add_seq_window)
	time.sleep(2)
	driver.find_element_by_id("seqid").send_keys("1743 pol6")
	time.sleep(2)
	driver.find_element_by_id("seqfield").send_keys(fasta)
	time.sleep(2)
	driver.find_element_by_xpath("/html/body/div[10]/div[10]/div/button[3]").click()
	time.sleep(2)
	driver.switch_to.window(main_page)
	time.sleep(2)
	driver.find_element_by_xpath("/html/body/div[3]/div[1]/div[2]/div[3]/div[2]/label[1]").click()
	time.sleep(2)
	
	#Now start the loop to get all ppsitions with 16 amino acids
	primer_dict = {}
	for start_pos, end_pos in nt_start_stop:
		time.sleep(2)
		nt_start = driver.find_element_by_id("selStart")
		nt_start.send_keys(str(start_pos))
		nt_end = driver.find_element_by_id("selEnd")
		nt_end.send_keys(str(end_pos))			
		for aa in library_aa:
			#time.sleep(2)
			nt_to_change_to = driver.find_element_by_id('repseq')
			nt_to_change_to.send_keys(codon_table[aa])
			time.sleep(1)
			f_oligo = driver.find_element_by_xpath("//table/tbody/tr[2]/td[2]").text
			r_oligo = driver.find_element_by_xpath("//table/tbody/tr[3]/td[2]").text
			f_anneal_temp = driver.find_element_by_xpath("//table/tbody/tr[2]/td[5]").text
			r_anneal_temp = driver.find_element_by_xpath("//table/tbody/tr[3]/td[5]").text
			rec_anneal_temp = driver.find_element_by_xpath("//table/tbody/tr[2]/td[6]").text
			#time.sleep(2)
			nt_to_change_to.clear()
			#time.sleep(2)
			aa_pos = end_pos/3
			aa_name_f = str(int(end_pos/3-92))+"_"+aa+"_fwd"
			aa_name_r = str(int(end_pos/3-92))+"_rev"
			primer_dict[aa_name_f] = [f_oligo, rec_anneal_temp]
			primer_dict[aa_name_r] = [r_oligo, rec_anneal_temp]
			#Write out a temp file every 100 primers
			if (len(primer_dict) % 100 )== 0:
				df = pd.DataFrame.from_dict(primer_dict, orient='index', columns=['Sequence', 'rec anneal temp'])
				df.drop_duplicates(subset='Sequence', keep='first', inplace = True)				
				df.to_csv('temp/'+aa_name_f+'_tmp.csv')
			
			else:
				continue
		nt_start.clear()
		nt_end.clear()
		
	primer_df = pd.DataFrame.from_dict(primer_dict, orient='index', columns=['Sequence', 'rec anneal temp'])
	primer_df.index.names = ['Primer name']
	primer_df = primer_df.reset_index()
	primer_df.drop_duplicates(subset='Sequence', keep='first', inplace = True)
	primer_df.to_csv('test-NEB_basechanger_primer_list.csv')
	return primer_df, library_aa
	
	#Make input for IDT Forward primers
def make_idt_input(primer_df, library_aa):
	reverse_primer_df = primer_df[primer_df['Primer name'].str.contains('_rev')]
	A = ['A1', 'A2','A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12']
	B = [x.replace('A', 'B') for x in A]
	C = [x.replace('A', 'C') for x in A]
	D = [x.replace('A', 'D') for x in A]
	E = [x.replace('A', 'E') for x in A]
	F = [x.replace('A', 'F') for x in A]
	G = [x.replace('A', 'G') for x in A]
	H = [x.replace('A', 'H') for x in A]
	rev_primer_labels = A+B+C+D+E+F+G+H
	rev_primer_labels = rev_primer_labels * (int(len(reverse_primer_df) / 96)) + rev_primer_labels[0:(len(reverse_primer_df) % 96)]
	reverse_primer_df.insert(loc=0, column='Plate Pos', value=rev_primer_labels)
		
	reverse_primer_df.to_csv('IDT-input-primers-rev.csv')
	#Make input for IDT reverse primers
	forward_primer_df = primer_df[primer_df['Primer name'].str.contains('_fwd')]
	forward_primer_df_set1 = forward_primer_df[forward_primer_df['Primer name'].apply(lambda x: x[x.find('_')+1:x.find('_f')]).isin(library_aa[0:8])]
	forward_primer_df_set2 = forward_primer_df[forward_primer_df['Primer name'].apply(lambda x: x[x.find('_')+1:x.find('_f')]).isin(library_aa[8:])]
	print(library_aa[0:9])
	print(library_aa[9:])
	print(forward_primer_df_set2.head())
	col1 = ['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1']
	col2 = [x.replace('1', '2') for x in col1]
	col3 = [x.replace('1', '3') for x in col1]
	col4 = [x.replace('1', '4') for x in col1]
	col5 = [x.replace('1', '5') for x in col1]
	col6 = [x.replace('1', '6') for x in col1]
	col7 = [x.replace('1', '7') for x in col1]
	col8 = [x.replace('1', '8') for x in col1]
	col9 = [x.replace('1', '9') for x in col1]
	col10 = [x.replace('1', '10') for x in col1]
	col11 = [x.replace('1', '11') for x in col1]
	col12 = [x.replace('1', '12') for x in col1]
	fwd_primer_labels = col1+col2+col3+col4+col5+col6+col7+col8+col9+col10+col11+col12
	fwd_primer_labels_set1 = fwd_primer_labels * (int(len(forward_primer_df_set1) / 96)) + fwd_primer_labels[0:(len(forward_primer_df_set1) % 96)]
	forward_primer_df_set1.insert(loc=0, column='Plate Pos', value=fwd_primer_labels_set1)
	fwd_primer_labels_set2 = fwd_primer_labels * (int(len(forward_primer_df_set2) / 96)) + fwd_primer_labels[0:(len(forward_primer_df_set2) % 96)]
	forward_primer_df_set2.insert(loc=0, column='Plate Pos', value=fwd_primer_labels_set2)
	forward_primer_df_set1.to_csv('IDT-input-primers-fwd-set1.csv')
	forward_primer_df_set2.to_csv('IDT-input-primers-fwd-set2.csv')

if __name__ == "__main__":
	# --COMMAND LINE ARGUMENTS
	parser = argparse.ArgumentParser(description='Automating Baschanger primer design')
	parser.add_argument('--fasta', '-f', required=True, help='fasta DNA file you want to mutate' )
	parser.add_argmuent('--aa_csv', '-aa', required=True, help='Amino acid CSV containing amino acids to mutate')
	parser.add_argument('--offset', '-o', default=0, required=False, help='offset in amino acids between bases in DNA file and protein')
	
	args = parser.parse_args()
	aa_csv = args.aa_csv()
	fasta = args.fasta()
	offset = int(args.offset())
	
	#--CALL FUNCTIONS
	nt_start_stop, fasta, library_aa, codon_table, library_aa = get_data(aa_csv, fasta, offset)
	primer_df, library_aa = get_primers(nt_start_stop, fasta, library_aa, codon_table, library_aa)
	make_idt_input(primer_df, library_aa)