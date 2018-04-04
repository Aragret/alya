#!/home/aragret/anaconda3/bin/python

import os
import re

path = '/home/aragret/Alina/other_projects/alya/1_CDS_by_species'
# path = '/home/aragret/Alina/other_projects/alya/noverlaps'

with open('genes.txt', 'w') as outfile:
	for file in os.listdir(path):
	    with open('1_CDS_by_species/{0}'.format(file)) as infile:
	    	for line in infile:
		    	if line[0] == '>':
		    		tab_split = line.split('\t')
		    		if tab_split[2] == 'translation':
		    			continue
		    		species = tab_split[1]
		    		gene = tab_split[0][1:]
		    	if re.match('[ATGC]+', line) != None:
		    		# print(line)
		    		outfile.write('{0}\t{1}\t{2}'.format(species, gene,line))
