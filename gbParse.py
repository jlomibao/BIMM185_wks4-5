#FILE NAME: gbParse.py 
#AUTHOR   : John Francis Lomibao
#PID      : A11591509

import os.path, sys, gzip, re

#Warning in case we don't want to overwrite files in directory
print 'WARNING! The following files will be overwritten:\n'
print ('genome.tbl\nreplicon.tbl\ngene.tbl\nxref.tbl\n'+
	   'exon.tbl\ngeneSynonym.tbl\nfunction.tbl\n')
print 'Is this okay? (RETURN -> yes / CTRL-C -> no)'
raw_input("")

infiles    = sys.argv[1:]
open_files = []

genome_tbl   = 'genome.tbl'
replicon_tbl = 'replicon.tbl'
gene_tbl     = 'gene.tbl'
xref_tbl     = 'xref.tbl'
exon_tbl     = 'exon.tbl'
gene_syn_tbl = 'geneSynonym.tbl'
function_tbl = 'function.tbl'

#variables for each column in genome.tbl
genome_id, genome_name, tax_id, domain              = '', '', '', ''
num_replicons, num_genes_gnm, genome_size, assembly = '', '', '', ''

#unique variables in columns of replicon.tbl
replicon_id, replicon_name, replicon_type, replicon_shape = '', '', '', ''
num_genes_rep, len_replicon, accession, release_date      = '', '', '', ''

#unique variables in columns of gene.tbl
gene_id, locus_tag, protein_id, gene_name = '', '', '', ''
strand, num_exons, len_gene, product      = '', '', '', ''

#unique variables in columns of xref.tbl
xdb, xid = '', ''
xrefs = []

#unique variables in columns of exon.tbl
exon, left_pos, right_pos, len_exon = '', '', '', ''

#unique variable in column of geneSynonym.tbl
gene_syn = []

#unique variable in column of function.tbl
funct = ''

#counters
'''keep count parsing through all genomes'''
cntGenome    = 0
cntReplicons = 0
cntCDS       = 0
'''reset after every genbank file'''
gnmSize      = 0
gnmReplicons = 0
gnmCDS       = 0
'''reset after each replicon'''
repCDS       = 0
'''reset after each gene'''
exonOrder    = 0

#function to empty files
def emptyFile(*filesToEmpty):
	for i in filesToEmpty:
		if os.path.exists(i):
			thisFile = open(i, 'w')
			thisFile.truncate()
			thisFile.close()

#function to write to file			
def writeToFile(fileToWrite, strToPutInFile):
	thisFile = open(fileToWrite, 'a')
	thisFile.write(strToPutInFile)
	strToPutInFile = ''
	thisFile.close()

#open all zipped or unzipped genbank files			
for gbff in infiles:
	if '.gz' in gbff:
		open_files.append(gzip.open(gbff, 'rb'))
	else:
		open_files.append(open(gbff, 'r'))

#empty these files if they exist, so we can write the data we want to them
emptyFile(genome_tbl, replicon_tbl, gene_tbl, xref_tbl,
		  exon_tbl, gene_syn_tbl, function_tbl)
		  
#Using Biopython
from Bio import SeqIO

#Loop through each genbank formatted file, 'gbff'
for gbff in open_files:
	cntGenome += 1
	genome_id  = str(cntGenome)

	#Loop through each replicon in current gbff
	for rec in SeqIO.parse(gbff, 'genbank'):
		cntReplicons  += 1
		gnmReplicons  += 1
		replicon_id    = str(cntReplicons)
		gnmSize       += len(rec.seq)
		len_replicon   = str(len(rec.seq))
		accession      = rec.name
		replicon_shape = rec.annotations["topology"]
		release_date   = rec.annotations["date"]
		replicon_name  = rec.annotations["organism"]
		genome_name    = rec.annotations["organism"]
		
		if (rec.description.find('chromosome') != -1 or
				rec.description.find('genome') != -1):
			replicon_type = 'chromosome'
		elif rec.description.find('plasmid') != -1:
			replicon_type = 'plasmid'

		#get assembly via regular expressions
		m = re.search(r'(Assembly:GCF_)([0-9]+.[0-9])',str(rec.dbxrefs))
		getAssembly = re.search(r'(GCF_[0-9]+.[0-9])', str(m.group()))
		assembly = str(getAssembly.group())
		
		if str(rec.annotations["taxonomy"]).find('Bacteria') != -1:
			domain = 'bacteria'
		elif str(rec.annotations["taxonomy"]).find('Eukaryota') != -1:
			domain = 'eukarya'
		elif str(rec.annotations["taxonomy"]).find('Archaea') != -1:
			domain = 'archaea'
		
		#Loop through each gene in replicon
		for feature in rec.features:
			#check each gene (CDS), then maake sure it is not fragmented
			#'<' indicates partial on 5', '>' indicates partial on 3'
			if (feature.type == "CDS" and str(feature.location).find('<') == -1
									  and str(feature.location).find('>') == -1):
				#increment counters for CDS
				cntCDS += 1
				gnmCDS += 1
				repCDS += 1
				gene_id = str(cntCDS)
				
				cntExon = 0
				cntExon += str(feature.location).count('+')
				cntExon += str(feature.location).count('-')
				num_exons = str(cntExon)
				exonList = re.findall(r"([0-9]+):([0-9]+)", str(feature.location))
				gen_length_ctnr = 0
				for i in exonList:
					#increment exonOrder
					exonOrder += 1
					exon      = str(exonOrder)
					left_pos  = i[0]
					right_pos = i[1]
					len_exon  = str(abs(int(right_pos)-int(left_pos)))
					gen_length_ctnr += int(len_exon)
					writeData = (gene_id+'\t'+exon+'\t'+left_pos+'\t'+right_pos
								 +'\t'+len_exon+'\n')
					writeToFile(exon_tbl, writeData)
				len_gene = str(gen_length_ctnr)
				exonOrder = 0
				try:
					#regex to find whether the strand is fwd or rev, and set that strand to that value
					m = re.search(r"([\-\+()]+)", str(feature.location))
					
					if m.group() == '(+)':
						strand = 'F'
					elif m.group() == '(-)':
						strand = 'R'
					
					try:
						gene_name = feature.qualifiers['gene'][0]
					except KeyError:
						gene_name = feature.qualifiers['locus_tag'][0]
					locus_tag = feature.qualifiers['locus_tag'][0]
					product   = feature.qualifiers['product'][0]
					
					#for protein id, we only need the data before the 'dot'
					m = re.search(r"([A-Z]+_[0-9]+)", feature.qualifiers['protein_id'][0])
					protein_id = m.group(1)
					
					#if there are any, gene_syn will contain all of that gene's synonyms
					gene_syn = feature.qualifiers['gene_synonym'][0]
					if gene_syn != '':	
						if gene_syn.find(';') != -1:
							synList = re.findall(r"([a-zA-Z0-9]+)", gene_syn)
							for j in range(len(synList)):
								writeData = (gene_id+'\t'+synList[j]+'\n')
								writeToFile(gene_syn_tbl, writeData)
						else:
							writeData = (gene_id+'\t'+gene_syn+'\n')
							writeToFile(gene_syn_tbl, writeData)
			
					funct = feature.qualifiers['function'][0]
					writeData = (gene_id+'\t'+funct+'\n')
					writeToFile(function_tbl, writeData)
					
					xrefs = feature.qualifiers["db_xref"]
				except KeyError:
					pass
				#xref table should start with reference sequence at the start of each gene
				#refer to it with protein_id
				xdb = 'refseq'
				xid = protein_id
					
				writeData = (gene_id+'\t'+xdb+'\t'+xid+'\n')
				writeToFile(xref_tbl, writeData)
				for i in xrefs:
					if i.find('GI:') != -1:
						m = re.match(r"(GI):([0-9]+)", str(i))
						xdb = 'gi'
						xid = str(m.group(2))
					if i.find('ASAP:') != -1:
						m = re.match(r"(ASAP):(ABE-[0-9]+)", str(i))
						xdb = 'asap'
						xid = str(m.group(2))
					if i.find('UniProtKB/Swiss-Prot:') != -1:
						m = re.match(r"([A-Za-z]+\/[A-Za-z-]+):([0-9A-Z]+)", str(i))
						xdb = 'uniprot'
						xid = str(m.group(2))
					if i.find('EcoGene:') != -1:
						m = re.match(r"(EcoGene):([A-Z]+[0-9]+)", str(i))
						xdb = 'ecogene'
						xid = str(m.group(2))
					if i.find('GeneID:') != -1:
						m = re.match(r"(GeneID):([0-9]+)", str(i))
						xdb = 'geneid'
						xid = str(m.group(2))
					writeData = (gene_id+'\t'+xdb+'\t'+xid+'\n')
					writeToFile(xref_tbl, writeData)
				
				
				writeData = (gene_id+'\t'+replicon_id+'\t'+genome_id+'\t'+locus_tag
							 +'\t'+protein_id+'\t'+gene_name+'\t'+strand+'\t'+
							 num_exons+'\t'+len_gene+'\t'+product+'\n')
				writeToFile(gene_tbl, writeData)
					
			elif feature.type == "source":
				m = re.search(r'(taxon:)([0-9]+)',str(feature.qualifiers['db_xref'][0]))
				getID = re.search(r'([0-9]+)', str(m.group()))
				tax_id = str(getID.group())
			
		num_genes_rep = str(repCDS)
		writeData = (replicon_id+'\t'+genome_id+'\t'+replicon_name+'\t'+
					 replicon_type+'\t'+replicon_shape+'\t'+num_genes_rep+'\t'
					 +len_replicon+'\t'+accession+'\t'+release_date+'\n')
		writeToFile(replicon_tbl, writeData)
		repCDS = 0
		
	num_replicons = str(gnmReplicons)
	num_genes_gnm = str(gnmCDS)
	writeData = (genome_id+'\t'+genome_name+'\t'+tax_id+'\t'+domain+'\t'+
		 num_replicons+'\t'+num_genes_gnm+'\t'+str(gnmSize)+'\t'+assembly+'\n')
	writeToFile(genome_tbl, writeData)
	gnmSize      = 0
	gnmReplicons = 0
	gnmCDS       = 0
	
'''
try:
	for i in open_files:
		print i
except IOError:
	pass
'''
