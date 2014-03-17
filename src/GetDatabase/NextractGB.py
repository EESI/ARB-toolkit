## Created:	E. Reichenberger
## Date:	10.30.2011

'''
Purpose: To extact information from GenBank file and place the information to a new csv file.

Steps:
1. Create text file (write mode)
2. Loop through directory to make list of all pertinant files 
3. Loop through all *.gb files
        i. Open GenBank files (read mode)
        ii. extract pertinent info from each SEQ record
		a. Use taxID to generate taxonomic assignment (separate module) avoids naming inconsistency found in genbank files
        iii. concatenate info to end of txt file (tab delimated form)
	iv. concatenate sequence of SEQ record to end of text file

Special Notes:
1. All variables should be initialized to 'NA'. 
	Otherwise, if a variable does not exist in a genbank file, it will be unassigned. Unassigned variables will not import into ARB properly causing fatal errors.
2. Variables should be stripped of all non-alphanumeric characters with the exception of the underscore '_' characters.
	All non-alphanumeric characters with the exception of '_' will cause fatal errors when imported into ARB.

Additional Notes:
GI: series of digits that are assigned consecutively to each sequence record processed by NCBI. The GI number bears 
no resemblance to the Accession number of the sequence record. Found in: 
	GI number found in VERSION field of the database record (nucleotide sequence)
	GI number found in CDS/db_xref field of a nucleotide database record, and the VERSION field of a protein database record

Region (associated with Accession):

References: A genbank may have multiple references which are all stored in a reference list [reference1, reference2, ... referenceN].
The first entry in the Reference list will have an index of zero. The last entry (N) in the Reference list will have an index of N-1.
Reference lists can be accessed as follows: 
record.annotations['references'][index]

TaxonID and taxonomy: Genbank files are often inconsistent with other genbank files. An example might be that rather than having a decending taxonomic hierarchy
(...Phylum, Class, Order...) it may have been entered in ascending order (...Order, Class, Phylum...). To work around this issue, the taxonID is extracted 
from genbank files. The taxonID is used in a separate module to provide the scientic name and it's corresponding rank (using tree files downloaded from NCBI).
However, recent entries may not have made their way into the tree files and will not exist there causing 'None' to be return for the name and the rank.
If either the rank or name returns 'None', the 'Taxonomic Assignment Loop' will not be entered (if spec != None and spec_rank != None:), and the 
taxonomic names will retain their 'NA' value. To ensure that some taxonomic information is extracted from the file, the taxonomy listed in the file is extracted. 

CDS: In our case, the cds has already been converted into amino acids and stored in a separate directory. The cds starts on the second line. A list of fasta files will be
created in the same order as the genbank files. The cds from each fasta file will be extracted and stored in a cdsMetaData.txt file.

SPECIAL NOTE: This script has been coupled with other scripts as well as file 93279115_MFS.gb. This file was originally downloaded directly from NCBI, but this script would 
crash when it came into contact with this file. The reason is that this file illegally combines join and bond together. According to https://redmine.open-bio.org/issues/3197,
"location operator "complement" can be used in combination with either "join" or "order" within the same location; combinations of "join" and "order" within the same location (nested operators) are illegal." The accompanying file has had all instances of join(bond()...) removed. If you overwrite this file, please note that this script will not function properly and you most remove all instances of join(bond(number), bond(number)) with join(number,number).
'''

import Bio
import os
import glob
import sys
import re #regular expressions
import csv #comma separated values format
import inspect
import collections
from Bio import GenBank
from Bio import SeqIO
from Bio import SeqFeature
from Bio.GenBank import Record
from Bio.GenBank.Record import Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqFeature import Reference, SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from StringIO import StringIO

lineStarter = '>' #Arb requires all non-sequence lines to begin with '>', this value will never change 
 
'''
########################TREE/NODE INFO###############################
# Not beholden to relying on naming inconsistency in GenBank files, use build_tax_tree() to extract taxonomy from NCBI Taxonomy.
# At least 6GB of RAM Required for this function
current = os.getcwd()
os.chdir('NCBI_Taxonomy')
from taxid2lineageModule import build_tax_tree #taxonomic assignment
from taxid2lineageModule import get_taxonomy #taxonomic assignment
node,Nname = build_tax_tree() 
getTaxonomy = 1
os.chdir(current)
'''

databases = ['MFS'] #this is list of directory(ies) containing gb files

for database in databases:
	########################Checking for and removing Old Files###############################
	#If this script is ran more than once, metaData & CDSmetaData file are appended not overwritten. Following lines delete previous files.
	# IMPORTANT, 'Output' must be a subdirectory of the database directory!!! E.g. 'mkdir Output' 
        if os.path.exists(database + '/Output/' + database + '_CDSmetaData.txt'):
                os.remove(database + '/Output/' + database + '_CDSmetaData.txt')
        if os.path.exists(database + '/Output/' + database + '_ExonmetaData.txt'):
                os.remove(database + '/Output/' + database + '_ExonmetaData.txt')                                                                 
        if os.path.exists(database + '/Output/' + database + '_metaData.txt'):
                os.remove(database + '/Output/' + database + '_metaData.txt')                                                                 
        if os.path.exists(database + '/Output/' + database + '_CDSsubmetaData.txt'):
                os.remove(database + '/Output/' + database + '_CDSsubmetaData.txt')                                                                 


	########################FILE INFORMATION###############################
	OutputMetaData = open(database + '/Output/' + database + '_metaData.txt','a') # Create a file object in append mode for metaData
	OutputExonData= open(database + '/Output/' + database + '_ExonmetaData.txt', 'a') #create file object ina ppend mode for CDS
	OutputCDSData= open(database + '/Output/' + database + '_CDSmetaData.txt', 'a') #create file object ina ppend mode for CDS
	OutputCDSsubfeatureData = open(database + '/Output/' + database + '_CDSsubmetaData.txt', 'a') #create file object ina ppend mode for CDSsubfeatures
	headerLabelFile = open(database + '/Output/' + database + '_Field_labels.txt', 'w')
	      	
	path = database + '/'

	fileList = [] #list for gb files in folder

        for infile in glob.glob( os.path.join(path, '*.gb') ): #file extension type
                fileList.append(infile)

	fileCounter = 0 #will increment +1 for each loop & will be the unique number
	uniqueID = str(fileCounter) #each record must have a unique ID

	######################Looping Information#############################
	for file in fileList: 
		uniqueID = str(fileCounter) #each record must have a unique ID
		stripFilename = file.replace(database + '/', '')
		#print stripFilename + '\t' + uniqueID #This can be uncommented to see which file is being analyzed
		record = SeqIO.parse(file, 'genbank').next() #use this to prevent exception errors for multiple handles.
		
		metaData = {'locus' : 'NA', 'definition' : 'NA', 'accession' : 'NA', 'accession_region' : 'NA', 'version' : 'NA', 'gi' : 'NA', 'Source' : 'NA', 'Organism' : 'NA', 'taxonomy' : 'NA', 'size' : len(record), 'date' : 'NA', 'authors' : 'NA', 'title' : 'NA', 'journal' : 'NA', 'pubmed' : 'NA', 'comments' : 'NA', 'source_mol_type' : 'NA', 'source_strain' : 'NA', 'source_variety' : 'NA', 'taxonID' : 'NA', 'source_note' : 'NA', 'source_chromosome' : 'NA', 'mRNA_transcript_id' : 'NA', 'mRNA_product' : 'NA', 'cds_instructions' : 'NA', 'cds_subfeature_instructions' : 'NA', 'CDS_locus_tag' : 'NA', 'CDS_protein_id' : 'NA', 'CDS_translation' : 'NA', 'CDS_transl_table' : 'NA', 'CDS_product' : 'NA', 'CDS_note' : 'NA', 'strand' : 'NA', 'superkingdom' : 'NA', 'kingdom' : 'NA', 'subkingdom' : 'NA', 'phylum' : 'NA', 'Class' : 'NA', 'order' : 'NA', 'family' : 'NA', 'genus' : 'NA', 'species' : 'NA', 'strain' : 'NA', 'subspecies' : 'NA', 'subclass' : 'NA', 'subphylum' : 'NA', 'subgenus' : 'NA', 'suborder' : 'NA', 'superclass' : 'NA', 'varietas' : 'NA', 'superfamily' : 'NA', 'subfamily' : 'NA', 'superphylum' : 'NA', 'subtribe' : 'NA', 'candidate' : 'NA', 'cultivar' : 'NA', 'domain' : 'NA', 'fileName': 'NA', 'exonInstructions' : 'NA'}

		#stripFilename = file.replace(database + '/', '')
		#print stripFilename
		metaData['fileName'] = stripFilename
		seq = 'NA' #record sequence	
		cds_ranges = [] #this is a list that will contain all exons ranges
		CDS_sequence = 'NA' #this is the concatenated CDS sequence
		CDS_subfeaturesequence = 'NA' #this is the concatenated CDS sequence
		cds_instructionList = []
		cds_subfeatureInstructionList = []
		strandList = []
		CDS_sequenceList = []
		CDS_subfeatureSequenceList = []
		starts = []
		stops = []
		strand = []
		exon = 'NA'
		exonInstructions = []
		exon_sequenceList = []
	######################################################################################################################
	#In top portion of genbank file, looking for Locus, Definition, Version, and Sequences. 
	#For a list of additional items: print dir(record)
	#To get values, use record.attribute_name
	######################################################################################################################

		for attributes in dir(record):#loops through all attributes in dir(record)
			attributes = attributes.lower() #converts to lower case to avoid potential differences in file styles
			if attributes == 'name':
				metaData['locus'] = record.name 
			if attributes == 'description':
				metaData['definition'] = record.description
			if attributes == 'id':
				metaData['version'] = record.id

	######################################################################################################################
	#In top portion of genbank file, looking for items Date, Source, Organism, Taxonomy, Reference (author, title, journal,
	# pubmedID) & comments. 
	#For a list of additional items: print dir(record.annotations) or record.annotations.keys()
	#To get values, use record.annotations[key_name]
	######################################################################################################################
		for key in record.annotations.keys(): #loops through all annotation keys
			key = key.lower() #converts to lower case to avoid potential differences in file styles
			if key == 'date':
				metaData['date'] = record.annotations['date']
			if key == 'source':
				metaData['Source'] = record.annotations['source']
			if key == 'organism':
				metaData['Organism'] = record.annotations['organism']
			if key == 'comment':
				metaData['comments'] = record.annotations['comment']
			if key == 'taxonomy':
				metaData['taxonomy'] = record.annotations['taxonomy']
			if key == 'gi':
				metaData['gi'] = record.annotations['gi']
			if key == 'references':
				#There may be more than reference set. Sets can be accessed with indices. The first set will have an index of 0. 
				for ref in dir(record.annotations['references'][0]): #Extracting info from 1st reference only 
					ref = ref.lower() #converts to lower case to avoid potential differences in file styles
					if ref == 'authors':
						metaData['authors'] = record.annotations['references'][0].authors #Authors 
					if ref == 'title':
						metaData['title'] =  record.annotations['references'][0].title #Title
					if ref == 'journal':
						metaData['journal'] = record.annotations['references'][0].journal #Journal
					if ref == 'pubmed_id':
						metaData['pubmed'] =  record.annotations['references'][0].pubmed_id #pubmed id
			if key == 'accessions': #accessions can also be a list and have more than one entry
				accessionLength = len(record.annotations['accessions'])	
				if accessionLength == 1: #first entry in list is actual number
					metaData['accession'] = record.annotations['accessions']
				if accessionLength > 1:
					i = 0
					for acc in record.annotations['accessions']:
						acc = acc.upper() #changes to upper case to avoid potential differences in files styles
						metaData['accession'] = record.annotations['accessions'][0] #first entry in list [0] is actual number
						if acc == 'REGION:': #this will not always exist
							metaData['accession_region'] = record.annotations['accessions'][i + 1] #if REGION exists @ index N, index N+1 will hold the REGION value. 
						i+=1

	######################################################################################################################
	#Looking for Feature items (source, mRNA, CDS)
	#source (mol_type, strain, variety, db_xref (taxon), note, chromosome, ) 
	#print source.qualifiers to see what is here once source has been defined
	#For a list of additional items: print record.features
	#To get values, print record.features[featureIndexNumber].qualifiers[qualifierName]
	######################################################################################################################
		featureCount = 0	
		for feat in record.features:
			if feat.type == 'source':
				source = record.features[featureCount] #local variable to hold position
				for qualifiers in source.qualifiers:
					if qualifiers == 'mol_type':
						metaData['source_mol_type'] = source.qualifiers['mol_type']
					if qualifiers == 'strain':
						metaData['source_strain'] = source.qualifiers['strain']
					if qualifiers == 'variety':
						metaData['source_variety'] = source.qualifiers['variety']
					if qualifiers == 'db_xref':
						taxCount = 0
						for taxons in source.qualifiers['db_xref']:
							taxons = str(taxons).lower()
							if taxons.startswith('taxon:') ==  True:
								metaData['taxonID'] = source.qualifiers['db_xref'][taxCount]
								metaData['taxonID'] = metaData['taxonID'].replace('taxon:','')
								metaData['taxonID'] = re.sub(r'[^\w ]', '', metaData['taxonID'])
							taxCount = taxCount + 1
					if qualifiers == 'note':
						metaData['source_note'] = source.qualifiers['note']
						metaData['source_note'] = str(metaData['source_note'])
						metaData['source_note'] = metaData['taxonID'].replace('note:','')
					if qualifiers == 'chromosome':
						metaData['source_chromosome'] = source.qualifiers['chromosome']
						metaData['source_chromosome'] = str(metaData['source_chromosome'])

	######################################################################################################################
	#mRNA (product, transcript_id)
	#print mRNA.qualifiers use this to see what is there once mRNA has been initialized
	######################################################################################################################

			if feat.type == 'mRNA':
				mRNA = record.features[featureCount] #local variable, place holder
				for mRNA_qualifiers in mRNA.qualifiers:
					if mRNA_qualifiers == 'product':
						metaData['mRNA_product'] = mRNA.qualifiers['product']
					if mRNA_qualifiers == 'transcript_id':
						metaData['mRNA_transcript_id'] = mRNA.qualifiers['transcript_id']

	######################################################################################################################
	#Exon (location, , protein_id, translation, translation table, product, note, instructions)
	# from http://biopython.org/DIST/docs/tutorial/Tutorial.html: "this gene/CDS has location string 4343..4780, or in 
	# Python counting 4342:4780."
	######################################################################################################################

		#print record.features  
			if feat.type == 'exon':
				exonLocation = str(feat.location)
				exonLocation = re.sub(r'[\'\]\[]', '', exonLocation)
				exons = exonLocation
				exonLocation = exonLocation.replace('>', 'gt')
				exonLocation = exonLocation.replace('<', 'lt')
				exonInstructions.append(exonLocation)
				exons = exons.replace('>','')
				exons = exons.replace('<','')
				splitLocation = exons.rpartition(':')
				exonStart = int(splitLocation[0])
				if exonStart != 0:
					exonStart -= 1 #python counts start at 0
				exonStop = int(splitLocation[2])
				exon = record.seq[exonStart:exonStop]
				if feat.strand == -1:
					exon = exon.complement()
				exon_sequenceList.append(str(exon)) 

				#exon_sequenceList.append(str(exon)) 

			metaData['exonInstructions'] = str(exonInstructions) 
			metaData['exonInstructions'] = re.sub(r'[ \'\]\[]', '', metaData['exonInstructions'])

	######################################################################################################################
	#CDS (locus_tag, protein_id, translation, translation table, product, note, instructions)
	#print CDS.qualifiers or print CDS to see what is there once CDS has been initialized
	######################################################################################################################

			if feat.type == 'CDS':
				CDS = record.features[featureCount]
				if fileCounter == 0:
					print CDS
				for CDS_qualifiers in CDS.qualifiers:
					if CDS_qualifiers == 'locus_tag':
						metaData['CDS_locus_tag'] = CDS.qualifiers['locus_tag']
					if CDS_qualifiers == 'protein_id':
						metaData['CDS_protein_id'] = CDS.qualifiers['protein_id']
					if CDS_qualifiers == 'translation':
						metaData['CDS_translation'] = CDS.qualifiers['translation']
					if CDS_qualifiers == 'transl_table':
						metaData['CDS_transl_table'] = CDS.qualifiers['transl_table']
					if CDS_qualifiers == 'product':
						metaData['CDS_product'] = CDS.qualifiers['product']
					if CDS_qualifiers == 'note':
						metaData['CDS_note'] = CDS.qualifiers['note']
				if len(CDS.sub_features) != 0:
					for subfeatures in CDS.sub_features:
						cds_subfeatureLocation = str(subfeatures.location)
						cds_subfeatureLocation = re.sub(r'[ \'\]\[]', '', cds_subfeatureLocation)
						cds_subfeatureLocationStr = cds_subfeatureLocation
						cds_SubFeatureLocation = cds_subfeatureLocation.replace('>','gt')
						cds_SubFeatureLocation = cds_SubFeatureLocation.replace('<','lt')
						cds_subfeatureInstructionList.append(cds_SubFeatureLocation)
						cds_subfeatureLocationStr = cds_subfeatureLocationStr.replace('>','')
						cds_subfeatureLocationStr = cds_subfeatureLocationStr.replace('<','')
						splitcdssubfeature = cds_subfeatureLocationStr.rpartition(':')
						subfeaturestart = int(splitcdssubfeature[0])
						if subfeaturestart != 0:
							subfeaturestart -= 1 #python counts start at 0, ref files start at 1
						subfeaturestop = int(splitcdssubfeature[2])
						subfeatureSequence = record.seq[subfeaturestart:subfeaturestop]
						if subfeatures.strand == -1:
							subfeatureSequence = subfeatureSequence.complement()
						sub_featureSequence = str(subfeatureSequence.tostring()) #str should not shorten sequence post vrs 1.45
						CDS_subfeatureSequenceList.append(sub_featureSequence)
							
				cds = CDS.location #local variable, place-holder
				cdsStr = str(cds)
				cdsStr = re.sub(r'[ \'\]\[]', '', cdsStr)					
				instructions = cdsStr.replace('>','gt')
				instructions = instructions.replace('<','lt')
				cds_instructionList.append(instructions)
				cdsStr = cdsStr.replace('>','')
				cdsStr = cdsStr.replace('<','')
				splitcds = cdsStr.rpartition(':')
				start = int(splitcds[0])
				if start != 0:
					start -= 1 #python counts start at 0, ref files start at 1
				stop = int(splitcds[2])
				cdsSequence = record.seq[start:stop]
				if CDS.strand == -1:
					metaData['strand'] = 'complement'
					#strand.append('complement')
					cdsSequence = cdsSequence.complement()
					#CDS_sequenceList.append(str(cdsSequence.complement()))
				else:
					metaData['strand'] = 'forward'
					#strand.append('forward')
					#CDS_sequenceList.append(str(cdsSequence)) #record.seq would also work if sequence is CDS 
				if not cds_instructionList:
					cds_instructionList.append('NA')
				metaData['cds_instructions'] = str(cds_instructionList)
				metaData['cds_instructions'] = re.sub(r'[ \'\]\[]', '', metaData['cds_instructions'])
				CDS_sequenceList.append(str(cdsSequence)) #record.seq would also work if sequence is CDS 

			featureCount+=1

		if not cds_subfeatureInstructionList:
			cds_subfeatureInstructionList.append('NA')
		metaData['cds_subfeature_instructions'] = str(cds_subfeatureInstructionList)
		metaData['cds_subfeature_instructions'] = re.sub(r'[ \'\]\[]', '', metaData['cds_subfeature_instructions'])

	######################TAXONOMIC ASSIGNMENT LOOP#############################
	# This uses an id to taxon module to get name and ranking assignment

		if 'getTaxonomy' in locals():
			spec,spec_rank = get_taxonomy(metaData['taxonID'],node,Nname) #calls module
			taxonName = ['superkingdom','kingdom','subkingdom','superphylum','phylum','subphylum','superclass','subclass','order','suborder','superfamily','family','subfamily','genus','subgenus','species','subspecies','strain','varietas','candidate','subtribe','cultivar','domain']
			rankCount = 0
			
			if spec != None and spec_rank != None: #Enter loop only if name and rank return actual values, otherwise, loop is skipped and variables retain 'NA' assignment. 
				for ranking in spec_rank:
					ranking = ranking.lower() #change to lower case in order to avoid potential differences in case.
	
					if ranking == 'class':
						metaData['Class'] = spec[rankCount] #this must remain 'Class', 'class' has meaning in python
					
					for name in taxonName:
						if ranking == name:
							metaData[name] = spec[rankCount]
					rankCount = rankCount + 1
		
	####################CLEANSING and FILE WRITING LOOP#########################################
	# Loop through list removing any non-alphanumeric characters (exceptions '_,:') that may have slipped by & make certain any unassigned variables are assigned 'NA' 
	#This first line is only used for testing puruposes (OutputMetatData) & should be commented out for actual usage or removed from output file before tree construction.

		keylist = metaData.keys() #dictionaries print in some order, but it is not in the order it was created nor is it ordered alphanumerically.
		keylist.sort()	 #this creates an ordered (alpha) list

		if fileCounter == 0:
			#When the output file is opened up excel, it will allow the user to see the header name for each column to know what variables they are looking at.
			#OutputMetaData.write('lineStarter' + '\t' + 'uniqueID' + '\t') #comment this out once testing is complete

			headerLabelFile.write('uniqueID' + '\n')
			for key in keylist:
				#OutputMetaData.write(key + '\t')
				headerLabelFile.write(key + '\n')
		#OutputMetaData.write('\n')
		headerLabelFile.close()
		
		OutputMetaData.write(lineStarter + '\t' + uniqueID + '\t')
		for values in keylist:
			values = metaData[values]
			values = str(values)
			if values == '' or not values:
				values = 'NA'
			values = values.replace('\n', '')
			values = values.strip( '\t\n' ) #removes leading/trailing whitespace, tab, comma, newline 
			values = values.replace("'", '') #replace whitespaces with underscores '_'
			values = re.sub('\s\s+' , ' ', values) #replaces multiple whitespaces with one whitespace.
			values = values.replace(' ', '_') #replace whitespaces with underscores '_'
			values = re.sub(r'[] []', '', values) #NEW 2.8.2012 remove [] from data (['text'] => 'text'
			OutputMetaData.write(values + '\t') #write variable to file (followed by tab)
		OutputMetaData.write('\n' + str(record.seq) + '\n')
			

	# Writing CDS Sequences to file
		OutputCDSData.write(lineStarter + '\t' + uniqueID + '\t' + file + '\n')
		if not CDS_sequenceList:
			CDS_sequenceList.append('NA')
		cdssequenceList = str(CDS_sequenceList)
		cdssequenceList = re.sub(r'[\',\]\[]', '', cdssequenceList)					
		OutputCDSData.write(cdssequenceList + '\n')
		
		OutputCDSsubfeatureData.write(lineStarter + '\t' + uniqueID + '\t' + file + '\n')
		if not CDS_subfeatureSequenceList:
			CDS_subfeatureSequenceList.append('NA')
		cdssubfeaturesequenceList = str(CDS_subfeatureSequenceList)
		cdssubfeaturesequenceList = re.sub(r'[\' ,\]\[]', '', cdssubfeaturesequenceList)					
		OutputCDSsubfeatureData.write(cdssubfeaturesequenceList + '\n')
		
		OutputExonData.write(lineStarter + '\t' + uniqueID + '\t' + file + '\n')
		if not exon_sequenceList:
			exon_sequenceList.append('NA')
		exonsequenceList = str(exon_sequenceList)
		exonsequenceList = re.sub(r'[ \',\]\[]', '', exonsequenceList)					
		OutputExonData.write(exonsequenceList + '\n')
		
		#augment counter by 1 at the end of each loop
		fileCounter+=1 

		#to re-use dictionary, contents are removed after each loop pass
		del metaData 
		del exon_sequenceList
		del CDS_sequenceList
		del CDS_subfeatureSequenceList

	OutputMetaData.close()
	OutputExonData.close()
	OutputCDSData.close()
	OutputCDSsubfeatureData.close()
