import os
import subprocess
import glob
import pickle
import gffutils
import numpy
import pandas

from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
from collections import Counter


# paths
genomes_dir = 'insert_path_here'
hog_fasta_dir = 'insert_path_here'
results_dir = 'insert_path_here'


# inputs
species_filename = dict([('DROME','GCF_000001215.4_Release_6_plus_ISO1_MT'), ('DROSI','GCF_000754195.2_ASM75419v2'), ('DROSE','GCF_000005215.3_dsec_caf1'), ('DROYA','GCF_000005975.2_dyak_caf1'), ('DROER','GCF_000005135.1_dere_caf1')])
species_sciname = dict([('DROME','Drosophila melanogaster'), ('DROSI','Drosophila simulans'), ('DROSE','Drosophila sechellia'), ('DROYA','Drosophila yakuba'), ('DROER','Drosophila erecta')])
clade = ['DROSI', 'DROSE', 'DROME']
clade_sci = ['Drosophila simulans', 'Drosophila sechellia', 'Drosophila melanogaster']
out_group = ['DROYA', 'DROER']
ordered_species = ['DROSI', 'DROSE', 'DROME', 'DROYA', 'DROER']
trash_dir_name = 'non_trgs_with_homologs'
trash_dir_name2 = 'non_trgs_with_orfs'

# blast settings
aa_blast_db = 'insert_path_here'
blast_out_format = '7 qseqid qstart qend evalue sseqid sstart send sstrand bitscore score length pident qcovhsp qcovs qseq sseq qlen slen'
cpu = '2'
ecutoff = 1e-03
lencutoff = 0.5






### create input files for oma from files downloaded from NCBI

os.chdir(genomes_dir)
oma_to_refseq = {}
refseq_to_oma = {}
list_of_identifiers = []

# for each of the species
for k,v in species_filename.items():
	print(k)
	# create gff db
	# db = gffutils.create_db('%s_genomic.gff'%(v), dbfn='%s_genomic.gffdb'%(v), force=False, keep_order=False, merge_strategy='merge', sort_attribute_values=False)
	# or open an existing one
	db = gffutils.FeatureDB('%s_genomic.gffdb'%(v), keep_order=True)
	# open protein fasta
	record_dict = SeqIO.index('%s_protein.faa'%(v), 'fasta')
	# add proteins to list to check for duplicates
	list_of_identifiers.extend(list(record_dict.keys()))
	# create files for OMA and record protein lengths
	fasta = ''
	splice_file = ''
	counter = 1
	genelens = []
	# loop through all genes and find it's proteins
	genes = db.features_of_type(featuretype='gene')
	for gene in tqdm(genes, miniters=0,mininterval=1,maxinterval=5):
		proteins = []
		# loop through mRNA associated with it and collect protein IDs for annotated CDSs
		for mrna in db.children(gene, featuretype='mRNA'):
			for cds in db.children(mrna, featuretype='CDS'):
				proteins.append(cds['protein_id'][0])
		proteins = list(set(proteins))
		# if this gene produces proteins, add them to files
		if len(proteins) > 0:
			protein_lens = []
			omaids = []
			# assign OMA IDs to these proteins
			for protein in proteins:
				omaid = '%s%s'%(k, "{:05}".format(counter))
				omaids.append(omaid)
				counter += 1
				# add them to fasta and to oma_to_refseq dictionary
				fasta += '>%s\n%s\n'%(omaid, str(record_dict[protein].seq))
				oma_to_refseq[omaid] = protein
				refseq_to_oma[protein] = omaid
				protein_lens.append(len(str(record_dict[protein].seq)))
			# record the mean lengths of a protein for this gene
			genelens.append(numpy.mean(protein_lens))
			# and add this gene to splice file if more than one mRNA associated with this gene
			if len(proteins) > 1:
				splice_file += ';'.join(omaids) + '\n'
	# write out OMA input files
	with open('%s.fa'%(k), 'w') as f:
		f.write(fasta)
	with open('%s.splice'%(k), 'w') as f:
		f.write(splice_file)
	print('%s - %s genes'%(k, str(len(genelens))))
# write out dictionaries
with open('oma_to_refseq.pickle', 'wb') as f:
	pickle.dump(oma_to_refseq, f)
with open('refseq_to_oma.pickle', 'wb') as f:
	pickle.dump(refseq_to_oma, f)
# check for identifier duplicates
duplicates = [x for x, y in Counter(list_of_identifiers).items() if y > 1]
print(duplicates)









### Run OMA

### find putative TRG HOGs based on OMA output and write out aa fasta with original identifiers

with open('%s/oma_to_refseq.pickle'%(genomes_dir), 'rb') as f:
	oma_to_refseq = pickle.load(f)
os.chdir(results_dir)
hogs = glob.glob('%s/*.fa'%(hog_fasta_dir))
putative_trg_hogs1 = []
# loop throug all HOGs
for hog in tqdm(hogs):
	records = list(SeqIO.parse(hog, 'fasta'))
	hog_members = []
	for record in records:
		hog_members.append(record.description.split()[-1][1:-1])
	hog_members = list(set(hog_members))
	outgroup_members = [x for x in hog_members if x not in clade]
	if len(outgroup_members) == 0:
		hog_name = hog.split('/')[-1].split('.')[0]
		putative_trg_hogs1.append(hog_name)
		subprocess.run("mkdir %s"%(hog_name), shell=True)
		os.chdir('%s/%s'%(results_dir,hog_name))
		faa = ''
		for record in records:
			faa += '>%s [%s]\n%s\n' %(oma_to_refseq[record.id], record.id[:5], str(record.seq))
		with open('%s.faa' %(hog_name), 'w') as f:
			f.write(faa)
with open('trg_families_based_on_oma_hogs.txt', 'w') as f:
	f.write(','.join(putative_trg_hogs1))


















### Run BLASTp

my_path = 'insert_path_here'

hogs = glob.glob('HOG*')
blast_out_format = '7 qseqid qstart qend evalue sseqid sstart send sstrand bitscore score length pident qcovhsp qcovs qseq sseq qlen slen sscinames'

for hog in hogs:
    with open('blastp.bash', 'r') as f:
        script = f.read()
    script += "cd %s/%s\n"%(my_path, hog)
    script += "module add Blast/ncbi-blast/2.7.1+;\n"
    script += "blastp -db refseq -query %s.faa -out %s_vs_refseq.blastp -num_threads 4 -evalue 1 -max_target_seqs 500 -word_size 3 -gapopen 11 -gapextend 1 -matrix BLOSUM62 -threshold 11 -comp_based_stats 2 -seg no -soft_masking false -xdrop_ungap 20 -xdrop_gap 30 -xdrop_gap_final 25 -window_size 40 -outfmt '%s';\n"%(hog, hog, blast_out_format)
    with open('blastp_%s.bash'%(hog), 'w') as f:
        f.write(script)
    cmd = 'bsub < ./blastp_%s.bash'%(hog)
    subprocess.run(cmd, shell=True)





### Validate putative TRG HOGs with blastp search against OMA


for hog in tqdm(putative_trg_hogs1, miniters=0,mininterval=1,maxinterval=5):
	os.chdir('%s/%s'%(results_dir,hog))
	if not os.path.isfile('%s_vs_oma.blastp'%(hog)):
		cmd = "blastp -db %s -query %s.faa -out %s_vs_oma.blastp -num_threads %s -evalue 1 -max_target_seqs 500 -word_size 3 -gapopen 11 -gapextend 1 -matrix BLOSUM62 -threshold 11 -comp_based_stats 2 -seg no -soft_masking false -xdrop_ungap 20 -xdrop_gap 30 -xdrop_gap_final 25 -window_size 40 -outfmt '%s'"%(aa_blast_db, hog, hog, cpu, blast_out_format)
		subprocess.run(cmd, shell=True)

putative_trg_hogs2 = []
for hog in putative_trg_hogs1:
	gene_hits ={}
	with open('%s/%s/%s_vs_oma.blastp'%(results_dir,hog,hog), 'r') as f:
		for line in f.readlines():
			if not line.startswith('#'):
				this_line = line.split()
				if this_line[0] not in gene_hits:
					gene_hits[this_line[0]] = []
				# if not in the clade
				if this_line[4][:5] not in clade:
					# if the query to target alignment covers most of the protein (50%)
					if float(this_line[10])/int(this_line[16]) > lencutoff:
						# if evalue is < 0.001
						if float(this_line[3]) < ecutoff:
							gene_hits[this_line[0]].append(this_line[4])
		number_of_hits_per_gene = [len(v) for k, v in gene_hits.items()]
		# if not every gene has a hit as discribed above, add hog to trg list
		if 0 in number_of_hits_per_gene:
			putative_trg_hogs2.append(hog)
with open('%s/putative_trg_families_without_blastp_hits_in_oma.txt'%(results_dir), 'w') as f:
	f.write(','.join(putative_trg_hogs2))

## move discarded HOGs to a folder
os.chdir('%s'%(results_dir))
cmd = 'mkdir putative_trg_families_with_blastp_hits_in_oma'
subprocess.run(cmd, shell=True)
non_trgs = [x for x in putative_trg_hogs1 if x not in putative_trg_hogs2]
for hog in non_trgs:
	cmd = 'mv %s %s/putative_trg_families_with_blastp_hits_in_oma'%(hog,results_dir)
	subprocess.run(cmd, shell=True)



### Validate putative TRG HOGs with blastp search against RefSeq (on the cluster)
cluster_path = 'insert_path_here'

for hog in putative_trg_hogs2:
	cmd = 'scp -r %s %s'%(hog, cluster_path)
	subprocess.run(cmd, shell=True)

for hog in putative_trg_hogs2:
	cmd = 'scp %s/%s/%s_vs_refseq.blastp %s'%(cluster_path, hog, hog, hog)
	subprocess.run(cmd, shell=True)



putative_trg_hogs3 = []
for hog in tqdm(putative_trg_hogs2):
	gene_hits ={}
	with open('%s/%s/%s_vs_refseq.blastp'%(results_dir,hog,hog), 'r') as f:
		for line in f.readlines():
			if not line.startswith('#'):
				this_line = line.split()
				if this_line[0] not in gene_hits:
					gene_hits[this_line[0]] = []
				# if not in the clade
				if this_line[4][:5] not in clade:
					# if the query to target alignment covers most of the protein (50%)
					if float(this_line[10])/int(this_line[16]) > lencutoff:
						# if evalue is < 0.001
						if float(this_line[3]) < ecutoff:
							gene_hits[this_line[0]].append(this_line[4].split('|')[-1])
	for k, v in gene_hits.items():
		cmd = 'efetch -format docsum -db protein -id %s | xtract -pattern DocumentSummary -element Title'%(','.join(list(set(v))))
		species = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
		species = species.decode('utf-8').split('\n')
		species = [x.split('[')[-1].split(']')[0] for x in species if len(x) > 0]
		gene_hits[k] = [x for x in species if x not in clade_sci]
	number_of_hits_per_gene = [len(v) for k, v in gene_hits.items()]
	if 0 in number_of_hits_per_gene:
		putative_trg_hogs3.append(hog)
with open('%s/putative_trg_families_without_blastp_hits_in_ref_seq.txt'%(results_dir), 'w') as f:
	f.write(','.join(putative_trg_hogs3))

## move discarded HOGs to a folder
os.chdir('%s'%(results_dir))
cmd = 'mkdir putative_trg_families_with_blastp_hits_in_ref_seq'
subprocess.run(cmd, shell=True)
non_trgs = [x for x in putative_trg_hogs2 if x not in putative_trg_hogs3]
for hog in non_trgs:
	cmd = 'mv %s %s/putative_trg_families_with_blastp_hits_in_oma'%(hog,results_dir)
	subprocess.run(cmd, shell=True)







### Collect sequences from genome fasta files

# creat df with TRG genes
trg_genes = pandas.DataFrame(columns=['HOG', 'Species', 'Protein', 'Gene', 'fna'])
# add genes from TRG HOGs to df
for hog in putative_trg_hogs3:
	records = list(SeqIO.parse('%s/%s/%s.faa'%(results_dir, hog, hog), 'fasta'))
	for record in records:
		trg_genes.loc[len(trg_genes)]=[hog, record.description.split('[')[1][:-1], record.id, '', '']
trg_genes.to_csv('%s/trg_genes.tsv'%(results_dir), sep='\t', index=False)

# add info from gff files
for key, value in species_filename.items():
	print(key)
	genome_fasta = SeqIO.index('%s/%s_genomic.fna'%(genomes_dir, value), 'fasta')
	db = gffutils.FeatureDB('%s/%s_genomic.gffdb'%(genomes_dir, value), keep_order=True)
	cds_list = list(db.features_of_type('CDS'))
	for index, row in tqdm(trg_genes.iterrows()):
		# for each TRG of this species
		if trg_genes.iloc[index]['Species'] == key:
			protein = trg_genes.iloc[index]['Protein']
			# find CDSs of this TRG
			cds_of_this_protein = [x for x in cds_list if x['protein_id'][0] == protein]
			# check that all it's CDSs belong to a single gene
			gene = []
			for cds in cds_of_this_protein:
				parents = list(db.parents(cds, featuretype='gene'))
				gene.extend(parents)
			gene = list(set(gene))
			# if no mistakes 
			if len(set(gene)) == 1:
				# add gene id to df
				trg_genes.at[index,'Gene'] = gene[0].id
				# find all CDS associated with this gene
				cds_of_this_gene = list(db.children(gene[0], featuretype='CDS', order_by='start'))
				# add DNA sequences of the gene and it's CDSs to fasta TAKING POSITIVE STRAND SEQ REGARDLESS OF GENE STRAND
				chromosome = gene[0].seqid
				fna = '>%s [%s] %s:%s-%s \n%s\n' %(protein, key, chromosome, gene[0].start, gene[0].end, str(genome_fasta[chromosome].seq)[gene[0].start-1:gene[0].end])
				cds_coordinates = []
				for cds in cds_of_this_gene:
					coords = str(cds.start)+str(cds.end)
					if coords not in cds_coordinates:
						cds_coordinates.append(coords)
						fna += '>%s_cds_%s [%s] %s:%s-%s \n%s\n' %(protein, len(cds_coordinates), key, chromosome, cds.start, cds.end, str(genome_fasta[chromosome].seq)[cds.start-1:cds.end])
				trg_genes.at[index,'fna'] = fna
			else:
				print('Error in %s!'%(protein))
				break	
trg_genes.to_csv('%s/trg_genes.tsv'%(results_dir), sep='\t', index=False)

# write out DNA fasta files
for hog in tqdm(putative_trg_hogs3):
	fna = ''
	for index, row in trg_genes.iterrows():
		if trg_genes.iloc[index]['HOG'] == hog:
			fna += trg_genes.iloc[index]['fna']
	with open('%s/%s/%s.fna' %(results_dir, hog, hog), 'w') as f:
		f.write(fna)






### BLASTN against all genomes and store identified regions in a dataframe

# Make blast db from 5 drosophila genome
os.chdir(genomes_dir)
for k,v in species_filename.items():
	print(k)
	cmd = 'makeblastdb -in %s_genomic.fna -dbtype nucl'%(v)
	subprocess.run(cmd, shell=True)

# blastn against each genome
for hog in tqdm(putative_trg_hogs3):
	os.chdir("%s/%s"%(results_dir, hog))
	for k,v in species_filename.items():
		cmd = "blastn -db %s/%s_genomic.fna -query %s.fna -out %s_vs_%s.blastn -num_threads %s -evalue 1 -word_size 7 -max_target_seqs 500 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 -strand both -perc_identity 0 -outfmt '%s'"%(genomes_dir, v, hog, hog, k, cpu, blast_out_format)
		subprocess.run(cmd, shell=True)

# add hits from blastn to trg_hits dictionary (including self-hits)
trg_hits = {}
for hog in tqdm(putative_trg_hogs3):
	for k,v in species_filename.items():
		with open('%s/%s/%s_vs_%s.blastn'%(results_dir, hog, hog, k), 'r') as f:
			for line in f.readlines():
				if not line.startswith('#'):
					tl = line.split()
					# filter out low confidence hits
					if float(tl[3]) <= ecutoff:
						# and hits that cover insufficient proportion of the query sequence
						if float(this_line[10])/int(this_line[16]) > lencutoff:
							start = str(min(int(tl[5]), int(tl[6])))
							end = str(max(int(tl[5]), int(tl[6])))
							trg_hits.setdefault('%s\t%s\t%s'%(hog, k, tl[4]), []).append([int(start), int(end), float(tl[3])])

# create df with genome regions
trg_df = pandas.DataFrame(columns=['HOG', 'Species', 'Chromosome', 'Start', 'End', 'E-value', 'Sequence'])
for k,v in tqdm(trg_hits.items()):
	# sort intervals by start coordinate
	locations = sorted(v, key=lambda x: x[0])
	# identify unique locations
	locations_unique = []
	# iterate over the list merging overlapping intervals
	current_location = locations[0]
	for i in range(1, len(locations)):
		# if next interval starts before previous ended, merge them
		if locations[i][0] < current_location[1]:
			current_location = [current_location[0], max(current_location[1], locations[i][1]), min(current_location[2], locations[i][2])]
		else:
		# else add it to unique locations
			locations_unique.append(current_location)
			current_location = locations[i]
	# add the final current interval to unique locations
	locations_unique.append(current_location)
	# sort intervals by e-value
	locations_unique = sorted(locations_unique, key=lambda x: x[2])
	for i in range(len(locations_unique)):
		trg_df.loc[len(trg_df)]=[k.split()[0], k.split()[1], k.split()[2], locations_unique[i][0], locations_unique[i][1], locations_unique[i][2], '']
trg_df.to_csv('%s/trg_hogs.tsv'%(results_dir), sep='\t', index=False)
# # Check how many genome regions each of the HOGs is associated with 
# trg_df.groupby('HOG').agg('count')

# look up sequences from the genome fasta files
genome_fasta = {}
for k,v in species_filename.items():
	genome_fasta[k] = SeqIO.index('%s/%s_genomic.fna'%(genomes_dir, v), 'fasta')

species = ''
chromosome = ''
sequence = ''
for index, row in tqdm(trg_df.iterrows()):
	if trg_df.iloc[index]['Chromosome'] == chromosome:
		if trg_df.iloc[index]['Species'] == species:
			# ensure that annotations don't spill out outside the sequence boundaries
			start = max(0, trg_df.iloc[index]['Start']-1)
			end = min(len(sequence), trg_df.iloc[index]['End'])
			trg_df.at[index,'Sequence'] = sequence[start:end]
	else:
		species = trg_df.iloc[index]['Species']
		chromosome = trg_df.iloc[index]['Chromosome']
		sequence = str(genome_fasta[species][chromosome].seq)
		# ensure that annotations don't spill out outside the sequence boundaries
		start = max(0, trg_df.iloc[index]['Start']-1)
		end = min(len(sequence), trg_df.iloc[index]['End'])
		trg_df.at[index,'Sequence'] = sequence[start:end]
trg_df.to_csv('%s/trg_hogs.tsv'%(results_dir), sep='\t', index=False)

# write homologous DNA fasta files
for hog in putative_trg_hogs3:
	print(hog)
	# create a dictionary with species and sequences
	species_chromosome_fna = {}
	for index, row in tqdm(trg_df.iterrows()):
		if trg_df.iloc[index]['HOG'] == hog:
			species_chromosome_fna.setdefault('%s-%s'%(trg_df.iloc[index]['Species'],trg_df.iloc[index]['Chromosome']), []).append('%s:%s\n%s\n' %(trg_df.iloc[index]['Start'], trg_df.iloc[index]['End'], trg_df.iloc[index]['Sequence']))
	# create fasta file of genome regions
	fna = ''
	for species in ordered_species:
		for key, value in species_chromosome_fna.items():
			if key.startswith(species):
				for entry in value:
					fna += '>%s-%s'%(key, entry)
	with open('%s/%s/%s_homolog.fna' %(results_dir, hog, hog), 'w') as f:
		f.write(fna)



## for hogs with more than 1000 hits, load tsv of all hits, select 5 best per species and create a new fasta file
trg_df = pandas.read_csv('%s/trg_hogs.tsv'%(results_dir), sep='\t')
for hog in ['HOG13860','HOG13765','HOG13746','HOG13754','HOG13766','HOG13767','HOG14054']:
	fna = ''
	hog_df = trg_df.loc[trg_df['HOG'] == hog]
	for species in ordered_species:
		sub_hog_df = hog_df.loc[hog_df['Species'] == species]
		sub_hog_df = sub_hog_df.sort_values('E-value')[:5]
		for index, row in sub_hog_df.iterrows():
			fna += '>%s-%s-%s:%s\n%s\n'%(trg_df.iloc[index]['Species'],trg_df.iloc[index]['Chromosome'],trg_df.iloc[index]['Start'], trg_df.iloc[index]['End'], trg_df.iloc[index]['Sequence'])
	with open('%s/%s/%s_top_5_homolog.fna' %(results_dir, hog, hog), 'w') as f:
		f.write(fna)








### After manual inspection...

## for two candidate TRGs run blastn for extra 4 outgroups
extra_outgroup_species_filename = dict([('DROAN','GCF_000005115.1_dana_caf1'), ('DROSU','GCF_000472105.1_Dsuzukii.v01'), ('DROPS','GCF_000001765.3_Dpse_3.0'), ('DROMI','GCF_000269505.1_DroMir_2.2')])
for hog in ['HOG14059', 'HOG14096']:
	print(hog)
	os.chdir("%s/%s"%(results_dir, hog))
	for k,v in extra_outgroup_species_filename.items():
		print(k)
		cmd = "blastn -db %s/%s_genomic.fna -query %s.fna -out %s_vs_%s.blastn -num_threads %s -evalue 1 -word_size 7 -max_target_seqs 500 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 -strand both -perc_identity 0 -outfmt '%s'"%(genomes_dir, v, hog, hog, k, cpu, blast_out_format)
		subprocess.run(cmd, shell=True)

## collect hits from 4 extra outgroups
trg_hits_2 = {}
for hog in ['HOG14059', 'HOG14096']:
	print(hog)
	for k,v in extra_outgroup_species_filename.items():
		with open('%s/%s/%s_vs_%s.blastn'%(results_dir, hog, hog, k), 'r') as f:
			for line in f.readlines():
				if not line.startswith('#'):
					tl = line.split()
					# filter out low confidence hits
					if float(tl[3]) <= ecutoff:
						# and hits that cover insufficient proportion of the query sequence
						if float(tl[10])/int(tl[16]) > lencutoff:
							start = str(min(int(tl[5]), int(tl[6])))
							end = str(max(int(tl[5]), int(tl[6])))
							trg_hits_2.setdefault('%s\t%s\t%s'%(hog, k, tl[4]), []).append([int(start), int(end), float(tl[3])])

## add those extra hits to trg_df
for k,v in tqdm(trg_hits_2.items()):
	# sort intervals by start coordinate
	locations = sorted(v, key=lambda x: x[0])
	# identify unique locations
	locations_unique = []
	# iterate over the list merging overlapping intervals
	current_location = locations[0]
	for i in range(1, len(locations)):
		# if next interval starts before previous ended, merge them
		if locations[i][0] < current_location[1]:
			current_location = [current_location[0], max(current_location[1], locations[i][1]), min(current_location[2], locations[i][2])]
		else:
		# else add it to unique locations
			locations_unique.append(current_location)
			current_location = locations[i]
	# add the final current interval to unique locations
	locations_unique.append(current_location)
	# sort intervals by e-value
	locations_unique = sorted(locations_unique, key=lambda x: x[2])
	for i in range(len(locations_unique)):
		trg_df.loc[len(trg_df)]=[k.split()[0], k.split()[1], k.split()[2], locations_unique[i][0], locations_unique[i][1], locations_unique[i][2], '']
trg_df.to_csv('%s/trg_hogs_2.tsv'%(results_dir), sep='\t', index=False)

# look up sequences from the genome fasta files
genome_fasta = {}
for k,v in species_filename.items():
	genome_fasta[k] = SeqIO.index('%s/%s_genomic.fna'%(genomes_dir, v), 'fasta')
for k,v in extra_outgroup_species_filename.items():
	genome_fasta[k] = SeqIO.index('%s/%s_genomic.fna'%(genomes_dir, v), 'fasta')

ordered_outgroups = ['DROAN', 'DROSU', 'DROPS', 'DROMI']
species = ''
chromosome = ''
sequence = ''
for index, row in tqdm(trg_df.iterrows()):
	if trg_df.iloc[index]['Species'] in ordered_outgroups:
		species = trg_df.iloc[index]['Species']
		chromosome = trg_df.iloc[index]['Chromosome']
		sequence = str(genome_fasta[species][chromosome].seq)
		# ensure that annotations don't spill out outside the sequence boundaries
		start = max(0, trg_df.iloc[index]['Start']-1)
		end = min(len(sequence), trg_df.iloc[index]['End'])
		trg_df.at[index,'Sequence'] = sequence[start:end]
trg_df.to_csv('%s/trg_hogs.tsv'%(results_dir), sep='\t', index=False)

# write homologous DNA fasta files for 2 candidate HOGs
for hog in ['HOG14059', 'HOG14096']:
	print(hog)
	# create a dictionary with species and sequences
	species_chromosome_fna = {}
	for index, row in tqdm(trg_df.iterrows()):
		if trg_df.iloc[index]['HOG'] == hog:
			species_chromosome_fna.setdefault('%s-%s'%(trg_df.iloc[index]['Species'],trg_df.iloc[index]['Chromosome']), []).append('%s:%s\n%s\n' %(trg_df.iloc[index]['Start'], trg_df.iloc[index]['End'], trg_df.iloc[index]['Sequence']))
	# create fasta file of genome regions
	fna = ''
	for species in ordered_outgroups:
		for key, value in species_chromosome_fna.items():
			if key.startswith(species):
				for entry in value:
					fna += '>%s-%s'%(key, entry)
	with open('%s/%s/%s_homolog_extra_outgroup.fna' %(results_dir, hog, hog), 'w') as f:
		f.write(fna)





## look up annotations associated with hit regions

for k,v in species_filename.items():
	db = gffutils.FeatureDB('%s/%s_genomic.gffdb'%(genomes_dir, v), keep_order=True)
	scaffold = 'NT_167065.1'
	start = 1351813
	end = 1352284
	annotations = list(db.region(region=(scaffold, start, end), completely_within=False))

	# db.children(annotation)





## write element hits to file
trgf_hog_dir = 'insert_path_here'
element_hit_file = 'HOG14096_elements_vs_genome.blastn'
genome_dna_seqs_file = 'HOG14096_homolog_for_image.fna'
element_dna_seqs_file = 'HOG14096_elements.fna'


os.chdir(trgf_hog_dir)

## create a dictionary of elements and their hits
element_hits = {}
added_hits = []
with open('%s/%s'%(trgf_hog_dir,element_hit_file), 'r') as f:
	for line in f.readlines():
		if line.startswith('# Query: '):
			element_name = line.strip()[9:]
		elif not line.startswith('#'):
			tl = line.split()
			hit_name = '%s %s'%(tl[4], element_name)
			if hit_name not in added_hits:
				added_hits.append(hit_name)
				element_hits.setdefault(tl[4], []).append([element_name, min([int(tl[5]), int(tl[6])]), max([int(tl[5]), int(tl[6])])])

genome_dna_seqs = SeqIO.index('%s/%s'%(trgf_hog_dir, genome_dna_seqs_file), 'fasta')
element_dna_seqs = SeqIO.index('%s/%s'%(trgf_hog_dir, element_dna_seqs_file), 'fasta')

file_name_list = []

for k,v in element_hits.items():
	hits = element_hits[k]
	hits.sort(key=lambda x: x[1])
	start = hits[0][1]
	end = hits[0][2]
	element_list = [hits[0][0]]
	for i in range(1,len(hits)):
		if hits[i][1] < end:
			end = max(end, hits[i][2])
			element_list.append(hits[i][0])
		else:
			file_name = '%s_%s_%s'%(k.split('-')[0], str(start), str(end))
			file_name_list.append(file_name)
			fasta = '>%s-%s-%s:%s\n%s\n'%(k.split('-')[0], k.split('-')[1], str(start), str(end), str(genome_dna_seqs[k].seq)[start-1:end])
			with open('%s/%s.fna'%(trgf_hog_dir, file_name), 'w') as f:
				f.write(fasta)
			fasta = ''
			for element in element_list:
				fasta += '>%s\n%s\n'%(element, str(element_dna_seqs[element.split()[0]].seq))
			with open('%s/%s_elements.fna'%(trgf_hog_dir, file_name), 'w') as f:
				f.write(fasta)
			start = hits[i][1]
			end = hits[i][2]
			element_list = [hits[i][0]]
	file_name = '%s_%s_%s'%(k.split('-')[0], str(start), str(end))
	file_name_list.append(file_name)
	fasta = '>%s-%s-%s:%s\n%s\n'%(k.split('-')[0], k.split('-')[1], str(start), str(end), str(genome_dna_seqs[k].seq)[start-1:end])
	with open('%s/%s.fna'%(trgf_hog_dir, file_name), 'w') as f:
		f.write(fasta)
	fasta = ''
	for element in element_list:
		fasta += '>%s\n%s\n'%(element, str(element_dna_seqs[element.split()[0]].seq))
	with open('%s/%s_elements.fna'%(trgf_hog_dir, file_name), 'w') as f:
		f.write(fasta)

for file_name in file_name_list:
	cmd = "mafft --genafpair --maxiterate 1000 --thread 3 --adjustdirectionaccurately --multipair --addfragments %s_elements.fna %s.fna > %s.msa"%(file_name, file_name, file_name)
	subprocess.run(cmd, shell=True)