
from ete3 import NCBITaxa
import argparse
import textwrap
import sys
import re
from collections import defaultdict
import logging


logging.basicConfig(
	format="%(asctime)s [%(levelname)s] %(message)s",
	datefmt="%Y-%m-%d %H:%M:%S",
	level=logging.INFO
)
logger = logging.getLogger(__name__)


def main(argv):
	parser = argparse.ArgumentParser(
		description=textwrap.dedent("""\
			Summarize and filter alignments by taxid with multi-genome support.
			Uses read-based abundance with length normalization for universal markers.
			Required arguments are --dbfile, --inherited_markers, --taxid_link, 
			--readcounts, --primarytab, --eukfrac, --alltab, --taxid_genelens
		"""),
		formatter_class = argparse.RawDescriptionHelpFormatter
	)
	parser.add_argument("--dbfile", type=str, required=True, help="Eukdetect database folder")
	parser.add_argument("--inherited_markers", type=str, required=True, help="Inherited markers file")
	parser.add_argument("--taxid_link", type=str, required=True, help="BUSCO to taxid to genome link file (3 columns)")
	parser.add_argument("--readcounts", type=str, required=True, help="Read counts and mismatches file")
	parser.add_argument("--eukfrac", type=str, required=True, help="Eukaryotic abundance & fraction output file")
	parser.add_argument("--primarytab", type=str, required=True, help="Table output of filtered hits")
	parser.add_argument("--alltab", type=str, required=True, help="Table output of all hits")
	parser.add_argument("--taxid_genelens", type=str, required=True, help="Cumulative gene length per taxid")
	
	files = parser.parse_args()
	
	try:
		ncbi = NCBITaxa(files.dbfile)
		logger.info("Successfully loaded NCBI taxonomy database")
	except Exception as e:
		logger.error(f"Failed to load NCBI taxonomy database from {files.dbfile}: {e}")
		sys.exit(1)
	
	#Load genome and taxid-specific gene lengths
	#File format: identifier\tcumulative_length
	#Identifier can be either genome name string or taxid

	taxid_genelen = {}  #taxid -> length (species-level)
	genome_genelen = {}  #genome_name -> length (genome-specific)
	
	try:
		with open(files.taxid_genelens) as f:
			for line in f:
				parts = line.strip().split('\t')
				if len(parts) >= 2:
					identifier = parts[0]
					length = int(parts[1])
					
					#Check if identifier is a taxid (numeric) or genome name (string)
					if identifier.isdigit():
						#Taxid
						taxid_genelen[identifier] = length
					else:
						#Genome name
						genome_genelen[identifier] = length
					
		logger.info(f"Loaded gene lengths for {len(taxid_genelen)} taxa and {len(genome_genelen)} genomes")
	except FileNotFoundError:
		logger.error(f"Taxid gene lengths file not found: {files.taxid_genelens}")
		sys.exit(1)
	except Exception as e:
		logger.error(f"Error reading taxid gene lengths: {e}")
		sys.exit(1)
	
	#Parse taxid_link file (busco, taxid, genome_id)
	taxid_seqs = defaultdict(list)
	seq_taxids = {}
	seq_genomes = {}
	taxid_genomes = defaultdict(set)
	
	try:
		with open(files.taxid_link) as f:
			for line in f:
				parts = line.strip().split('\t')
				if len(parts) < 2:
					continue
				seq = parts[0]
				taxid = parts[1]
				genome_id = parts[2] if (len(parts) > 2) else "NA"
				
				taxid_seqs[taxid].append(seq)
				seq_taxids[seq] = taxid
				seq_genomes[seq] = genome_id
				if genome_id != "NA":
					taxid_genomes[taxid].add(genome_id)
		logger.info(f"Loaded {len(seq_taxids)} sequence-to-taxid mappings")
	except FileNotFoundError:
		logger.error(f"Taxid link file not found: {files.taxid_link}")
		sys.exit(1)
	except Exception as e:
		logger.error(f"Error reading taxid link file: {e}")
		sys.exit(1)
	
	#Parse read counts file
	taxid_counts = defaultdict(list)
	genome_counts = defaultdict(list)
	genuses = defaultdict(list)
	above_species = []
	counter = 0
	
	try:
		with open(files.readcounts) as countfile:
			countfile.readline()  #Skip header
			
			for line in countfile:
				counter += 1
				parts = line.strip().split('\t')
				if len(parts) < 8:
					logger.warning(f"Skipping malformed line {counter} in read counts file")
					continue
				
				seq = parts[0]
				
				#Check if sequence is in our mapping
				if seq not in seq_taxids:
					logger.warning(f"Sequence {seq} not found in taxid_link, skipping")
					continue
				
				try:
					count = int(parts[1])
					correct_bases = int(parts[2])
					incorrect_bases = int(parts[3])
					total_bases = int(parts[4])
					subjlen = int(parts[5])
					coverage = float(parts[6])
					pid = float(parts[7])
				except (ValueError, IndexError) as e:
					logger.warning(f"Error parsing line {counter}: {e}, skipping")
					continue
				
				taxid = seq_taxids[seq]
				genome = seq_genomes.get(seq, "NA")
				
				#Extract BUSCO identifier
				if "Collapse" not in seq:
					busco_match = re.findall(r'-\d+at\d+-', seq)
					busco = busco_match[0].strip('-') if busco_match else "Unknown"
				else:
					busco = "Collapsed"
				
				#Determine genus
				try:
					lineage = ncbi.get_lineage(int(taxid))
					ranks = {value: key for (key, value) in ncbi.get_rank(lineage).items()}
					
					if 'genus' in ranks and 'Collapse' not in seq and "species" in ranks:
						genus = ranks['genus']
						if taxid not in genuses[genus]:
							genuses[genus].append(taxid)
					elif "SSCollapse" not in seq:
						above_species.append(taxid)
				except Exception as e:
					logger.warning(f"Could not get lineage for taxid {taxid}: {e}")
					continue
				
				seq_data = [seq, count, correct_bases, total_bases, subjlen, coverage, pid, busco]
				taxid_counts[taxid].append(seq_data)
				
				#Also track by genome
				if genome != "NA":
					genome_counts[(taxid, genome)].append(seq_data)
		
		logger.info(f"Processed {counter} alignment records")
	
	except FileNotFoundError:
		logger.error(f"Read counts file not found: {files.readcounts}")
		sys.exit(1)
	except Exception as e:
		logger.error(f"Error reading read counts file: {e}")
		sys.exit(1)
	
	#End if empty read count file
	if counter == 0:
		message = "Empty read count file. Likely no aligned reads in sample."
		logger.warning(message)
		for outfile in [files.eukfrac, files.alltab, files.primarytab]:
			try:
				with open(outfile, 'w') as f:
					f.write(message + '\n')
			except Exception as e:
				logger.error(f"Error writing to {outfile}: {e}")
		sys.exit()
	
	#Calculate initial stats for each taxid
	logger.info("Calculating initial taxon coverage statistics.")
	taxon_coverage = {}
	seen_taxids = []
	taxon_observed_genomes = defaultdict(set)
	
	for tax in taxid_counts:
		mc = len(taxid_counts[tax])
		counts = sum(entry[1] for entry in taxid_counts[tax])
		correct = sum(entry[2] for entry in taxid_counts[tax])
		total_bases = sum(entry[3] for entry in taxid_counts[tax])
		subj_len = sum(entry[4] for entry in taxid_counts[tax])
		
		buscos = [entry[7] for entry in taxid_counts[tax] if len(entry[7]) > 1]
		
		#Collect genomes that had alignments for this taxid
		for entry in taxid_counts[tax]:
			seq = entry[0]
			genome_id = seq_genomes.get(seq, "NA")
			if genome_id != "NA":
				taxon_observed_genomes[tax].add(genome_id)
		
		percent_identity = round((correct / total_bases) * 100, 2) if total_bases > 0 else 0
		overall_coverage = round((total_bases / subj_len) * 100, 2) if subj_len > 0 else 0
		total_markers = len(taxid_seqs.get(tax, []))
		marker_percentage = round(mc / total_markers * 100, 2) if total_markers > 0 else 0
		
		taxid_len = taxid_genelen.get(tax, 0)
		seen_taxids.append(tax)
		
		taxon_coverage[tax] = [mc, counts, total_bases, marker_percentage, overall_coverage, 
							   percent_identity, subj_len, buscos, taxid_len, total_markers]
	
	
	#Store original (pre-reassignment) statistics for reporting
	#Format: [mc, counts, total_bases, marker_percentage, overall_coverage, 
	#		 percent_identity, subj_len, buscos, taxid_len, total_markers]
	original_taxon_stats = {tax: list(stats) for tax, stats in taxon_coverage.items()}
	
	logger.info("Starting genome-level disambiguation.")
	
	species_genome_mapping = defaultdict(lambda: defaultdict(dict))
	
	#For each species, collect stats per genome
	for (taxid, genome), seq_data_list in genome_counts.items():
		total_reads = sum(seq[1] for seq in seq_data_list)
		total_bases = sum(seq[3] for seq in seq_data_list)
		num_markers = len(seq_data_list)
		buscos = [seq[7] for seq in seq_data_list if len(seq[7]) > 1]
		
		#Calculate overall PID (weighted average by total bases)
		total_correct = sum(seq[2] for seq in seq_data_list)
		total_aligned = sum(seq[3] for seq in seq_data_list)
		overall_pid = (total_correct / total_aligned * 100) if total_aligned > 0 else 0.0
		
		species_genome_mapping[taxid][genome] = {
			'reads': total_reads,
			'bases': total_bases,
			'markers': num_markers,
			'buscos': buscos,
			'seq_data': seq_data_list,
			'overall_pid': overall_pid
		}
	
	species_primary_genomes = defaultdict(list)
	genome_secondary_to_primary = {}
	
	#For each species with multiple genomes, determine primary genome (can be multiple primary)
	for taxid, genomes_dict in species_genome_mapping.items():
		if len(genomes_dict) > 1:
			genome_list = list(genomes_dict.keys())
			reads_list = [genomes_dict[g]['reads'] for g in genome_list]
			bases_list = [genomes_dict[g]['bases'] for g in genome_list]
			
			maxreads = max(reads_list)
			maxbases = max(bases_list)
			primary_genomes = []
			
			#Determine primary genome
			if (reads_list.count(maxreads) == 1 and bases_list.count(maxbases) == 1 and
				reads_list.index(maxreads) == bases_list.index(maxbases)):
				max_genome = genome_list[reads_list.index(maxreads)]
				primary_genomes.append(max_genome)
			else:
				for i, g in enumerate(genome_list):
					if reads_list[i] == maxreads or bases_list[i] == maxbases:
						primary_genomes.append(g)
			
			species_primary_genomes[taxid] = primary_genomes
			
			#Check remaining genomes for secondary status
			secondary_genomes = [g for g in genome_list if g not in primary_genomes]
			
			for sec_genome in secondary_genomes:
				is_secondary = False
				responsible_primaries = []
				
				sec_buscos = genomes_dict[sec_genome]['buscos']
				sec_seq_data = genomes_dict[sec_genome]['seq_data']
				sec_overall_pid = genomes_dict[sec_genome]['overall_pid']
				
				for pri_genome in primary_genomes:
					pri_buscos = genomes_dict[pri_genome]['buscos']
					pri_seq_data = genomes_dict[pri_genome]['seq_data']
					pri_overall_pid = genomes_dict[pri_genome]['overall_pid']
					overlap = [b for b in sec_buscos if b in pri_buscos]
					
					#Compare overall PIDs regardless of overlap. If secondary has worse overall PID, it should be reassigned
					
					if len(overlap) == 0:
						#No shared BUSCOs
						#Check if secondary has worse overall PID than primary
						if sec_overall_pid <= pri_overall_pid:
							#Secondary has worse PID, reassign to primary
							is_secondary = True
							responsible_primaries.append(pri_genome)

					elif len(overlap) > 0:
						#Has overlap - check both overlapping PIDs and overall PID
						if sec_overall_pid <= pri_overall_pid:
							#Secondary has worse overall PID → reassign
							is_secondary = True
							responsible_primaries.append(pri_genome)
						else:
							#Secondary has better overall PID - check specific overlaps
							buscos_with_worse_or_equal_pid = []
							
							for busco in overlap:
								#Get PIDs for this BUSCO
								sec_pids = [s[6] for s in sec_seq_data if s[7] == busco]
								pri_pids = [s[6] for s in pri_seq_data if s[7] == busco]
								
								if len(pri_pids) > 0 and len(sec_pids) > 0:
									#Secondary has worse or equal PID on this BUSCO
									if sec_pids[0] <= pri_pids[0]:
										buscos_with_worse_or_equal_pid.append(busco)
							
							#Threshholding
							if len(sec_buscos) < 5:
								#Few markers: strict threshold
								if len(buscos_with_worse_or_equal_pid) >= 1:
									is_secondary = True
									responsible_primaries.append(pri_genome)
							else:
								#Many markers: lenient threshold
								if len(buscos_with_worse_or_equal_pid) > len(sec_buscos) / 2:
									is_secondary = True
									responsible_primaries.append(pri_genome)
				
				if is_secondary:
					genome_secondary_to_primary[(taxid, sec_genome)] = responsible_primaries
		else:
			#Only one genome for this species
			genome = list(genomes_dict.keys())[0]
			species_primary_genomes[taxid] = [genome]
	
	logger.info(f"Identified {len(genome_secondary_to_primary)} secondary genomes to reassign")
	
	#Calculate genome-specific marker lengths for each taxid based on primary genomes
	taxid_primary_genome_length = {}
	
	for taxid, primary_genomes in species_primary_genomes.items():
		total_length = 0
		found_genomes = []
		
		for genome in primary_genomes:
			#Look up genome-specific length by genome name
			genome_length = genome_genelen.get(genome, 0)
			
			if genome_length > 0:
				total_length += genome_length
				found_genomes.append(genome)
			else:
				logger.warning(f"No genome-specific length found for {genome} (taxid {taxid})")
		
		if total_length > 0:
			#Use genome-specific lengths
			taxid_primary_genome_length[taxid] = total_length
			logger.debug(f"Taxid {taxid}: Using {len(found_genomes)} genome(s) with total {total_length} bp")
		else:
			#Fallback to taxid-level length if no genome-specific lengths found
			taxid_level_length = taxid_genelen.get(taxid, 0)
			if taxid_level_length > 0:
				taxid_primary_genome_length[taxid] = taxid_level_length
				logger.debug(f"Taxid {taxid}: Falling back to taxid-level length {taxid_level_length} bp")
			else:
				logger.warning(f"No marker length found for taxid {taxid} (neither genome nor taxid level)")
	
	logger.info(f"Calculated genome-specific marker lengths for {len(taxid_primary_genome_length)} taxa")
	
	#Reassign reads from secondary genomes to primary genomes
	genome_reassignment_details = defaultdict(lambda: {
		'reassigned_reads': 0,
		'reassigned_correct_bases': 0,
		'reassigned_total_bases': 0,
		'reassigned_marker_length': 0,
		'reassigned_from_genomes': set()
	})
	
	for (taxid, sec_genome), pri_genomes in genome_secondary_to_primary.items():
		sec_seq_data_list = species_genome_mapping[taxid][sec_genome]['seq_data']
		num_primaries = len(pri_genomes)
		
		#Build busco mapping for primary genomes
		primary_genome_buscos = {}
		for pri_genome in pri_genomes:
			primary_genome_buscos[pri_genome] = {}
			pri_seq_data_list = species_genome_mapping[taxid][pri_genome]['seq_data']
			for seq_data in pri_seq_data_list:
				marker_seq = seq_data[0]
				busco = seq_data[7]
				if busco not in primary_genome_buscos[pri_genome]:
					primary_genome_buscos[pri_genome][busco] = []
				primary_genome_buscos[pri_genome][busco].append(marker_seq)
		
		#Process each marker from secondary genome
		for seq_data in sec_seq_data_list:
			marker_seq = seq_data[0]
			busco = seq_data[7]
			count = seq_data[1]
			correct_bases = seq_data[2]
			total_bases = seq_data[3]
			marker_len = seq_data[4]
			
			reads_per_primary = count / num_primaries
			correct_per_primary = correct_bases / num_primaries
			total_per_primary = total_bases / num_primaries
			
			for pri_genome in pri_genomes:
				#Check if this busco exists in the primary genome
				busco_exists = (busco != "Collapsed" and len(busco) > 1 and 
							   busco in primary_genome_buscos[pri_genome])
				
				if busco_exists:
					#reassign reads to existing marker, dont add length
					target_marker = primary_genome_buscos[pri_genome][busco][0]
					
					marker_found = False
					for existing_seq_data in taxid_counts[taxid]:
						if existing_seq_data[0] == target_marker:
							existing_seq_data[1] += reads_per_primary
							existing_seq_data[2] += correct_per_primary
							existing_seq_data[3] += total_per_primary
							marker_found = True
							
							genome_reassignment_details[taxid]['reassigned_reads'] += reads_per_primary
							genome_reassignment_details[taxid]['reassigned_correct_bases'] += correct_per_primary
							genome_reassignment_details[taxid]['reassigned_total_bases'] += total_per_primary
							genome_reassignment_details[taxid]['reassigned_from_genomes'].add(sec_genome)
							break
					
					if not marker_found:
						new_seq = list(seq_data)
						new_seq[1] = reads_per_primary
						new_seq[2] = correct_per_primary
						new_seq[3] = total_per_primary
						taxid_counts[taxid].append(new_seq)
						
						genome_reassignment_details[taxid]['reassigned_reads'] += reads_per_primary
						genome_reassignment_details[taxid]['reassigned_correct_bases'] += correct_per_primary
						genome_reassignment_details[taxid]['reassigned_total_bases'] += total_per_primary
						genome_reassignment_details[taxid]['reassigned_marker_length'] += marker_len
						genome_reassignment_details[taxid]['reassigned_from_genomes'].add(sec_genome)
				else:
					#busco doesn't exist in primary genome, add full marker to total genome length
					new_seq = list(seq_data)
					new_seq[1] = reads_per_primary
					new_seq[2] = correct_per_primary
					new_seq[3] = total_per_primary
					taxid_counts[taxid].append(new_seq)
					
					genome_reassignment_details[taxid]['reassigned_reads'] += reads_per_primary
					genome_reassignment_details[taxid]['reassigned_correct_bases'] += correct_per_primary
					genome_reassignment_details[taxid]['reassigned_total_bases'] += total_per_primary
					genome_reassignment_details[taxid]['reassigned_marker_length'] += marker_len
					genome_reassignment_details[taxid]['reassigned_from_genomes'].add(sec_genome)
	

	logger.info("Recalculating taxon coverage after genome reassignment.")
	
	taxon_coverage = {}
	for tax in taxid_counts:
		mc = len(taxid_counts[tax])
		counts = sum(entry[1] for entry in taxid_counts[tax])
		correct = sum(entry[2] for entry in taxid_counts[tax])
		total_bases = sum(entry[3] for entry in taxid_counts[tax])
		subj_len = sum(entry[4] for entry in taxid_counts[tax])
		
		buscos = [entry[7] for entry in taxid_counts[tax] if len(entry[7]) > 1]
		
		percent_identity = round((correct / total_bases) * 100, 2) if total_bases > 0 else 0
		overall_coverage = round((total_bases / subj_len) * 100, 2) if subj_len > 0 else 0
		total_markers = len(taxid_seqs.get(tax, []))
		marker_percentage = round(mc / total_markers * 100, 2) if total_markers > 0 else 0
		
		taxid_len = taxid_genelen.get(tax, 0)
		
		taxon_coverage[tax] = [mc, counts, total_bases, marker_percentage, overall_coverage, 
							   percent_identity, subj_len, buscos, taxid_len, total_markers]

	logger.info("Building taxonomic tree.")
	
	try:
		tree = ncbi.get_topology(seen_taxids)
		tree_root = tree.get_tree_root().name
		lineage = ncbi.get_lineage(tree_root)
		tree_taxids = seen_taxids + lineage
		full_tree = ncbi.get_topology(tree_taxids, intermediate_nodes=True)
		full_taxid_lineage = [node.name for node in full_tree.traverse()]
		logger.info(f"Built tree with {len(full_taxid_lineage)} nodes")
	except Exception as e:
		logger.error(f"Error building taxonomic tree: {e}")
		sys.exit(1)
	

	logger.info("Parsing markers.")
	full_seq_taxids = {}
	try:
		with open(files.inherited_markers) as f:
			for line in f:
				parts = line.strip().split('\t')
				if len(parts) < 3:
					continue
				taxid = parts[0]
				if taxid in full_taxid_lineage:
					buscos = []
					for seq in parts[1].split(','):
						busco_match = re.findall(r'-\d+at\d+-', seq)
						if busco_match:
							busco = busco_match[0].strip('-')
							if busco not in buscos:
								buscos.append(busco)
					
					specific_count = len(parts[1].split(','))
					sp_and_inherited_count = len(parts[2].split(','))
					full_seq_taxids[taxid] = [buscos, specific_count, sp_and_inherited_count]
		logger.info(f"Loaded inherited markers for {len(full_seq_taxids)} taxa")
	except FileNotFoundError:
		logger.error(f"Inherited markers file not found: {files.inherited_markers}")
		sys.exit(1)
	except Exception as e:
		logger.error(f"Error reading inherited markers: {e}")
		sys.exit(1)
	
	logger.info("Starting genus-level disambiguation.")
	
	primary = {}
	secondary = {}
	genus_secondary_to_primary = {}
	
	for genus, taxids in genuses.items():
		if len(taxids) > 1:
			reads = [taxon_coverage[taxid][1] for taxid in taxids]
			bases = [taxon_coverage[taxid][2] for taxid in taxids]
			
			maxreads = max(reads)
			maxbases = max(bases)
			ptaxids = []
			
			#Determine primary taxids
			if (reads.count(maxreads) == 1 and bases.count(maxbases) == 1 and
				reads.index(maxreads) == bases.index(maxbases)):
				maxtax = taxids[reads.index(maxreads)]
				primary[maxtax] = taxon_coverage[maxtax][0:5]
				ptaxids.append(maxtax)
			else:
				for t in taxids:
					if taxon_coverage[t][1] == maxreads or taxon_coverage[t][2] == maxbases:
						ptaxids.append(t)
						primary[t] = taxon_coverage[t][0:5]
			
			#Check remaining taxids for secondary status
			unsorted_ataxids = [t for t in taxids if t not in ptaxids]
			ataxids = sorted(unsorted_ataxids, key=lambda x: taxon_coverage[x][1], reverse=True)
			
			for ataxid in ataxids:
				is_secondary = False
				responsible_primaries = []
				
				#Get secondary's overall PID
				sec_overall_pid = taxon_coverage[ataxid][5]
				
				for ptaxid in ptaxids:
					#Get primary's overall PID
					pri_overall_pid = taxon_coverage[ptaxid][5]
					
					p_buscos = full_seq_taxids.get(ptaxid, [[], 0, 0])[0]
					a_buscos = taxon_coverage[ataxid][7]
					a_remain = [b for b in a_buscos if b in p_buscos]

					
					if len(a_remain) == 0:
						#Check if secondary has worse overall PID than primary
						if sec_overall_pid <= pri_overall_pid:
							#Secondary has worse PID, reassign to primary
							is_secondary = True
							responsible_primaries.append(ptaxid)
					
					elif len(a_remain) > 0:
						#Has overlap - check both overall PID and overlapping PIDs
						if sec_overall_pid <= pri_overall_pid:
							#Secondary has worse overall PID, reassign
							is_secondary = True
							responsible_primaries.append(ptaxid)
						else:
							#Secondary has better overall PID - check specific overlaps
							a_above = []
							for b in a_remain:
								apid = [seq[6] for seq in taxid_counts[ataxid] if seq[7] == b]
								ppid = [seq[6] for seq in taxid_counts[ptaxid] if seq[7] == b]
								
								if len(ppid) > 0 and len(apid) > 0 and apid[0] >= ppid[0]:
									a_above.append(b)
								elif len(ppid) == 0:
									a_above.append(b)
							
							if len(a_buscos) < 5:
								if len(a_above) < len(a_buscos):
									is_secondary = True
									responsible_primaries.append(ptaxid)
							else:
								if len(a_above) <= len(a_buscos) / 2:
									is_secondary = True
									responsible_primaries.append(ptaxid)
				
				if is_secondary:
					secondary[ataxid] = taxon_coverage[ataxid][0:5] + [responsible_primaries]
					genus_secondary_to_primary[ataxid] = responsible_primaries
				else:
					primary[ataxid] = taxon_coverage[ataxid][0:5]
		else:
			taxid = taxids[0]
			primary[taxid] = taxon_coverage[taxid][0:5]
	
	#Add above-species taxids
	for t in above_species:
		if t not in primary:
			primary[t] = taxon_coverage[t][0:5]
	
	logger.info(f"Genus-level: {len(primary)} primary taxa, {len(secondary)} secondary taxa")

	filter_passing_taxids = []
	filter_failing_taxids = []
	for tax in primary:
		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		if int(mc) >= 2 and int(counts) >= 4:
			filter_passing_taxids.append(tax)
		else:
			filter_failing_taxids.append(tax)
	
	logger.info(f"Filtering: {len(filter_passing_taxids)} taxa passed, {len(filter_failing_taxids)} taxa failed")
	
	#Build genus to primary species mapping
	genus_to_primary_species = defaultdict(list)
	genus_to_all_species = defaultdict(list)
	
	for genus, taxids in genuses.items():
		genus_to_all_species[genus] = taxids
		for taxid in taxids:
			if taxid in filter_passing_taxids:
				genus_to_primary_species[genus].append(taxid)
	
	#Extend secondary mapping to include all non-passing species in genus
	for genus in genuses:
		if genus in genus_to_primary_species and len(genus_to_primary_species[genus]) > 0:
			primary_species_in_genus = genus_to_primary_species[genus]
			all_species_in_genus = genus_to_all_species[genus]
			
			for species in all_species_in_genus:
				if species not in filter_passing_taxids:
					if species not in genus_secondary_to_primary:
						genus_secondary_to_primary[species] = primary_species_in_genus
	

	logger.info("Writing full read table.")
	
	marker_sorted = sorted(taxon_coverage.keys(), reverse=True, 
						  key=lambda x: taxon_coverage[x][3])
	
	try:
		with open(files.alltab, 'w') as dest:
			dest.write("Name\tTaxid\tRank\tObserved_markers\tRead_counts\t"
					  "Total_marker_coverage\tPercent_identity\tTotal_marker_length\t"
					  "Genomes\tFiltered\n")
			
			for tax in marker_sorted:
				try:
					rank = list(ncbi.get_rank([tax]).values())[0]
					name = list(ncbi.get_taxid_translator([tax]).values())[0]
				except Exception as e:
					logger.warning(f"Could not get rank/name for taxid {tax}: {e}")
					continue
				
				if rank == "no rank":
					try:
						parent = ncbi.get_lineage(tax)[-2]
						rank = list(ncbi.get_rank([parent]).values())[0]
					except:
						rank = "no rank"
				
				mc, counts, _, marker_percentage, overall_coverage, percent_identity, _, _, blen, _ = taxon_coverage[tax]
				
				genomes = sorted(list(taxon_observed_genomes.get(tax, set())))
				genome_str = ",".join(genomes) if genomes else "NA"
				
				is_filtered = "No" if tax in filter_passing_taxids else "Yes"
				
				dest.write(f"{name}\t{tax}\t{rank}\t{mc}\t{counts}\t"
						  f"{overall_coverage}%\t{percent_identity}%\t{blen}\t"
						  f"{genome_str}\t{is_filtered}\n")
		
		logger.info(f"Wrote all taxa table to {files.alltab}")
	except Exception as e:
		logger.error(f"Error writing all taxa table: {e}")
		sys.exit(1)

	#if nothing passed filter, exit
	if len(filter_passing_taxids) == 0:
		message = "No taxa passing filter requirements."
		logger.warning(message)
		for outfile in [files.primarytab, files.eukfrac]:
			try:
				with open(outfile, 'w') as f:
					f.write(message + '\n')
			except Exception as e:
				logger.error(f"Error writing to {outfile}: {e}")
		sys.exit()
	

	logger.info("Starting species-level reassignment.")
	
	filtered_out_genomes = defaultdict(set)
	species_reassignment_details = defaultdict(lambda: {
		'reassigned_reads': 0,
		'reassigned_correct_bases': 0,
		'reassigned_total_bases': 0,
		'reassigned_marker_length': 0,
		'reassigned_genomes': set()
	})
	
	#Build complete copy of taxid_counts before reassignment
	original_taxid_counts = {taxid: list(seqs) for taxid, seqs in taxid_counts.items()}
	
	for sec_taxid, pri_taxids in genus_secondary_to_primary.items():
		#Track genomes from secondary (filtered) taxids
		if sec_taxid in taxon_observed_genomes:
			for genome in taxon_observed_genomes[sec_taxid]:
				filtered_out_genomes[sec_taxid].add(genome)
				for pri_taxid in pri_taxids:
					species_reassignment_details[pri_taxid]['reassigned_genomes'].add(genome)
		
		#Get all sequences from secondary hit (use original counts)
		if sec_taxid not in original_taxid_counts:
			continue
			
		sec_seqs = original_taxid_counts[sec_taxid]
		
		#Calculate how to split reads among primaries
		num_primaries = len(pri_taxids)
		
		#Process each marker from the secondary taxid
		for seq_data in sec_seqs:
			marker_seq = seq_data[0]
			busco = seq_data[7]
			count = seq_data[1]
			correct_bases = seq_data[2]
			total_bases = seq_data[3]
			marker_len = seq_data[4]
			
			reads_per_primary = count / num_primaries
			correct_per_primary = correct_bases / num_primaries
			total_per_primary = total_bases / num_primaries
			
			#Distribute reads to primary taxids
			for pri_taxid in pri_taxids:
				#Initialize taxid_counts if needed
				if pri_taxid not in taxid_counts:
					taxid_counts[pri_taxid] = []
				
				#Check if this exact marker sequence already exists in the primary
				marker_found = False
				for existing_seq_data in taxid_counts[pri_taxid]:
					if existing_seq_data[0] == marker_seq:
						#Same marker sequence exists - add reads, not marker length
						existing_seq_data[1] += reads_per_primary
						existing_seq_data[2] += correct_per_primary
						existing_seq_data[3] += total_per_primary
						marker_found = True
						
						species_reassignment_details[pri_taxid]['reassigned_reads'] += reads_per_primary
						species_reassignment_details[pri_taxid]['reassigned_correct_bases'] += correct_per_primary
						species_reassignment_details[pri_taxid]['reassigned_total_bases'] += total_per_primary
						break
				
				if not marker_found:
					#Add reads and marker length
					reassigned_seq = list(seq_data)
					reassigned_seq[1] = reads_per_primary
					reassigned_seq[2] = correct_per_primary
					reassigned_seq[3] = total_per_primary
					taxid_counts[pri_taxid].append(reassigned_seq)
					
					species_reassignment_details[pri_taxid]['reassigned_reads'] += reads_per_primary
					species_reassignment_details[pri_taxid]['reassigned_correct_bases'] += correct_per_primary
					species_reassignment_details[pri_taxid]['reassigned_total_bases'] += total_per_primary
					species_reassignment_details[pri_taxid]['reassigned_marker_length'] += marker_len
	

	logger.info("Recalculating taxon coverage after species reassignment.")
	
	taxon_coverage = {}
	for tax in taxid_counts:
		mc = len(taxid_counts[tax])
		counts = sum(entry[1] for entry in taxid_counts[tax])
		correct = sum(entry[2] for entry in taxid_counts[tax])
		total_bases = sum(entry[3] for entry in taxid_counts[tax])
		subj_len = sum(entry[4] for entry in taxid_counts[tax])
		
		buscos = [entry[7] for entry in taxid_counts[tax] if len(entry[7]) > 1]
		
		percent_identity = round((correct / total_bases) * 100, 2) if total_bases > 0 else 0
		overall_coverage = round((total_bases / subj_len) * 100, 2) if subj_len > 0 else 0
		total_markers = len(taxid_seqs.get(tax, []))
		marker_percentage = round(mc / total_markers * 100, 2) if total_markers > 0 else 0
		
		taxid_len = taxid_genelen.get(tax, 0)
		
		taxon_coverage[tax] = [mc, counts, total_bases, marker_percentage, overall_coverage, 
							   percent_identity, subj_len, buscos, taxid_len, total_markers]
	

	logger.info("Building tree for filtered taxa.")
	
	try:
		tree = ncbi.get_topology(filter_passing_taxids)
		tree_root = tree.get_tree_root().name
		lineage = ncbi.get_lineage(tree_root)
		primary_tree_taxids = [int(e) for e in filter_passing_taxids] + lineage
		primary_tree = ncbi.get_topology(primary_tree_taxids, intermediate_nodes=True)
		logger.info("Tree built successfully")
	except Exception as e:
		logger.error(f"Error building filtered tree: {e}")
		sys.exit(1)

	logger.info("Calculating lineages and node statistics...")
	
	relab_levels = {'species': [], 'genus': [], 'family': [], 'order': [], 'class': [], 'phylum': []}
	ordered_labels = ["phylum", "class", "order", "family", "genus", "species"]
	lineages = {}
	
	node_total_reads = defaultdict(float)
	node_total_marker_length = defaultdict(float)
	node_marker_sequences = defaultdict(list)
	
	#First pass: build lineages and collect data for all nodes in tree
	for node in primary_tree.traverse():
		try:
			currname = list(ncbi.get_taxid_translator([node.name]).values())[0]
			lineage_list = ncbi.get_lineage(node.name)
			names = ncbi.get_taxid_translator(lineage_list)
			ranks = ncbi.get_rank(lineage_list)
			ranks_rev = {ranks[e]: e for e in ranks}
			
			lin_name = ""
			for label in ordered_labels:
				if label in ranks_rev:
					lin_name += f"{label}-{names[ranks_rev[label]]}|"
			lin_name = lin_name.strip('|').replace(' ', "_")
			lineages[node.name] = lin_name
			
			rank = list(ncbi.get_rank([node.name]).values())[0]
			
			if rank in relab_levels:
				relab_levels[rank].append(node.name)
		except Exception as e:
			logger.warning(f"Error processing node {node.name}: {e}")
			continue
		
		#For all filter-passing nodes, collect marker data post reassignment
		if node.name in filter_passing_taxids:
			if node.name in taxid_counts:
				for seq_data in taxid_counts[node.name]:
					marker_seq = seq_data[0]
					marker_len = seq_data[4]
					read_count = seq_data[1]
					
					node_total_reads[node.name] += read_count
					node_total_marker_length[node.name] += marker_len
					node_marker_sequences[node.name].append([marker_seq, marker_len, read_count])
	
	#Track primary genome marker lengths from database for RPKS calculation
	#For species with multiple genomes, use only the detected primary genome
	node_original_marker_length = {}
	for tax in filter_passing_taxids:
		#Use genome-specific length (based on primary genome detected). Defaults to species-level if genome-specific not available
		node_original_marker_length[tax] = taxid_primary_genome_length.get(tax, original_taxon_stats.get(tax, [0]*9 + [0])[8])
	
	#Propagate markers and reads up the tree
	node_observed_genomes = defaultdict(set)
	
	for node in primary_tree.traverse("postorder"):
		#Only add genomes if this node passed filtering
		if node.name in filter_passing_taxids and node.name in taxon_observed_genomes:
			node_observed_genomes[node.name].update(taxon_observed_genomes[node.name])
		
		if not node.is_leaf():

			for child in node.children:
				#Propagate genomes up only from children that passed filtering
				if child.name in filter_passing_taxids and child.name in node_observed_genomes:
					node_observed_genomes[node.name].update(node_observed_genomes[child.name])
				
				#Sum reads and marker lengths from direct children
				if child.name in node_total_reads:
					node_total_reads[node.name] += node_total_reads[child.name]
					node_total_marker_length[node.name] += node_total_marker_length[child.name]
				
				#Sum original marker lengths from direct children
				if child.name in node_original_marker_length:
					if node.name not in node_original_marker_length:
						node_original_marker_length[node.name] = 0
					node_original_marker_length[node.name] += node_original_marker_length[child.name]
				
				#Add all marker sequences from direct children
				if child.name in node_marker_sequences:
					node_marker_sequences[node.name].extend(node_marker_sequences[child.name])
	
	#Remove incomplete levels
	levels_to_remove = []
	for tax in filter_passing_taxids:
		lin = lineages.get(tax, "")
		groups = [l.split('-')[0] for l in lin.split('|')]
		levels = ordered_labels[0:len(groups)]
		if levels != groups:
			for l in levels:
				if l not in groups and l not in levels_to_remove:
					levels_to_remove.append(l)
	
	for l in levels_to_remove:
		relab_levels.pop(l, None)

	logger.info("Calculating normalized abundance (RPKS).")
	
	relabs = {}
	
	#Calculate RPKS for the lowest level with data
	lowest_rank = None
	for rank in reversed(ordered_labels):
		if rank in relab_levels and len(relab_levels[rank]) > 0:
			lowest_rank = rank
			break
	
	if lowest_rank:
		lowest_level_rpks = {}
		for tax in relab_levels[lowest_rank]:
			total_reads = node_total_reads.get(tax, 0)
			#Use ORIGINAL marker length to avoid penalty from sparse reassigned markers
			original_marker_length = node_original_marker_length.get(tax, 0)
			
			if original_marker_length > 0:
				rpks = total_reads / (original_marker_length / 1000)
			else:
				rpks = 0.0
			
			lowest_level_rpks[tax] = rpks
		
		#Calculate total RPKS across all taxa at the lowest level for normalization
		total_rpks_lowest = sum(lowest_level_rpks.values())
		
		#Calculate relative abundance at lowest level
		for tax in relab_levels[lowest_rank]:
			rpks = lowest_level_rpks.get(tax, 0.0)
			
			if total_rpks_lowest > 0:
				rel_abundance = (rpks / total_rpks_lowest) * 100
			else:
				rel_abundance = 0.0
			
			total_reads = node_total_reads.get(tax, 0)
			original_marker_length = node_original_marker_length.get(tax, 0)
			num_sequences = len(node_marker_sequences.get(tax, []))
			
			relabs[tax] = [rpks, original_marker_length, rel_abundance, total_reads, num_sequences]
	
	#For all higher levels, sum the relative abundances of their children
	rank_order = [r for r in ordered_labels if r in relab_levels]
	
	for rank_idx in range(len(rank_order) - 2, -1, -1):
		rank = rank_order[rank_idx]
		child_rank = rank_order[rank_idx + 1] if rank_idx + 1 < len(rank_order) else None
		
		for tax in relab_levels[rank]:
			#Find all children of this taxon at the child rank
			if child_rank:
				children = []
				for child_tax in relab_levels[child_rank]:
					try:
						child_lineage = ncbi.get_lineage(int(child_tax))
						if int(tax) in child_lineage:
							children.append(child_tax)
					except:
						continue
				
				#Sum relative abundances from children
				rel_abundance = sum(relabs[child][2] for child in children if child in relabs)
			else:
				rel_abundance = 0.0
			
			total_reads = node_total_reads.get(tax, 0)
			#Use original marker length (sum of children's original lengths)
			original_marker_length = node_original_marker_length.get(tax, 0)
			
			if original_marker_length > 0:
				rpks = total_reads / (original_marker_length / 1000)
			else:
				rpks = 0.0
			
			num_sequences = len(node_marker_sequences.get(tax, []))
			
			relabs[tax] = [rpks, original_marker_length, rel_abundance, total_reads, num_sequences]

	logger.info("Writing primary table.")
	
	try:
		with open(files.primarytab, 'w') as dest:
			dest.write("Name\tRank\tLineage\tTaxid\tTotal_reads\tTotal_marker_length\t"
					  "RPKS\tReads_aligned\tPID_aligned\tGenomes\t"
					  "Reads_reassigned\tPID_reassigned\tReassigned_genomes\n")
			
			for tax in filter_passing_taxids:
				try:
					lin = lineages.get(tax, list(ncbi.get_rank([tax]).values())[0])
					rank = list(ncbi.get_rank([tax]).values())[0]
					
					if rank == "no rank":
						parent = ncbi.get_lineage(tax)[-2]
						prevrank = list(ncbi.get_rank([parent]).values())[0]
						if prevrank == "species":
							rank = "species"
					
					name = list(ncbi.get_taxid_translator([tax]).values())[0]
				except Exception as e:
					logger.warning(f"Error getting info for taxid {tax}: {e}")
					continue
				
				#Get original alignment stats (before any reassignment)
				mc_original = original_taxon_stats[tax][0]  #marker count
				reads_original = original_taxon_stats[tax][1]  #read count
				pid_original = original_taxon_stats[tax][5]  #percent identity
				
				total_genome_marker_length = taxid_primary_genome_length.get(tax, original_taxon_stats[tax][8])
				
				#Combine genome-level and species-level reassignment details
				genome_reassigned_reads = genome_reassignment_details[tax]['reassigned_reads']
				genome_reassigned_correct = genome_reassignment_details[tax]['reassigned_correct_bases']
				genome_reassigned_total = genome_reassignment_details[tax]['reassigned_total_bases']
				genome_reassigned_genomes = genome_reassignment_details[tax]['reassigned_from_genomes']
				
				species_reassigned_reads = species_reassignment_details[tax]['reassigned_reads']
				species_reassigned_correct = species_reassignment_details[tax]['reassigned_correct_bases']
				species_reassigned_total = species_reassignment_details[tax]['reassigned_total_bases']
				species_reassigned_genomes = species_reassignment_details[tax]['reassigned_genomes']
				
				#Total reassigned reads and PID
				total_reassigned_reads = round(genome_reassigned_reads + species_reassigned_reads, 2)
				total_reassigned_correct = genome_reassigned_correct + species_reassigned_correct
				total_reassigned_total = genome_reassigned_total + species_reassigned_total
				
				reassigned_pid = round((total_reassigned_correct / total_reassigned_total) * 100, 2) if total_reassigned_total > 0 else 0.0
				
				#Build full list of reassigned genomes
				all_reassigned_genomes = sorted(list(genome_reassigned_genomes | species_reassigned_genomes))
				reassigned_genome_str = ",".join(all_reassigned_genomes) if all_reassigned_genomes else "None"
				
				#Calculate RPKS using total reads with full genome marker length from database
				total_reads_all = reads_original + total_reassigned_reads
				rpks = round(total_reads_all / (total_genome_marker_length / 1000), 4) if total_genome_marker_length > 0 else 0.0
				
				#Get genomes detected for this taxid
				genomes = sorted(list(node_observed_genomes.get(tax, set())))
				genome_str = ",".join(genomes) if genomes else "NA"
				
				dest.write(f"{name}\t{rank}\t{lin}\t{tax}\t"
						  f"{round(total_reads_all, 2)}\t{int(total_genome_marker_length)}\t{rpks}\t"
						  f"{round(reads_original, 2)}\t{pid_original}%\t{genome_str}\t"
						  f"{total_reassigned_reads}\t{reassigned_pid}%\t{reassigned_genome_str}\n")
		
		logger.info(f"Wrote primary table to {files.primarytab}")
	except Exception as e:
		logger.error(f"Error writing primary table: {e}")
		sys.exit(1)
	
	logger.info("Writing relative abundance table.")
	
	try:
		with open(files.eukfrac, 'w') as dest:
			dest.write("Lineage\tRank\tName\tTaxID\tRPKS\t"
					  "Relative_abundance\tReads_total\tTotal_marker_length\tGenomes\n")
			
			for node in primary_tree.traverse("preorder"):
				try:
					rank = list(ncbi.get_rank([node.name]).values())[0]
					name = list(ncbi.get_taxid_translator([node.name]).values())[0]
					lin = lineages.get(node.name, rank)
				except Exception as e:
					logger.warning(f"Error getting info for node {node.name}: {e}")
					continue
				
				if rank == "no rank" and node.is_leaf():
					continue
				
				if node.name in relabs:
					#Only report RPKS for species-level taxa
					if rank == "species":
						rpks = round(relabs[node.name][0], 4)
					else:
						rpks = "NA"
					
					total_marker_length = int(relabs[node.name][1])
					rel_abundance = round(relabs[node.name][2], 4)
					total_reads = round(relabs[node.name][3], 2)
					
					genomes = sorted(list(node_observed_genomes.get(node.name, set())))
					genome_str = ",".join(genomes) if genomes else "NA"
					
					dest.write(f"{lin}\t{rank}\t{name}\t{node.name}\t{rpks}\t"
							  f"{rel_abundance}\t{total_reads}\t{total_marker_length}\t{genome_str}\n")
		
		logger.info(f"Wrote relative abundance table to {files.eukfrac}")
	except Exception as e:
		logger.error(f"Error writing relative abundance table: {e}")
		sys.exit(1)
	
	logger.info("Analysis complete!")


if __name__ == "__main__":
	main(sys.argv)
