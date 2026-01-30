#! /usr/bin/env python
from ete3 import NCBITaxa
import argparse
import textwrap
import sys
import re
from collections import defaultdict

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
    
    # Initialize NCBI taxdb
    ncbi = NCBITaxa(files.dbfile)
    
    # Load taxid gene lengths
    taxid_genelen = {}
    for line in open(files.taxid_genelens):
        parts = line.strip().split('\t')
        taxid_genelen[parts[0]] = int(parts[1])
    
    # Parse taxid_link file (now 3 columns: busco, taxid, genome_id)
    # taxid_seqs: {taxid: [seq1, seq2, ...]}
    # seq_taxids: {seq: taxid}
    # seq_genomes: {seq: genome_id or "NA"}
    # taxid_genomes: {taxid: set of genome_ids}
    taxid_seqs = defaultdict(list)
    seq_taxids = {}
    seq_genomes = {}
    taxid_genomes = defaultdict(set)
    
    for line in open(files.taxid_link):
        parts = line.strip().split('\t')
        seq = parts[0]
        taxid = parts[1]
        genome_id = parts[2] if len(parts) > 2 else "NA"
        
        taxid_seqs[taxid].append(seq)
        seq_taxids[seq] = taxid
        seq_genomes[seq] = genome_id
        if genome_id != "NA":
            taxid_genomes[taxid].add(genome_id)
    
    # Parse read counts file
    # taxid_counts: {taxid: [[seq, count, correct_bases, total_bases, seqlen, coverage, pid, busco]]}
    # genome_counts: {(taxid, genome): [[seq, count, correct_bases, total_bases, seqlen, coverage, pid, busco]]}
    taxid_counts = defaultdict(list)
    genome_counts = defaultdict(list)
    genuses = defaultdict(list)
    above_species = []
    counter = 0
    
    with open(files.readcounts) as countfile:
        countfile.readline()  # Skip header
        
        for line in countfile:
            counter += 1
            parts = line.strip().split('\t')
            seq = parts[0]
            count = int(parts[1])
            correct_bases = int(parts[2])
            incorrect_bases = int(parts[3])
            total_bases = int(parts[4])
            subjlen = int(parts[5])
            coverage = float(parts[6])
            pid = float(parts[7])
            
            taxid = seq_taxids[seq]
            genome = seq_genomes.get(seq, "NA")
            
            # Extract BUSCO identifier
            if "Collapse" not in seq:
                busco_match = re.findall(r'-\d+at\d+-', seq)
                busco = busco_match[0].strip('-') if busco_match else "Unknown"
            else:
                busco = "Collapsed"
            
            # Determine genus
            lineage = ncbi.get_lineage(int(taxid))
            ranks = {value: key for (key, value) in ncbi.get_rank(lineage).items()}
            
            if 'genus' in ranks and 'Collapse' not in seq and "species" in ranks:
                genus = ranks['genus']
                if taxid not in genuses[genus]:
                    genuses[genus].append(taxid)
            elif "SSCollapse" not in seq:
                above_species.append(taxid)
            
            seq_data = [seq, count, correct_bases, total_bases, subjlen, coverage, pid, busco]
            taxid_counts[taxid].append(seq_data)
            
            # Also track by genome
            if genome != "NA":
                genome_counts[(taxid, genome)].append(seq_data)
    
    if counter == 0:
        message = "Empty read count file. Likely no aligned reads in sample."
        for outfile in [files.eukfrac, files.alltab, files.primarytab]:#, files.marker_list]:
            with open(outfile, 'w') as f:
                f.write(message + '\n')
        sys.exit()
    
    # Calculate stats for each taxid
    taxon_coverage = {}
    seen_taxids = []
    taxon_observed_genomes = defaultdict(set)  # Track which genomes had alignments
    
    for tax in taxid_counts:
        mc = len(taxid_counts[tax])
        counts = sum(entry[1] for entry in taxid_counts[tax])
        correct = sum(entry[2] for entry in taxid_counts[tax])
        total_bases = sum(entry[3] for entry in taxid_counts[tax])
        subj_len = sum(entry[4] for entry in taxid_counts[tax])
        
        buscos = [entry[7] for entry in taxid_counts[tax] if len(entry[7]) > 1]
        
        # Collect genomes that had alignments for this taxid
        for entry in taxid_counts[tax]:
            seq = entry[0]
            genome_id = seq_genomes.get(seq, "NA")
            if genome_id != "NA":
                taxon_observed_genomes[tax].add(genome_id)
        
        percent_identity = round((correct / total_bases) * 100, 2) if total_bases > 0 else 0
        overall_coverage = round((total_bases / subj_len) * 100, 2) if subj_len > 0 else 0
        total_markers = len(taxid_seqs[tax])
        marker_percentage = round(mc / total_markers * 100, 2) if total_markers > 0 else 0
        
        taxid_len = taxid_genelen.get(tax, 0)
        seen_taxids.append(tax)
        
        taxon_coverage[tax] = [mc, counts, total_bases, marker_percentage, overall_coverage, 
                               percent_identity, subj_len, buscos, taxid_len, total_markers]
    
    # Build tree structure
    tree = ncbi.get_topology(seen_taxids)
    tree_root = tree.get_tree_root().name
    lineage = ncbi.get_lineage(tree_root)
    tree_taxids = seen_taxids + lineage
    full_tree = ncbi.get_topology(tree_taxids, intermediate_nodes=True)
    full_taxid_lineage = [node.name for node in full_tree.traverse()]
    
    # Parse inherited markers
    full_seq_taxids = {}
    for line in open(files.inherited_markers):
        parts = line.strip().split('\t')
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
    
    # Determine primary and secondary hits
    primary = {}
    secondary = {}
    genus_secondary_to_primary = {}  # Maps secondary taxid -> primary taxid(s)
    
    for genus, taxids in genuses.items():
        if len(taxids) > 1:
            reads = [taxon_coverage[taxid][1] for taxid in taxids]
            bases = [taxon_coverage[taxid][2] for taxid in taxids]
            
            maxreads = max(reads)
            maxbases = max(bases)
            ptaxids = []
            
            # Determine primary taxids
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
            
            # Check remaining taxids for secondary status
            unsorted_ataxids = [t for t in taxids if t not in ptaxids]
            ataxids = sorted(unsorted_ataxids, key=lambda x: taxon_coverage[x][1], reverse=True)
            
            for ataxid in ataxids:
                is_secondary = False
                responsible_primaries = []
                
                for ptaxid in ptaxids:
                    p_buscos = full_seq_taxids.get(ptaxid, [[], 0, 0])[0]
                    a_buscos = taxon_coverage[ataxid][7]
                    a_remain = [b for b in a_buscos if b in p_buscos]
                    
                    if len(a_remain) > 0:
                        a_above = []
                        for b in a_remain:
                            apid = [seq[6] for seq in taxid_counts[ataxid] if seq[7] == b]
                            ppid = [seq[6] for seq in taxid_counts[ptaxid] if seq[7] == b]
                            
                            if len(ppid) > 0 and apid[0] >= ppid[0]:
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
    
    # Add above-species taxids
    for t in above_species:
        if t not in primary:
            primary[t] = taxon_coverage[t][0:5]
    
    # If > 1 genome per species, attempt read reassignment as for genuses
    
    species_genome_mapping = defaultdict(lambda: defaultdict(dict))  # species_genome_mapping[sptaxid] = {genomeid: {reads, bases, markers, buscos, seq_data}}
    
    # For each species, collect stats per genome
    for (taxid, genome), seq_data_list in genome_counts.items():
        total_reads = sum(seq[1] for seq in seq_data_list)
        total_bases = sum(seq[3] for seq in seq_data_list)
        num_markers = len(seq_data_list)
        buscos = [seq[7] for seq in seq_data_list if len(seq[7]) > 1]
        
        species_genome_mapping[taxid][genome] = {
            'reads': total_reads,
            'bases': total_bases,
            'markers': num_markers,
            'buscos': buscos,
            'seq_data': seq_data_list
        }
    species_primary_genomes = defaultdict(list)  # species_primary_genomes[taxid] = [genome1, ..., genomeN]
    genome_secondary_to_primary = {}  # genome_secondary_to_primary[(taxid, secondary_genome)] = [primary_genome]
    # For each species with multiple genomes, determine primary genome/genomes
    for taxid, genomes_dict in species_genome_mapping.items():
        if len(genomes_dict) > 1:
            genome_list = list(genomes_dict.keys())
            reads_list = [genomes_dict[g]['reads'] for g in genome_list]
            bases_list = [genomes_dict[g]['bases'] for g in genome_list]
            
            maxreads = max(reads_list)
            maxbases = max(bases_list)
            primary_genomes = []
            
            # Determine primary genome
            if (reads_list.count(maxreads) == 1 and bases_list.count(maxbases) == 1 and
                reads_list.index(maxreads) == bases_list.index(maxbases)):
                max_genome = genome_list[reads_list.index(maxreads)]
                primary_genomes.append(max_genome)
            else:
                for i, g in enumerate(genome_list):
                    if reads_list[i] == maxreads or bases_list[i] == maxbases:
                        primary_genomes.append(g)
            
            species_primary_genomes[taxid] = primary_genomes
            
            # Check remaining genomes for secondary status
            secondary_genomes = [g for g in genome_list if g not in primary_genomes]
            
            for sec_genome in secondary_genomes:
                is_secondary = False
                responsible_primaries = []
                
                sec_buscos = genomes_dict[sec_genome]['buscos']
                
                for pri_genome in primary_genomes:
                    pri_buscos = genomes_dict[pri_genome]['buscos']
                    overlap = [b for b in sec_buscos if b in pri_buscos]
                    
                    if len(overlap) > 0:
                        # Use same logic as species-level filtering
                        if len(sec_buscos) < 5:
                            if len(overlap) >= len(sec_buscos):
                                is_secondary = True
                                responsible_primaries.append(pri_genome)
                        else:
                            if len(overlap) > len(sec_buscos) / 2:
                                is_secondary = True
                                responsible_primaries.append(pri_genome)
                
                if is_secondary:
                    genome_secondary_to_primary[(taxid, sec_genome)] = responsible_primaries
        else:
            # Only one genome for this species
            genome = list(genomes_dict.keys())[0]
            species_primary_genomes[taxid] = [genome]
    
    #Reassign reads from secondary genomes to primary genomes
    genome_reassignment_details = defaultdict(lambda: {
        'reassigned_reads': 0,
        'reassigned_correct_bases': 0,
        'reassigned_total_bases': 0,
        'reassigned_marker_length': 0,
        'reassigned_from_genomes': set()
    })
    
    # Store original taxid_counts before genome reassignment
    pre_genome_reassignment_counts = {taxid: list(seqs) for taxid, seqs in taxid_counts.items()}
    
    for (taxid, sec_genome), pri_genomes in genome_secondary_to_primary.items():
        sec_seq_data_list = species_genome_mapping[taxid][sec_genome]['seq_data']
        num_primaries = len(pri_genomes)
        
        # Build BUSCO mapping for primary genomes
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
        
        # Process each marker from secondary genome
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
                # Check if this BUSCO exists in the primary genome
                busco_exists = (busco != "Collapsed" and len(busco) > 1 and 
                               busco in primary_genome_buscos[pri_genome])
                
                if busco_exists:
                    # BUSCO exists - reassign reads to existing marker, don't add length
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
                        # Shouldn't happen, but add it
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
                    # BUSCO doesn't exist in primary genome - add full marker
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
    
    # Filter primary hits
    filter_passing_taxids = []
    filter_failing_taxids = []
    for tax in primary:
        mc = taxon_coverage[tax][0]
        counts = taxon_coverage[tax][1]
        if int(mc) >= 2 and int(counts) >= 4:
            filter_passing_taxids.append(tax)
        else:
            filter_failing_taxids.append(tax)
    
    # Build a mapping of genus -> primary species within that genus
    genus_to_primary_species = defaultdict(list)
    genus_to_all_species = defaultdict(list)
    
    for genus, taxids in genuses.items():
        genus_to_all_species[genus] = taxids
        for taxid in taxids:
            if taxid in filter_passing_taxids:
                genus_to_primary_species[genus].append(taxid)
    
    # Extend genus_secondary_to_primary to include all non-passing species in each genus
    for genus in genuses:
        if genus in genus_to_primary_species and len(genus_to_primary_species[genus]) > 0:
            primary_species_in_genus = genus_to_primary_species[genus]
            all_species_in_genus = genus_to_all_species[genus]
            
            # Find all species that should reassign to primaries
            for species in all_species_in_genus:
                if species not in filter_passing_taxids:
                    # This species either failed filters or was marked secondary - reassign to all primary species
                    if species not in genus_secondary_to_primary:
                        genus_secondary_to_primary[species] = primary_species_in_genus
    
    # Write full read table
    marker_sorted = sorted(taxon_coverage.keys(), reverse=True, 
                          key=lambda x: taxon_coverage[x][3])
    
    with open(files.alltab, 'w') as dest:
        dest.write("Name\tTaxid\tRank\tObserved_markers\tRead_counts\t"
                  "Percent_observed_markers\tTotal_marker_coverage\t"
                  "Percent_identity\tAmount of marker length in EukDetect db\tGenomes\tFiltered\n")
        
        for tax in marker_sorted:
            rank = list(ncbi.get_rank([tax]).values())[0]
            name = list(ncbi.get_taxid_translator([tax]).values())[0]
            
            if rank == "no rank":
                parent = ncbi.get_lineage(tax)[-2]
                rank = list(ncbi.get_rank([parent]).values())[0]
            
            mc, counts, _, marker_percentage, overall_coverage, percent_identity, _, _, blen, _ = taxon_coverage[tax]
            
            # Get genome list - show all genomes that had alignments
            genomes = sorted(list(taxon_observed_genomes.get(tax, set())))
            genome_str = ",".join(genomes) if genomes else "NA"
            
            # Determine if this taxid was filtered out
            is_filtered = "No" if tax in filter_passing_taxids else "Yes"
            
            dest.write(f"{name}\t{tax}\t{rank}\t{mc}\t{counts}\t{marker_percentage}%\t"
                      f"{overall_coverage}%\t{percent_identity}%\t{blen}\t{genome_str}\t{is_filtered}\n")
    
    # Check if any taxa passed filtering
    if len(filter_passing_taxids) == 0:
        message = "No taxa passing filter requirements."
        for outfile in [files.primarytab, files.eukfrac]:#, files.marker_list]:
            with open(outfile, 'w') as f:
                f.write(message + '\n')
        sys.exit()
    
    # Reassign reads from secondary species to primary species - must occur after species level reassignment
    filtered_out_genomes = defaultdict(set)
    
    # Track reassignment details for reporting species level reassignment
    # reassignment_details[pri_taxid] = {
    #   'reassigned_reads': total reads reassigned,
    #   'reassigned_bases': total correct bases from reassigned reads,
    #   'reassigned_total_bases': total bases from reassigned reads,
    #   'reassigned_marker_length': total marker length added from reassignments,
    #   'reassigned_genomes': set of genomes from reassigned species
    # }
    species_reassignment_details = defaultdict(lambda: {
        'reassigned_reads': 0,
        'reassigned_correct_bases': 0,
        'reassigned_total_bases': 0,
        'reassigned_marker_length': 0,
        'reassigned_genomes': set()
    })
    
    # Build a complete copy of taxid_counts
    original_taxid_counts = {taxid: list(seqs) for taxid, seqs in taxid_counts.items()}
    
    for sec_taxid, pri_taxids in genus_secondary_to_primary.items():
        # Track genomes from secondary (filtered) taxids
        if sec_taxid in taxon_observed_genomes:
            for genome in taxon_observed_genomes[sec_taxid]:
                filtered_out_genomes[sec_taxid].add(genome)
                # Add to reassigned genomes for all primaries
                for pri_taxid in pri_taxids:
                    species_reassignment_details[pri_taxid]['reassigned_genomes'].add(genome)
        
        # Get all sequences from secondary hit (use original counts)
        if sec_taxid not in original_taxid_counts:
            continue
            
        sec_seqs = original_taxid_counts[sec_taxid]
        
        # Calculate how to split reads among primaries
        num_primaries = len(pri_taxids)
        
        # Process each marker from the secondary taxid
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
            
            # Distribute reads to primary taxids
            for pri_taxid in pri_taxids:
                # Initialize taxid_counts if needed
                if pri_taxid not in taxid_counts:
                    taxid_counts[pri_taxid] = []
                
                # Check if this exact marker sequence already exists in the primary
                marker_found = False
                for existing_seq_data in taxid_counts[pri_taxid]:
                    if existing_seq_data[0] == marker_seq:
                        # Same marker sequence exists add reads not markerlength
                        existing_seq_data[1] += reads_per_primary
                        existing_seq_data[2] += correct_per_primary
                        existing_seq_data[3] += total_per_primary
                        marker_found = True
                        
                        species_reassignment_details[pri_taxid]['reassigned_reads'] += reads_per_primary
                        species_reassignment_details[pri_taxid]['reassigned_correct_bases'] += correct_per_primary
                        species_reassignment_details[pri_taxid]['reassigned_total_bases'] += total_per_primary
                        break
                
                if not marker_found:
                    # Add reads and markerlenght
                    reassigned_seq = list(seq_data)
                    reassigned_seq[1] = reads_per_primary
                    reassigned_seq[2] = correct_per_primary
                    reassigned_seq[3] = total_per_primary
                    taxid_counts[pri_taxid].append(reassigned_seq)
                    
                    species_reassignment_details[pri_taxid]['reassigned_reads'] += reads_per_primary
                    species_reassignment_details[pri_taxid]['reassigned_correct_bases'] += correct_per_primary
                    species_reassignment_details[pri_taxid]['reassigned_total_bases'] += total_per_primary
                    species_reassignment_details[pri_taxid]['reassigned_marker_length'] += marker_len
    
    # Build tree for filtered taxids (only if there are passing taxa)
    tree = ncbi.get_topology(filter_passing_taxids)
    tree_root = tree.get_tree_root().name
    lineage = ncbi.get_lineage(tree_root)
    primary_tree_taxids = [int(e) for e in filter_passing_taxids] + lineage
    primary_tree = ncbi.get_topology(primary_tree_taxids, intermediate_nodes=True)
    
    # Calculate lineages
    relab_levels = {'species': [], 'genus': [], 'family': [], 'order': [], 'class': [], 'phylum': []}
    ordered_labels = ["phylum", "class", "order", "family", "genus", "species"]
    lineages = {}
    
    # Calculate reads and marker lengths for each node. Includes all marker sequences (even if same BUSCO, different sequences count)
    node_total_reads = defaultdict(float)
    node_total_marker_length = defaultdict(float)
    node_marker_sequences = defaultdict(list)  # Track all marker sequences
    
    # First pass: build lineages and collect data for all nodes in tree
    for node in primary_tree.traverse():
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
        
        # For all filter-passing nodes, collect marker data post reassignemnt
        if node.name in filter_passing_taxids:
            if node.name in taxid_counts:
                for seq_data in taxid_counts[node.name]:
                    marker_seq = seq_data[0]
                    marker_len = seq_data[4]
                    read_count = seq_data[1]
                    
                    node_total_reads[node.name] += read_count
                    node_total_marker_length[node.name] += marker_len
                    node_marker_sequences[node.name].append([marker_seq, marker_len, read_count])
    
    # Propagate markers and reads up the tree
    node_observed_genomes = defaultdict(set)  # Track genomes for filtered nodes only
    
    for node in primary_tree.traverse("postorder"):
        # Only add genomes if this node passed filtering
        if node.name in filter_passing_taxids and node.name in taxon_observed_genomes:
            node_observed_genomes[node.name].update(taxon_observed_genomes[node.name])
        
        if not node.is_leaf():
            # Sum all marker data from descendants
            for desc in node.iter_descendants():
                # Propagate genomes up only from descendants that passed filtering
                if desc.name in filter_passing_taxids and desc.name in node_observed_genomes:
                    node_observed_genomes[node.name].update(node_observed_genomes[desc.name])
                
                # Sum reads and marker lengths from all descendants
                if desc.name in node_total_reads:
                    node_total_reads[node.name] += node_total_reads[desc.name]
                    node_total_marker_length[node.name] += node_total_marker_length[desc.name]
                
                # Add all marker sequences from descendants
                if desc.name in node_marker_sequences:
                    node_marker_sequences[node.name].extend(node_marker_sequences[desc.name])
    
    # Remove incomplete levels
    levels_to_remove = []
    for tax in filter_passing_taxids:
        lin = lineages[tax]
        groups = [l.split('-')[0] for l in lin.split('|')]
        levels = ordered_labels[0:len(groups)]
        if levels != groups:
            for l in levels:
                if l not in groups and l not in levels_to_remove:
                    levels_to_remove.append(l)
    
    for l in levels_to_remove:
        relab_levels.pop(l, None)
    
    # Calculate normalized abundance for each taxonomic level
    relabs = {}
    
    # Calculate RPKS for the lowest level with data
    lowest_rank = None
    for rank in reversed(ordered_labels):  # Start from species, go up
        if rank in relab_levels and len(relab_levels[rank]) > 0:
            lowest_rank = rank
            break
    
    if lowest_rank:
        lowest_level_rpks = {}
        for tax in relab_levels[lowest_rank]:
            total_reads = node_total_reads.get(tax, 0)
            total_marker_length = node_total_marker_length.get(tax, 0)
            
            if total_marker_length > 0:
                rpks = total_reads / (total_marker_length / 1000)
            else:
                rpks = 0.0
            
            lowest_level_rpks[tax] = rpks
        
        # Calculate total RPKS across all taxa at the lowest level for normalization
        total_rpks_lowest = sum(lowest_level_rpks.values())
        
        # Calculate relative abundance at lowest level
        for tax in relab_levels[lowest_rank]:
            rpks = lowest_level_rpks.get(tax, 0.0)
            
            if total_rpks_lowest > 0:
                rel_abundance = (rpks / total_rpks_lowest) * 100
            else:
                rel_abundance = 0.0
            
            total_reads = node_total_reads.get(tax, 0)
            total_marker_length = node_total_marker_length.get(tax, 0)
            num_sequences = len(node_marker_sequences.get(tax, []))
            
            relabs[tax] = [rpks, total_marker_length, rel_abundance, total_reads, num_sequences]
    
    # For all higher levels, sum the relative abundances of their children
    rank_order = [r for r in ordered_labels if r in relab_levels]
    
    for rank_idx in range(len(rank_order) - 2, -1, -1):  # Go from second-lowest up to highest
        rank = rank_order[rank_idx]
        child_rank = rank_order[rank_idx + 1] if rank_idx + 1 < len(rank_order) else None
        
        for tax in relab_levels[rank]:
            # Find all children of this taxon at the child rank
            if child_rank:
                children = []
                for child_tax in relab_levels[child_rank]:
                    child_lineage = ncbi.get_lineage(int(child_tax))
                    if int(tax) in child_lineage:
                        children.append(child_tax)
                
                # Sum relative abundances from children
                rel_abundance = sum(relabs[child][2] for child in children if child in relabs)
            else:
                # No children at a defined rank
                rel_abundance = 0.0
            
            total_reads = node_total_reads.get(tax, 0)
            total_marker_length = node_total_marker_length.get(tax, 0)
            
            if total_marker_length > 0:
                rpks = total_reads / (total_marker_length / 1000)
            else:
                rpks = 0.0
            
            num_sequences = len(node_marker_sequences.get(tax, []))
            
            relabs[tax] = [rpks, total_marker_length, rel_abundance, total_reads, num_sequences]
    
    # Write primary table
    with open(files.primarytab, 'w') as dest:
        dest.write("Name\tRank\tLineage\tTaxid\tObserved_markers\tRead_counts\t"
                  "Percent_observed_markers\tTotal_marker_coverage\tPercent_identity\t"
                  "Reassigned_reads\tReassigned_PID\tReassigned_genomes\t"
                  "Total_reads_with_reassigned\tTotal_marker_length_with_reassigned\tGenomes\n")
        
        for tax in filter_passing_taxids:
            lin = lineages.get(tax, list(ncbi.get_rank([tax]).values())[0])
            rank = list(ncbi.get_rank([tax]).values())[0]
            
            if rank == "no rank":
                parent = ncbi.get_lineage(tax)[-2]
                prevrank = list(ncbi.get_rank([parent]).values())[0]
                if prevrank == "species":
                    rank = "species"
            
            name = list(ncbi.get_taxid_translator([tax]).values())[0]
            mc, counts, _, marker_percentage, overall_coverage, percent_identity, *_ = taxon_coverage[tax]
            
            # Combine genome-level and species-level reassignment details
            genome_reassigned_reads = genome_reassignment_details[tax]['reassigned_reads']
            genome_reassigned_correct = genome_reassignment_details[tax]['reassigned_correct_bases']
            genome_reassigned_total = genome_reassignment_details[tax]['reassigned_total_bases']
            genome_reassigned_marker_length = genome_reassignment_details[tax]['reassigned_marker_length']
            genome_reassigned_genomes = genome_reassignment_details[tax]['reassigned_from_genomes']
            
            species_reassigned_reads = species_reassignment_details[tax]['reassigned_reads']
            species_reassigned_correct = species_reassignment_details[tax]['reassigned_correct_bases']
            species_reassigned_total = species_reassignment_details[tax]['reassigned_total_bases']
            species_reassigned_marker_length = species_reassignment_details[tax]['reassigned_marker_length']
            species_reassigned_genomes = species_reassignment_details[tax]['reassigned_genomes']
            
            # Total reassignment from both genome and species levels
            total_reassigned_reads = round(genome_reassigned_reads + species_reassigned_reads, 2)
            total_reassigned_correct = genome_reassigned_correct + species_reassigned_correct
            total_reassigned_total = genome_reassigned_total + species_reassigned_total
            total_reassigned_marker_length = genome_reassigned_marker_length + species_reassigned_marker_length
            
            # Calculate combined reassigned PID
            if total_reassigned_total > 0:
                reassigned_pid = round((total_reassigned_correct / total_reassigned_total) * 100, 2)
            else:
                reassigned_pid = 0.0
            
            # Combine reassigned genomes from both levels
            all_reassigned_genomes = sorted(list(genome_reassigned_genomes | species_reassigned_genomes))
            reassigned_genome_str = ",".join(all_reassigned_genomes) if all_reassigned_genomes else "NA"
            
            # Calculate total reads and marker length (original + reassigned)
            total_reads_with_reassigned = round(node_total_reads.get(tax, 0), 2)
            total_marker_length_with_reassigned = int(node_total_marker_length.get(tax, 0))
            
            # Get genome list (only from filter-passing taxa)
            genomes = sorted(list(node_observed_genomes.get(tax, set())))
            genome_str = ",".join(genomes) if genomes else "NA"
            
            dest.write(f"{name}\t{rank}\t{lin}\t{tax}\t{mc}\t{counts}\t"
                      f"{marker_percentage}%\t{overall_coverage}%\t{percent_identity}%\t"
                      f"{total_reassigned_reads}\t{reassigned_pid}%\t{reassigned_genome_str}\t"
                      f"{total_reads_with_reassigned}\t{total_marker_length_with_reassigned}\t{genome_str}\n")
    
    # Write relative abundance table
    with open(files.eukfrac, 'w') as dest:
        dest.write("Lineage\tRank\tName\tTaxID\tRPKS\tTotal_marker_length\t"
                  "Num_marker_sequences\tRelative_abundance\tTotal_reads\tGenomes\n")
        
        for node in primary_tree.traverse("preorder"):
            rank = list(ncbi.get_rank([node.name]).values())[0]
            name = list(ncbi.get_taxid_translator([node.name]).values())[0]
            lin = lineages.get(node.name, rank)
            
            if rank == "no rank" and node.is_leaf():
                continue
            
            if node.name in relabs:
                rpks = round(relabs[node.name][0], 4)
                total_marker_length = int(relabs[node.name][1])
                rel_abundance = round(relabs[node.name][2], 4)
                total_reads = round(relabs[node.name][3], 2)
                num_sequences = relabs[node.name][4]
                
                # Get genome list for this node
                genomes = sorted(list(node_observed_genomes.get(node.name, set())))
                genome_str = ",".join(genomes) if genomes else "NA"
                
                dest.write(f"{lin}\t{rank}\t{name}\t{node.name}\t{rpks}\t"
                          f"{total_marker_length}\t{num_sequences}\t{rel_abundance}\t{total_reads}\t{genome_str}\n")
    

if __name__ == "__main__":
    main(sys.argv)
