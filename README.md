# A metagenomic framework for rapid Listeria monocytogenes surveillance in food production environments
Listeria monocytogenes remains a major foodborne pathogen with high mortality and costly persistence in food-processing environments. Established diagnostics rely on selective enrichment and single-colony isolation, which can introduce biases by favouring fast-growing or more tolerant strains while suppressing slow-growing or viable-but-nonculturable strains.

To quantify the bias introduced by culture-based diagnostics and assess the advantage of metagenomics- based pathogen detection, we developed and evaluated a rapid nanopore sequencing–based metagenomic framework. We designed a defined mock community comprising hypervirulent,persistent, and low-virulence L. monocytogenes strains, other Listeria species, and a background microbiome representative of food-processing environments,and used this to compare (i) nanopore metagenomic sequencing of the environmental surveillance matrices (ii) quasimetagenomic sequencing after 4, 12, 24, or 48h of selective enrichment, and (iii) established culture-based whole-genome sequencing. We found that culture-based diagnostics recovered a limited subset of strains and failed to detect multi-strain contamination, while metagenomics captured a variety of spiked Listeria strains. 

Quasimetagenomics indicated selective enrichment of specific L. monocytogenes strains while only leading to a relative enrichment of L. monocytogenes after at least 12 h. We finally found that the sensitivity of metagenomics can be improved by utilising different environmental sampling materials, whereas adaptive sampling-based enrichment of L. monocytogenes had no effect. We conclude that integrating nanopore metagenomics into routine surveillance shows great promise to improve detection and source attribution in food systems.

The manuscript is available here

Basecalled sample bam files in ENA - PRJEB129024.

# ONT read based analysis for Listeria detection

The following pipeline is based on read based minimap2 alignment. The pipeline follows step alignment process (i) aligning read length filtered reads against the GTDB database(v232) for overall taxa identification and (ii) extracting all reads on mapped on the Listeria genus level and mapping it against all the available Listeria genomes (6607) in the GTDB database(v232). The second alignment step also includes genome from Listeriaceae family (Brochothrix thermosphacta) as decoy sequences, which is included in the Listeria genomes. Samtools coverage option is used to obtain the alignment summary stats. The mapped bases are calculated as follows: mapped bases = meandepth × (endpos − startpos + 1) (reference genome length corrected).

# Indirect ST inference

MLST was further performed on all L. monocytogenes genomes (n=5,748) obtained from GTDB using the mlst tool (v2.35) against the listeria_2 PubMLST scheme. Genomes were supplied as a list of assembly file paths, and sequence types (STs) were determined using the --full option to report the complete allelic profile for each locus in the scheme in addition to the assigned ST. To link the coverage table obtained for each L. monocytogenes genome to a known ST, the coverage tables were annotated with the ST of the corresponding reference assembly using a custom Python script. For each row, the genome assembly accession was extracted from the relevant column and cross-referenced against the MLST results generated above. 
Note: Reads and bases assigned to ST9 and ST35 were merged and reported as ST35 in downstream analyses.

# Scripts

[`align.sh`](align.sh) - Slurm script for preprocessing and alignment.

add_taxa.py - GTDB annotation of the samtools coverage tsv files.

background_genus_summary.py - relative abundance (reads and bases) along with absolute read counts and mapped bases for background genus.

listeria_annotate.py - add ST and Listeria species annotation.

stats_lis.py - relative abundance (reads and bases) along with absolute read counts and mapped bases for Listeria species and ST.
