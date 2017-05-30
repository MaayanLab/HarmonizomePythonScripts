=================================
Mouse dorsal LGN Cell Types study 
=================================

rnaseq_profile_id

RNA sequencing data of single cells isolated from mouse dorsal lateral geniculate nucleus (dLGN) of the thalamus.

The data set includes 1772 single cells collected from dLGN in adult mouse and transcriptionally profiled with RNA sequencing. Additional donor meta-data and broad and subclass labels are available for each cell.

The sequencing results were aligned and aggregated at the gene level using the RSEM algorithm, and FPKM values were calculated.

For more details, please see the Documentation tab in the Cell Types web application.


fpkm_table.csv
	Contains the (row, column) matrix of fpkm values obtained for each (gene, cell).
	The first row contains the unique identifiers of the RNA-seq profiles of the cells (rnaseq_profile_id)
	The first column contains the gene unique identifiers (gene_id)


columns-cells.csv 
	Contains information about the cells profiled with RNA sequencing
	
	rnaseq_profile_id
		Expression profile obtained from aligning the RNA-Seq data to the GRCm38.p3 reference genome.
	
	donor_id
		Donor from which the cells were obtained.
		
	donor_age
		Donor age.

	genotype_driver
		Transgenic mouse driver line.

	genotype_reporter
		Transgenic mouse reporter line.

	cell_reporter
		Detection status of reporter in cell.

	sampling_region
		Brain region targeted for cell sampling.
		
	cell_prep_sample_id
		Unique identifier of dissected brain sample.

	slice_id and slice_name
		Unique identifier and name of brain slice that was sampled.

	slice_index
		Brain slice number (from posterior to anterior) within each donor.

	cell_id
		Unique identifier of cell.
		
	rnaseq_profile_total_reads
		Total number of reads from RNA-sequencing.
		
	rnaseq_profile_percent_reads_aligned_to_mrna and rnaseq_profile_percent_reads_aligned_to_ncrna
		Percentage of total reads mapping to coding and non-coding RNA.

	rnaseq_profile_percent_reads_aligned_to_genome_only
		Percentage of total reads mapping to the genome but not transcriptome.
		
	subclass, subclass_order, and subclass_color
		Label that groups cells by putative type. Names are composed of marker genes for each subclass.

	broad_class, broad_class_order, broad_class_color
		Label that groups subclasses into broader cell types. Names are composed of marker genes for each broad class.
	
	
rows-genes.csv
	Contains information about the genes for which fpkm values were calculated. 

	gene_id
		Unique identifier for the gene.

	chromosome
		Chromosome associated with the gene.

	gene_entrez_id, gene_symbol, gene_name
		NCBI Entrez ID, gene symbol, and gene name.
