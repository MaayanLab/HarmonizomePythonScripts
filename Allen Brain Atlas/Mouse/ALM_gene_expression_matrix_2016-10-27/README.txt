=================================
Mouse ALM
=================================

RNA sequencing data of single cells isolated from mouse anterior lateral motor cortical area (ALM).

The data set includes 5992 single cells collected from various cortical layers of ALM.

The sequencing results were aligned and aggregated at the gene level using the RSEM algorithm, and FPKM values were calculated.

For more details, please see the Documentation tab in the Cell Types web application.


Gene expression data matrix
	fpkm_table.csv
		Contains the (row, column) matrix of fpkm values obtained for each (gene, cell).
		The first row contains the unique identifiers of the RNA-seq profiles of the cells (rnaseq_profile_id)
		The first column contains the gene unique identifiers (gene_id)

		
Sample information (columns-cells.csv)
	rnaseq_profile_id
		Expression profile obtained from aligning the RNA-Seq data to the GRCm38.p3 reference genome.
	donor_id
		Donor from which the cells were obtained.
	organism
		Mouse cells for thi data set
	donor_age
		Donor age.
	genotype_driver
		Transgenic mouse driver line.
	genotype_reporter
		Transgenic mouse reporter line.
	cell_reporter
		Detection status of reporter in cell.
	cell_prep_sample_id
		Unique identifier of dissected brain sample.
	sampling_region
		Brain region targeted for cell sampling.
	sample_id
		Sample unique identifier
	sample_type
		Sample type
	facs_well
		FACS well
	facs_container
		FACS container unique identifier
	rna_amplification
		Amplificaiton well
	rna_amplification_set
		Amplificaiton plate
	amplified_quantity_ng
		Amplificaiton cDNA yield in ng
	library_prep
		Library well
	library_prep_set
		Library plate
	library_prep_avg_size_bp
		Average size of library in base pairs
	library_prep_quantification_fmol
		Library yield in fmol
	tube_sent_for_sequencing
		Library pool barcode
	batch_sent_for_sequencing
		Sequencing Batch
	rnaseq_profile_total_reads
		Total reads
	rnaseq_profile_percent_reads_aligned_to_mrna
		% reads aligned to mRNA
	rnaseq_profile_percent_reads_aligned_to_ncrna
		% reads aligned to non-coding RNA
	rnaseq_profile_percent_reads_aligned_to_genome_only
		% reads aligned to gDNA
	rnaseq_profile_percent_reads_aligned_to_synthetic_constructs
		% reads aligned to ERCC spike-in transcripts
	rnaseq_profile_percent_reads_aligned_to_ecoli
		% reads aligned to E. coli
	rnaseq_profile_percent_reads_not_aligned
		% unmapped reads
	fpkm>0_gene_count
		# of genes with fpkm values greater than 0
	fpkm>1_gene_count
		# of genes with fpkm values greater than 1
	fpkm>4_gene_count
		# of genes with fpkm values greater than 4
	fpkm>8_gene_count
		# of genes with fpkm values greater than 8
	fpkm>16_gene_count
		# of genes with fpkm values greater than 16
	fpkm>32_gene_count
		# of genes with fpkm values greater than 32
	fpkm>64_gene_count
		# of genes with fpkm values greater than 64

		
Gene information (rows-genes.csv)
	gene_id
		Unique identifier for the gene.
	chromosome
		Chromosome associated with the gene.
	gene_entrez_id
		NCBI Entrez ID
	gene_symbol
		Gene symbol
	gene_name
		Gene name
