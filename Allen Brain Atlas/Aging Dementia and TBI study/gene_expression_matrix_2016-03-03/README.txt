==============================
Aging, Dementia, and TBI study 
==============================

RNA sequencing data of hippcampus and cortex in an aged cohort isolated by macrodissection. 

The data set includes 377 RNA-Seq samples collected from hippocampus, temporal cortex, and parietal cortex (both grey and white matter) in 55 aged donors with TBI and their matched controls (107 donors total after QC).  Additional donor meta-data, neuropathology metrics, and IHC image quantifications for each sample are available.    

The sequencing results were aligned and aggregated at the gene level using the RSEM algorithm, and the resulting fpkm values were normalized across all samples within each brain region to account for processing batch and RNA quality (RIN).  Data matrices for fpkm values before and after normalization are included here.

For more details, please see the Documentation tab in the Aging, Dementia, and TBI study web application.


fpkm_table_unnormalized.csv
	Contains the (row, column) matrix of fpkm values obtained for each (gene, sample) from the RSEM analysis pipeline.
	The first row contains the unique identifiers of the RNA-seq profiles of the samples (rnaseq_profile_id)
	The first column contains the gene unique identifiers (gene_id)

	
fpkm_table_normalized.csv
	Contains the (row, column) matrix of fpkm values obtained for each (gene, sample) after correcting for RIN and batch effects. These are the data displayed in the RNA-Seq page heatmap
	The first row contains the unique identifiers of the RNA-seq profiles of the samples (rnaseq_profile_id)
	The first column contains the gene unique identifiers (gene_id)

	
columns-samples.csv 
	Contains information about the samples profiled with RNA sequencing
	
	rnaseq_profile_id
		Expression profile obtained from aligning the RNA-Seq data to the GRCh38.p2 reference genome.
	
	donor_id and donor_name
		Donor from which the sample was dissected
		
	specimen_id and specimen_name
		Specimen from which the sample was dissected (i.e., a particular brain structure from a particular donor)
		
	rna_well_id
		Unique identifier of the sample.

	polygon_id
		Unique identifier of an avg_graphic_object that outlines where the sample was cut from.

	structure_id, structure_abbreviation, structure_color, structure_name
		Label that groups samples by brain region (hippocampus, temporal cortex, parietal cortex, and forebrain white matter).

	hemisphere
		Hemisphere from which the processed sample was collected


rows-genes.csv
	Contains information about the genes for which fpkm values were calculated. 

	gene_id
		Unique identifier for the gene.

	chromosome
		Chromosome associated with the gene.

	gene_entrez_id, gene_symbol, gene_name
		entrez_id, NCBI symbol, and name of the gene.
