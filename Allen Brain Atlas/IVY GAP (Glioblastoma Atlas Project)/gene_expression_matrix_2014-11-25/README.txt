======================
Ivy Glioblastoma Atlas 
======================

RNA sequencing data for anatomic structures and putative cancer stem cell clusters isolated by laser microdissection. 

The data set includes 122 RNA samples of 5 anatomic structures (Leading Edge, Infiltrating Tumor, Cellular Tumor, Microvascular Proliferation, and Pseudopalisading Cells Around Necrosis) identified by H&E staining in 10 tumors, and 148 RNA samples of putative cancer stem cell clusters identified by ISH with 18 probe reference set in 34 tumors.

The sequencing results were aligned and aggregated at the gene level using the RSEM algorithm, and the resulting fpkm values were normalized across all samples based on genes not enriched in particular anatomic structures. 

For more details, please see the Downloads tab in the Ivy Glioblastoma Atlas web application.


fpkm_table.csv
	Contains the (row, column) matrix of fpkm values obtained for each (gene, sample).
	The first row contains the sample unique identifiers (rna_well_id)
	The first column contains the gene unique identifiers (gene_id)


columns-samples.csv 
	Contains information about the samples profiled with RNA sequencing
		
	tumor_id, block_id, specimen_id and corresponding names
		Specimen from which the sample was dissected and the specimen's parent block and tumor.

	rna_well_id
		Unique identifier of the sample.

	polygon_id
		Unique identifier of an avg_graphic_object that outlines where the sample was cut from.
			

	structure_id, structure_abbreviation, structure_color, structure_name
		Label that groups samples obtained from laser micro-dissected anatomic structures or putative cancer stem cell clusters.


rows-genes.csv
	Contains information about the genes for which fpkm values were calculated. 

	gene_id
		Unique identifier for the gene.

	chromosome
		Chromosome associated with the gene.

	gene_entrez_id, gene_symbol, gene_name
		entrez_id, NCBI symbol, and name of the gene.


