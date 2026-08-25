The GFF3 files in this directory are the same as those downloaded from Ensembl, with additional 
column 9 attributes appended to some gene-level feature rows. Specifically, where they can be inferred,
we have added "projection_parent_gene", which specifies the Ensembl gene model id for the corresponding gene in 
the reference genome (C57BL/6J), "mgp" which specifies the original Mouse Genomes Project ids from earlier releases.

Backgound/motivation: Prior to Ensembl release-114, the mouse strain genome annotations (GFF3 files) contained
projection_parent_gene attributes, which we used to link strain gene models with the canonical gene from MGI and 
with corresponding gene models in other strains. For example, the A/J gene ID=MGP_AJ_G0026191, had 
projection_parent_gene=ENSMUSG00000027168, which we use to link this gene with MGI:97490 (Pax6) and with equivalent
genes in other strains, e.g., ID=MGP_CBAJ_G0025928 in CBA/J which also has projection_parent_gene=ENSMUSG00000027168.
Starting with version 114, there were two changes that affected (i.e. broke) a number of our pipelines.
(1) projection_parent_gene is missing in all protein coding genes. They are present for non-coding genes, 
but not protein coding.
We asked several times to have them restored, but to no avail. 
(2) Gene models were assgned all new IDs. A/J gene MGP_AJ_G0026191 is now ENSMUSG00195028435. 
No mapping is provided between the old and new IDs.

The goal of patching is to infer (where possible) projection_parent_gene and add it to the record. 
We can then feed the patched files to our existing downstream pipelines. At the same time, we are often
able to recover the prior MGP ids and when successful, we add them as a new attribute "mgp". 

The patched files are being made public in case others find them useful.

Example: patched line for A/J Pax6 gene:
2	ensembl	gene	102484499	102513003	.	+	.	ID=gene:ENSMUSG00195028435;Name=Pax6;biotype=protein_coding;description=paired box 6 [Source:NCBI gene (formerly Entrezgene)%3BAcc:18508];gene_id=ENSMUSG00195028435;logic_name=ensembl;version=1;projection_parent_gene=ENSMUSG00000027168;mgp=MGP_AJ_G0026191

Complication:
In the majority of cases, the associations between strain genes, projection parents, and MGP ids is nicely 1:1.
But not aways. Multiple strain genes in a given genome may have the same projection parent. This was true before 
release 114 and is true after. The bottom line is that both projection_parent_gene and mgp attributes can be multivalued.

Example: gene Pakap (MGI:5141924). Prior to release-114, there were three A/J genes with the same associated with 
Pakap, and in the current release, there are two. The patched lines from the current release:
4	ensembl	gene	53890198	54168706	.	+	.	ID=gene:ENSMUSG00195007517;Name=Pakap;biotype=protein_coding;description=paralemmin A kinase anchor protein [Source:NCBI gene (formerly Entrezgene)%3BAcc:677884];gene_id=ENSMUSG00195007517;logic_name=ensembl;version=1;projection_parent_gene=ENSMUSG00000038729,ENSMUSG00000089945,ENSMUSG00000090053;mgp=MGP_AJ_G0028278,MGP_AJ_G0028279,MGP_AJ_G0028280
4	ensembl	gene	54198792	54354284	.	+	.	ID=gene:ENSMUSG00195007527;Name=Pakap;biotype=protein_coding;description=paralemmin A kinase anchor protein [Source:NCBI gene (formerly Entrezgene)%3BAcc:677884];gene_id=ENSMUSG00195007527;logic_name=ensembl;version=1;projection_parent_gene=ENSMUSG00000038729,ENSMUSG00000089945,ENSMUSG00000090053;mgp=MGP_AJ_G0028278,MGP_AJ_G0028279,MGP_AJ_G0028280


More details:

1. Patching the files (adding "projection_parent_gene=").
Patching uses other information often included in the gene's col 9 attributes:
- NCBI Gene ID: The feature's description attribute may contain a substring like this:
  "[Source:NCBI gene (formerly Entrezgene);Acc:226304]"
- MGI ID: The description may contain a substring like this:
  "[Source:MGI Symbol;Acc:MGI:3801960]"
- Mouse symbol: The feature may have a Name attribute containing a mouse symbol.

2. Limiting association counts.
Another (separate) issue with v116 data is that there are cases of multiple strain genes (sometimes many) 
within a single genome that refer to the same projection parent gene. The sequence group discussed this and 
decided to cap (at 3) the max number of strain genes (in one genome) that can be associated with 
one projection parent. If there are more than the limit, then none of the strain gene is associated with 
the projection parent. 

3. Old MGP ids (adding "mgp=").
The patching process uses a snapshot of data from MGI to recover the MGPs for many genes. 
The archived data is a table, one row per strain gene, with columns that include the
MGP id and the Ensembl id (projection parent) of the associated MGI gene. The Ensembl id is 
used to match against the projection_parent_genes (either supplied or inferred in part 1). 
Matched MGP ids are added as a new attribute ("mgp") in column 9.

