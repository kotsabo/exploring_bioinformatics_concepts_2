# Bioinformatics 1 - Assignment 2

## Tasks

### Q1 
Use a BLAST query to find potential orthologues of your gene. Describe what this query returned, and how you 
selected potential orthologues. Compare your result with a query for homologues, to find these, select 
'HomoloGene' in NCBI instead of 'Gene' as search option. How well do the results match when you consider 
e-value, % overlap or score from your BLAST search? List those results from your initial query that you think 
are true orthologues, and explain why. Note that NCBI BLAST has a useful view called 'Taxonomy Reports'. 
In addition to the NCBI service, you can also try the Ensembl database at 
http://www.ensembl.org/Homo_sapiens/Tools/Blast?db=core.

### Q2
Now consider the phylogenetic implications of your results. First, create a BLAST tree (using the function 
'Distance tree of results'). Describe the result and its relevance. How can this be used to study the disease 
in another model organism? Are there known orthologues in mouse, fruit fly or yeast? Using the information under 
'General gene information', discuss whether mouse, fruit fly or yeast (or any of these for which an orthologue 
exists) are suitable models to study the disease. Note: To give an example, if the gene is implicated in cell 
cycle, yeast may well be a good model because it is easier to study than mice. But if the gene is relevant for
brain function, yeast may of course be less relevant even if an orthologue exists.

### Q3
Select four orthologues, and create a rooted phylogenetic tree using the UPGMA algorithm as shown in class. 
Does it match the BLAST result from Q1/2?