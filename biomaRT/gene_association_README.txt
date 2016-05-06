

gene_association.cgd.gz    Contains all GO annotations for CGD (protein and RNA)

The gene_association.cgd.gz file uses the standard file format for
gene_association files of the Gene Ontology (GO) Consortium.  A more
complete description of the file format is found here:

http://www.geneontology.org/doc/GO.annotation.html#file

Columns are:					Contents:

 1) DB						- database contributing the file (always "CGD" for this file)
 2) DB_Object_ID				- CGDID
 3) DB_Object_Symbol				- see below
 4) Qualifier 			(optional)	- 'NOT', 'contributes_to', or 'colocalizes_with' qualifier 
							 for a GO annotation, when needed
 5) GO ID					- unique numeric identifier for the GO term
 6) DB:Reference(|DB:Reference)			- the reference associated with the GO annotation
 7) Evidence					- the evidence code for the GO annotation
 8) With (or) From 		(optional)	- any With or From qualifier for the GO annotation
 9) Aspect					- which ontology the GO term belongs in (see note below)
10) DB_Object_Name(|Name) 	(optional)	- a name for the gene product in words, e.g. 'acid phosphatase'
11) DB_Object_Synonym(|Synonym) (optional)	- see below
12) DB_Object_Type				- type of object annotated, e.g. gene, protein, etc.
13) taxon(|taxon)				- taxonomic identifier of species encoding gene product
14) Date					- date GO annotation was made
15) Assigned_by					- source of the annotation


Notes on particular columns:

Column 3 - When a Standard Gene Name (e.g. C. albicans CDC28 or COX2) has been
        conferred, it will be present in Column 3. When no Gene Name
        has been conferred, the ORF Name (e.g., C. albicans orf19.6632
        or C. glabrata CAGL0K12694g) will be present in column 3.

Column 9 - Aspect
	  C = Cellular Component
	  F = Molecular Function
	  P = Biological Process

Column 11 - The ORF Name (e.g., orf19.6632, CAGL0K12694g) will be the first name present 
	in Column 11. Any other names (except the Standard Name, which will 
	be in Column 3 if one exists), including Aliases used for the gene, 
	will also be present in this column.

Column 15 - GO annotations in CGD are either assigned by CGD curators
       or assigned computationally (see
       http://www.candidagenome.org/cgi-bin/reference/reference.pl?dbid=CAL0121033 and
       http://www.candidagenome.org/cgi-bin/reference/reference.pl?dbid=CAL0142013).
       Previously, CGD also included IEA annotations based on
       predictions made by the Annotation Working Group (see Braun BR, et
       al. (2005) A human-curated annotation of the Candida albicans
       genome. PLoS Genet 1(1):36-57. PMID: 16103911).  Many of these
       predictions were replaced with literature-based annotations during reference
       curation at CGD, and the remainder were archived on 10 Nov 2010.

Note:

This file contains ALL of the GO curation at CGD, whereas the
gene_association file that is available on the GO consortium (GOC) web
site, http://www.geneontology.org/, has been filtered according to GOC
guidelines, which are discussed in more detail at
http://www.geneontology.org/GO.annotation.shtml.
Before the Annotation Working Group annotations were archived in
November 2010 (see above), these annotations were being filtered from
the GOC version of the file because they had an IEA evidence code and
were over a year old.

Note:

There was an error in the "Qualifier" column of the gene_association
file between April 23, 2008 and October 14, 2008. The 'contributes_to'
and 'colocalizes_with' qualifiers were noted incorrectly as the 'NOT' qualifier
in this file between those dates. This error has been corrected as of October 15, 2008.

Note:

The files are gzip compressed tab-delimited text files. There are
several freely available software options for decompressing
gzipped files using Windows.  The software and other useful
information is available on these web sites:

- WinZip (http://www.winzip.com/)
- Stuffit (http://www.stuffit.com/)
- Gzip (http://www.gzip.org/

and the gzip user's manual:
http://www.math.utah.edu/docs/info/gzip_toc.html


