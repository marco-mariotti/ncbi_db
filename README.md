# ncbi_db
Collection of commands to query or process NCBI data

## Installation
conda install -c mmariotti -c conda-forge ncbi_db

## Tools
These tools are available:

- **ncbi_assembly**       search and download assemblies/genomes for any species/lineage, or its annotation/proteome
- **ncbi_sequences**      search and download nucleotide/protein sequences or their metadata
- **ncbi_pubmed**         search and format ncbi pubmed entries
- **ncbi_taxonomy**       search ncbi taxonomy for species or lineages
- **ncbi_taxonomy_tree**  obtain a tree from ncbi taxonomy for a set of input species
- **ncbi_search**         generic search tool for any ncbi DB
- **parse_genbank**       parse a genbank flat file; requires installation of GBParsy


Run any tool with option -h to display its usage.

Most tools require internet, as they connect online to ncbi.

