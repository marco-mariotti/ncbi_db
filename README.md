# ncbi_db
Collection of scripts that wrap the Entrez Bio python module.
Mainly to build your own copy of selected entries form the ncbi databases

Right now these tools are available: 

- ncbi_lib.py          library of functions; not to be executed
- ncbi_db_info.py      simply to display the available fields in any ncbi DB, browsable through the BioPython tools
- ncbi_assembly.py     to browse the available assemblies/genomes for any species or lineage, including their annotation/proteome, and possibly download them in a local folder
- parse_genbank.py     parse a gbf file to extract selected entries. Requires GBParsy module
- ncbi_pubmed.py       wrapper to ncbi pubmed search engine; it can format results for bibliography
