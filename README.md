# ncbi_db
Collection of commands to query or process NCBI data

## Installation
conda install -c mmariotti -c conda-forge ncbi_db

## Tools
These command line tools are available:

- **ncbi_assembly**       search and download assemblies/genomes for any species/lineage, or its annotation/proteome
- **ncbi_sequences**      search and download nucleotide/protein sequences or their metadata
- **ncbi_pubmed**         search and format ncbi pubmed entries
- **ncbi_taxonomy**       search ncbi taxonomy for species or lineages
- **ncbi_taxonomy_tree**  obtain a tree from ncbi taxonomy for a set of input species
- **ncbi_search**         generic search tool for any ncbi DB
- **parse_genbank**       parse a genbank flat file; requires installation of GBParsy

Run any tool with option -h to display its usage.

Most tools require internet, as they connect online to ncbi.

## Using ncbi_db as module
To use these functionalities from another python module, import them from ncbi_db and
run their "main" function providing the same arguments as you would on the command line,
but in form of dictionary. Use option 'silent' to avoid printing results on screen. For example:

```
from ncbi_db import ncbi_sequences
arguments={'m':'P', 'f':1, 'I':'AAB88790', 'silent':1}
results=ncbi_sequences.main(arguments)
print(results)
```

```
{'AAB88790':
 ['AAB88790.1 gi|2411487|gb|AAB88790.1| selenophosphate synthetase [Drosophila melanogaster]',
 'MSYAADVLNSAHLELHGGGDAELRRPFDPTAHDLDASFRLTRFADLKGRGCKVPQDVLSKLVSALQQDYSAQDQEPQFLNVAIPRIGIGLDCSVIPLRHGGLCLVQTTDFFYPIVDDPYMMGKIACANVLSDLYAMGVTDCDNMLMLLAVSTKMTEKERDVVIPLIMRGFKDSALEAGTTVTGGQSVVNPWCTIGGVASTICQPNEYIVPDNAVVGDVLVLTKPLGTQVAVNAHQWIDQPERWNRIKLVVSEKNVRKAYHRAMNSMARLNRVAARLMHKYNAHGATDITGFGLLGHAQTLAAHQKKDVSFVIHNLPVIAKMAAVAKACGNMFQLLQGHSAETSGGLLICLPREQAAAYCKDIEKQEGYQAWIIGIVEKGNKTARIIDKPRVIEVPAKD']}
```

## Developers
Marco Mariotti https://github.com/marco-mariotti

Didac Santesmasses https://github.com/didacs
