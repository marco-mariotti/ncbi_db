#! /usr/bin/python -u
from string import *
import sys
from commands import *
sys.path.append('/home/mmariotti/scripts')
from MMlib import *

help_msg="""This program parse a feature_table file download from NCBI and extract its features 
Usage:   $  ncbi_feature_table.py   [-option  value]

### Input:
-i   NCBI feature_table    *** required ***
-g   fasta genome file;    this is required if you request any sequence output file

## Possible output files:
-gff    write gff output to this file
-nt     write this nucleotide sequence fasta file for all features
-pep    write this protein sequence fasta file (only for CDS fields)

# Options:
-x      skip the lines with any of these comma-sep features (DEFAULT: "gene")
-f      select only these comma-sep features  (overrides -x)
-u      convert nt sequence output to uppercase
-c      check if protein sequences are correct, given this proteome in fasta format
-v      verbose; print explicitly every attribute for each entry in the file

-print_opt      print currently active options
-h OR --help    print this help and exit"""

command_line_synonyms={}

def_opt= {#'temp':'/home/mmariotti/temp', 
'i':0, 'g':0, 
'gff':0, 'pep':0, 'nt':0,  'o':0,
'c':0,'v':0,  'u':0,
'x':'gene', 'f':0,
}
## -o      provide an output prefix, so you can avoid specifying an argument for the other output options


#########################################################
###### start main program function

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms )
  else:  opt=args
  set_MMlib_var('opt', opt)
  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input
  input_file=opt['i'];   check_file_presence(input_file, 'input_feature_table_file')

  if opt['pep'] or opt['nt'] or opt['c']:
    genome_file=opt['g'];  check_file_presence(genome_file, 'input_genome_file')
    load_sequence_db(genome_file)

  if opt['c']:
    control_id2sequence={t.split()[0]:upper(s)   for t, s in parse_fasta(opt['c'])}

  if opt['gff']:  out_gff_fh=open(opt['gff'], 'w')
  if opt['pep']:  out_pep_fh=open(opt['pep'], 'w')
  if opt['nt']:   out_nt_fh= open(opt['nt'],  'w')

  exclude_fields={ f: None for f in opt['x'].split(',') }
  if opt['f']: accept_only={ f: None for f in opt['f'].split(',') }
  else:        accept_only=None

  locus_tag_named=0

  fields="# feature	class	assembly	assembly_unit	seq_type	chromosome	genomic_accession	start	end	strand	product_accession	non-redundant_refseq	related_accession	name	symbol	GeneID	locus_tag	feature_interval_length	product_length	attributes"[2:-1].split('\t')    ### this makes sense only if one wants to parse a partial feature table file
  #### opening input file now
  for line_index, line in enumerate(open(input_file)):
    if line_index==0 and line.startswith('#'):    fields=line[2:-1].split('\t'); continue
    splt=line[:-1].split('\t')
    if not splt: continue
    assert len(splt) == len(fields)
    x={ fields[i]: v for i, v in enumerate(splt) }

    # feature filter
    if accept_only:
      if not x['feature'] in accept_only:    continue
    elif     x['feature'] in exclude_fields: continue  

    if opt['v']: write('#'*60+'\n{0}'.format( join([str(f).ljust(25)+': '+str(x[f]) for f in x if x[f]], '\n')), 1)
    if not x['genomic_accession']: raise Exception, 'ERROR no chromosome/contig name in this object: '+str(x)

    the_strand=x['strand']
    if the_strand=='?': the_strand='+'
    g=gene(chromosome=x['genomic_accession'], strand=the_strand, tag=x['feature'])
    g.add_exon(int(x['start']), int(x['end'])  )
    g.id= x['product_accession']
    if not g.id:   
      if x['locus_tag']: g.id=x['locus_tag'];   locus_tag_named+=1
      else:              raise Exception, "ERROR cannot find a valid product_accession in this line: "+line
    g.title=g.id
    if x['name']:  g.title+=' ' + x['name']

    if opt['gff']:      print >> out_gff_fh, g.gff(program='annotation')

    if opt['pep'] or opt['nt'] or opt['c']:
      nt_seq=g.fast_sequence()
      if opt['u']:      nt_seq=upper(nt_seq)
      if opt['nt']:     print >> out_nt_fh,  ">{0}\n{1}".format(g.title, nt_seq)
      if x['feature']=='CDS' and (opt['c'] or opt['pep']): 
        cds_seq =  nt_seq[   :len(nt_seq)/3 *3  ] # some seqs are annotated not multiple of 3 and this is not considered in translation
        pep_seq=transl(cds_seq).rstrip('*')
        if opt['c']:    
          proteome_seq=None
          try:               proteome_seq=control_id2sequence[g.id]
          except KeyError:   write('sequenceId not found in proteome: {0}'.format(g.id), 1)
          if not proteome_seq is None and proteome_seq != pep_seq: 
            write('Sequences differing for {0}!\nProteome:  {1}\nFeatTable: {2}'.format(g.id, proteome_seq, pep_seq), 1)
            write('#'*60+'\n{0}'.format( join([str(f).ljust(25)+': '+str(x[f]) for f in x if x[f]], '\n')), 1)

        if opt['pep']:  print >> out_pep_fh, ">{0}\n{1}".format(g.title, pep_seq)
    #raw_input('...')

  if locus_tag_named:
    write('WARNING {0} entries were named with a locus_tag instead of a product_accession, since absent'.format(locus_tag_named), 1)
  ###############



#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:
    pass

  if 'log_file' in globals(): log_file.close()


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
    raise 
