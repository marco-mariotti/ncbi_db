#! /usr/bin/env python3
from string import *
import sys
from subprocess import *
from .MMlib3 import *
#import gbparsy # import GBParsy module


help_msg="""Extract specific feature from a GenBank flat file, using the GBParsy module.

Usage: $ parse_genbank.py   gbff_file    [other_options]   > output_file.fa

## Options:
-f   feature to extract (default: "gene")
-q   qualifiers to add to the fasta header  (special values accepted: "all", "fasta_title", "header")
-s   filter; word to search in values of any qualifier, only if found the entry will be printed
-t   try and use this qualifier as title; comma-separated list can be provided (default: auto; which means gene->locus_tag; cds->protein_id,locus_tag)
-p   for coding sequences extraction, output peptide sequence instead of nucleotide (implies -f CDS and -q product,locus_tag)
-l   long output; qualifiers are added in the format "qualifier=value;" instead of just "value"

-x   output 'source' entries instead of features (i.e. contig sequences); disregard all other options
-d   do not output sequence; the headers will also lack the leading ">"
-g   suppress normal output; output gff format instead  

-print_opt      print currently active options
-h OR --help    print this help and exit

#### Examples:
Extract proteome:  # parse_genbank.py genbank.gbff -p      > proteome.fa
Extract genome:    # parse_genbank.py genbank.gbff -x      > genome.fa
Extract mRNA:      # parse_genbank.py genbank.gbff -f mRNA > mRNAs.fa"""

command_line_synonyms={}

def_opt= { 
'i':0, 'q':0, 's':0, 'p':0, 'd':0, 'l':0,
'f':'gene', 't':'auto', 
'x':0, 'g':0,
}


#########################################################
###### start main program function

def main(args={}):
#########################################################
############ loading options
  import gbparsy # import GBParsy module
  

  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms )
  else:  opt=args
  set_MMlib_var('opt', opt)
  #checking input

  
  check_file_presence(opt['i'], 'Genbank gbff input file')
  sFileName = opt['i'];  sFeature = opt['f'];  
  sQualifiers = []  if not opt['q'] else opt['q'].split(',')
  sWord =      None if not opt['s'] else opt['s']
  add_fields = []   if not opt['a'] else opt['a'].split(',')
  title_fields= [] if not opt['t'] else opt['t'].split(',')
  get_protein= opt['p']
  
  if get_protein: 
    sFeature = "CDS"
    if not sQualifiers: sQualifiers=['product', 'locus_tag']

  if sQualifiers==['all']: sQualifiers='ALL'

  if title_fields ==['auto']:
    if    sFeature=='CDS':  title_fields=['protein_id','locus_tag']
    else  :                 title_fields=['locus_tag']

  

  id_dictionary={} # just to keep track of ids already used
  ids_changed_since_duplicated=0

  oGBParsy = gbparsy.gbparsy()
  service('parsing {f}'.format(f=sFileName))
  oGBParsy.parse(sFileName) # parse a GBF file which contains more than one GBF sequence data

#  printerr('something', 1)

  # If you want to get a list of SeqRecord instances, use getBioData method instead of getRawData.
  lSeqData = oGBParsy.getRawData()

  failed_parse_titles=0   
  lack_translation_warnings=0
  for dSeqData in lSeqData: # dSeqData has a parsed data of a GBF sequence data
      # start of user process
      if opt['x']:
        # for k in dSeqData:
        #   if k!='sequence':
        #     write( k, how='red')
        #     write(dSeqData[k], 1)
        # raw_input('..')
        
        accession = dSeqData['version']  if 'version' in dSeqData and dSeqData['version'] else    dSeqData['accession'] 
        source_title='{id} {desc}'.format(id=accession, desc=dSeqData['definition'])   ### you also have organism, lineage ...
        write(">{tit}\n{seq}".format(tit=source_title, seq=fasta(   dSeqData['sequence'].upper())  ),1  )
        continue

      for dFeature in dSeqData["features"]:
          if dFeature["feature"] == sFeature or sFeature=='*':
              #write(dFeature, 1)

              #for k in dSeqData.keys():
                #write( k+'   --> '+str(dSeqData[k]), 1, how='red')
                #raw_input('...')

              print_this=True if not sWord else bool( searchQualValue(sWord, dFeature) )  ## we print only if we find this word as value of any qualifiers
              title='{feat}_{start}-{end}_{strand}'.format(feat=dFeature["feature"], start=dFeature["start"], end=dFeature["end"], strand=dFeature["direction"])

              #printerr(title, 1)
              chromosome_name= dSeqData["version"] if 'version' in dSeqData else dSeqData["accession"]
              g=gene(chromosome=chromosome_name   , strand={'C':'-', 'N':'+'}[dFeature["direction"]])
              for start, end in dFeature["locations"]:          g.add_exon(start, end)
              
              for title_field in title_fields:
                q=getQualValue(title_field, dFeature)
                if not q:
                  if title_field == title_fields[-1]:             failed_parse_titles+=1
                else:                 title=q; break

              added_suffix=0
              while title in id_dictionary:
                if not added_suffix:  title=title+'-'+str(added_suffix+1)
                else:                 title=   '-'.join(title.split('-')[:-1]) +'-'+str(added_suffix+1)
                added_suffix+=1

              id_dictionary[title]=True
              if added_suffix: ids_changed_since_duplicated+=1

              g.id=title

              description=' '
              if sQualifiers != 'ALL': 
                for sQualifier in sQualifiers:                    
                  if   sQualifier == 'header':      value=g.header(no_id=True)
                  elif sQualifier == 'fasta_title': value=g.fasta_title()
                  else:                             value= getQualValue(sQualifier, dFeature)
                  if not opt['l']: description+=value+' '
                  else:            description+=sQualifier+'='+value+'; '
              else: 
                for sQualifier, value in dFeature["qualifiers"]:
                  if not opt['l']: description+=value+' '
                  else:            description+=sQualifier+'='+value+'; '

              if   opt['d']:
                write('{title}{desc}'.format(title=title, desc=description.rstrip()), 1)
              elif opt['g']:
                write(g.gff(program='Genbank', tag=dFeature["feature"], comment=description.rstrip()), 1)
              else:
                if get_protein: 
                  sequence=getQualValue('translation',  dFeature)
                  if not sequence: 
                    lack_translation_warnings+=1
                    if lack_translation_warnings<10:
                      printerr('WARNING lacking translation field for entry: {t} ; obtaining from raw sequence'.format(t=title), 1)
                    if lack_translation_warnings==10:
                      printerr('LAST WARNING: lacking translation field for entry: {t}; obtaining from raw sequence. Next warnings of this type will not be printed'.format(t=title), 1)                      
                    codon_start=getQualValue('codon_start',  dFeature)
                    cds=getSequence(dSeqData["sequence"], dFeature)
                    if codon_start and codon_start!='1':  cds=cds[int(codon_start)-1:]
                    sequence=transl( cds ).rstrip('*').rstrip('X')
                else:  sequence=getSequence(dSeqData["sequence"], dFeature) 
                write('>{title}{desc}\n{seq}'.format(title=title, desc=description.rstrip(), seq=sequence), 1)


  if failed_parse_titles: 
    printerr('WARNING: could not get desired field title for {n} entries; used a generic title instead'.format(n=failed_parse_titles), 1)

  if ids_changed_since_duplicated: 
    printerr('WARNING: {n} titles were not unique, and were changed'.format(n=ids_changed_since_duplicated), 1)

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



sNorBase = "ACGTRYMKWSBDHVNacgtrymkwsbdhvn";
sComBase = "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn";

def getRevCom(sSequence):
    lSequence = [sComBase[sNorBase.find(x)] for x in sSequence]
    lSequence.reverse()
    return "".join(lSequence)

def getQualValue(sQualifier, dFeature):
    if not sQualifier: return ""
    for tQualifier in dFeature["qualifiers"]:
        if tQualifier[0] == sQualifier: return tQualifier[1]
    return ""

def searchQualValue(sWord, dFeature):
    if not sWord: return ""
    for tQualifier in dFeature["qualifiers"]:
        if sWord in tQualifier[1]: return tQualifier[1]
    return ""

def getSequence(sWholeSequence, dFeature):
    sSequence = sWholeSequence[dFeature["start"] - 1: dFeature["end"]]
    if dFeature["direction"] == "C":        sSequence = getRevCom(sSequence)
    return sSequence


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
    raise 
