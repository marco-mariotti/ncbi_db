#!/usr/bin/env python3
from string import *
import sys
from subprocess import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
from Bio import Entrez
from time import sleep
import  search_ncbi_entries
Entrez.email = "dpto.bioinformatics@crg.es"

help_msg="""Search or fetch entries (species or lineages) from the ncbi taxonomy database through internet and print information for each hit. It wraps functions from Bio->Entrez.

ncbi_taxonomy.py  [define taxids || text search]   [other options]

## Define taxids (alternative, mutually exclusive):
-i   +     define inputfile  with taxids, one per line
-I  +,+,.. define taxonomy taxids, comma separated. 

## or search database:
-S  "+"    search a single string into the database, print information for each hit

### Options:
# these options allow to print additional information. See examples at the bottom
-r         print rank 
-c         print genetic code id
-m         print mitochondrial genetic code
-n         print folder-like name (with chars masked as {chXX} and spaces replaced by underscores)

-L         legacy mode. Some taxids were merged into others by ncbi. When requesting these, normally the entry for the new taxid is reported and a warning printed to stderr. If this option is active, the old id is used in output instead, so that there would be a perfect correspondance between what was asked for and what is in output (expect for taxids not found at all)
-D         deep search: if no result if found with the string provided (-S), it tries to trim it to get a partial result. E.g.  Homo sapiens mdeuslls  -> it will ignore mdeuslls and return the entry for Homo sapiens. 
-s         print summary instead of full lineage (condensate usual output to have it smaller)
-a     +   max attempts numbers. Sometimes unstability of network results in failed attempts to reach ncbi. The program tries by default this number of times before giving up. There's a wait time of 1 sec between attempts. Set -a to 0 to never give up. Careful: if this number is high (or 0), it means it takes a long time before it realizes there's a problem, like: no internet connection.
-b         batch size. default is 200, but it is lowered by 15% if the following error is found: IOError("Requested URL too long (try using EPost?)") 
-v         verbose mode
-h OR --help    print this help and exit

## Output:
ID : lineage (; separated) -> Scientific name [ additional fields (# separated) ]

### Example:   ncbi_taxonomy.py -S "human" -r -c -m -n 
9606 : cellular organisms; Eukaryota; Fungi/Metazoa group; Metazoa; Eumetazoa; Bilateria; Coelomata; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Euarchontoglires; Primates; Haplorrhini; Simiiformes; Catarrhini; Hominoidea; Hominidae; Homininae; Homo -> Homo sapiens # Rank: species # GeneticCode: 1=Standard # MitoGeneticCode: 2=Vertebrate Mitochondrial # MaskedName: Homo_sapiens

If used as a function (imported in another python program), the main() function must be called with a opt dictionary as argument. opt['silent'] should be 1. An dictionary is returned. For each entry found, its taxid is a key in the dictionary, the corresponding value is the string that would be printed (so, possibly affected by options -r, -c, -m)
"""

command_line_synonyms={}

def_opt= { #'temp':'/users-d3/mmariotti/temp', 
'i':0, 'I':'', 'L':0,
'S':"", 'c':0, 'm':0, 'r':0, 's':0, 'D':0, 'n':0,
'v':0, 'Q':0,
'b':200, 'silent':0, 'a':100,
}



def lineage_string_to_abstract(lineage):
  """ lineage is a string which is usually returned by this program. This function condensate it keeping the most interesting classes. """
  splt=lineage.split('; ')
  if "Bacteria; " in lineage  :    return 'B; '+'; '.join(splt[2:min(5, len(splt))])
  elif 'Archaea; ' in lineage :    return 'A; '+'; '.join(splt[2:min(5, len(splt))])
  elif 'Eukaryota; ' in lineage:   
    out='E; '
    if 'Metazoa; ' in lineage:
      out+='M; '
      if 'Deuterostomia; ' in lineage:
        out+='Deuterostomia; '
        if 'Vertebrata; ' in lineage:
          out+='Vertebrata; '      
          if 'Mammalia; ' in lineage:           out+='Mammalia; '
          elif 'Sauropsida; ' in lineage:       out+='Sauropsida; '
          elif 'Amphibia; ' in lineage:         out+='Amphibia; '
          elif 'Actinopterygii; ' in lineage:   out+='Actinopterygii; '
          elif 'Chondrichthyes; ' in lineage:   out+='Chondrichthyes; '
        elif 'Tunicata; ' in lineage:          
          out+='Tunicata; '
          if 'Ascidiacea; ' in lineage:           out+='Ascidiacea; '         
        elif 'Branchiostomidae; ' in lineage:  out+='Branchiostomidae; '
        elif 'Echinodermata; ' in lineage:     out+='Echinodermata; '
      elif 'Protostomia; ' in lineage:
        out+='Protostomia; '
        if  'Arthropoda; ' in lineage:      
          out+='Arthropoda; '
          if 'Insecta; ' in lineage:        out+='Insecta; '
          elif 'Crustacea; ' in lineage:    out+='Crustacea; '
          elif 'Myrapoda; ' in lineage:     out+='Myrapoda; '
          elif 'Arachnida; ' in lineage:     out+='Arachnida; '
          elif 'Merostomata; ' in lineage:     out+='Merostomata; '
        elif 'Nematoda; ' in lineage:     out+='Nematoda; '
        elif 'Mollusca; ' in lineage:           
          out+='Mollusca; '
          if 'Gastropoda; ' in lineage: out+='Gastropoda; '
          elif 'Bivalvia; ' in lineage: out+='Bivalvia; '
        elif 'Annelida; ' in lineage:           out+='Annelida; '
        else:      out+=     lineage.split('Protostomia; ')[1].split(';')[0]+'; '
      else: #basal metazoan
        if 'Cnidaria; ' in lineage:        out+='Cnidaria; '
        elif 'Porifera; ' in lineage:      out+='Porifera; '          
        elif 'Ctenophora; ' in lineage:    out+='Ctenophora; '
        elif 'Placozoa; ' in lineage:      out+='Placozoa; '
        elif 'Platyhelminthes; ' in lineage: out+='Platyhelminthes; '

    else:      out+= '; '.join(splt[2:min(4, len(splt))])+'; '
    return out[:-2]
  else:      return '; '.join(splt[0:min(4, len(splt))])
            
  
#########################################################
###### start main program function

def verbose(msg):
  if opt['v']:      print(msg)

def id_is_included_in_list(main_id, id_list):
  for a_id in id_list:
    if main_id in a_id:
      return True
  return False


def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms )
  else:  
    opt=args
    for k in def_opt: 
      if k not in opt: opt[k]=def_opt[k]
    
  #determining input ids
  id_list=[]
  if opt['i']:
    check_file_presence(opt['i'], 'input_file_with_ids')
    for line in open(opt['i']):       id_list.append(line[:-1])
  elif opt['I']:    id_list=str(opt['I']).split(',')
  elif not opt['S']:  raise Exception("ERROR must provide ids either with option -i FILE or option -I id1,id2,...,idN, or a search string with opt -S.  Run with option -h for help")
  
  ########## search string
  if opt['S']:
    if "_"in opt['S'] and not ' ' in opt['S']:      opt['S']=replace_chars(opt['S'], '_', ' ')
    opt['S']=unmask_characters(opt['S'])
    search_ncbi_entries_opt={'silent':1, 's':opt['S'], 'm':'T', 'a':opt['a']}
    id_list= search_ncbi_entries.main(search_ncbi_entries_opt)    
      
    # searching again changing characters such as "_", ":"
    if not id_list:
      search_ncbi_entries_opt['s']=replace_chars(opt['S'], '_', ' ')
      id_list= search_ncbi_entries.main(search_ncbi_entries_opt)
    if not id_list:
      search_ncbi_entries_opt['s']=replace_chars(opt['S'], ':', ' ')
      id_list= search_ncbi_entries.main(search_ncbi_entries_opt)
    if not id_list:
      search_ncbi_entries_opt['s']=replace_chars(opt['S'], '()', ' ')
      id_list= search_ncbi_entries.main(search_ncbi_entries_opt)
    
    if opt['D'] and not id_list:
      # approximate match      
      search_string= opt['S'];       search_string= ' '.join(search_string.split()[:len(search_string.split())-1])
      while not id_list and search_string:
        search_ncbi_entries_opt={'silent':1, 's':search_string, 'm':'T', 'a':opt['a']}
        id_list= search_ncbi_entries.main(search_ncbi_entries_opt)
        if not id_list: search_string= ' '.join(search_string.split()[:len(search_string.split())-1])        

      if not search_string:
        # approximate match  replacing _ with " "
        search_string= replace_chars(opt['S'], '_', ' ');       search_string= ' '.join(search_string.split()[:len(search_string.split())-1])
        while not id_list and search_string:
          search_ncbi_entries_opt={'silent':1, 's':search_string, 'm':'T', 'a':opt['a']}
          id_list= search_ncbi_entries.main(search_ncbi_entries_opt)
          if not id_list: search_string= ' '.join(search_string.split()[:len(search_string.split())-1])        

      if id_list and search_string: printerr('WARNING returning the approximate match reducing: "'+opt['S']+'"  to:  "'+search_string+'"', 1)
      
  verbose('tot hits: '+str(len(id_list)))
  #batch preparation
  batch_size=opt['b']
  batch_index=0
  id_retrieved_until_now=0

  output_hash={}

  while id_retrieved_until_now!=len(id_list):
    try:
      batch_id_list=id_list[id_retrieved_until_now:id_retrieved_until_now+batch_size]
      if len(batch_id_list)!=len(id_list) and not opt['silent']:
        verbose('Batch n'+str(batch_index+1)+' --- entries:'+str(len(batch_id_list)))
      
      number_of_records=0
      batch_id_list_index=0

      #### fetching!
      fetch_successful=False
      attempts=0
      while not fetch_successful and (  attempts < opt['a'] or opt['a']==0):
        try:
          fetch_handle=Entrez.efetch(db='taxonomy', id=','.join(batch_id_list), retmode='xml')    ## connecting to ncbi via internet
          parsed_results=[r for r in Entrez.parse(fetch_handle)]
          fetch_successful=True
        except KeyboardInterrupt:          raise
        except:
          service('Something failed (network problem?). Trying again... attempt number: '+str(attempts)+ '(will keep indefinitely'*int(opt['a']==0)+int(opt['a']!=0)*' (will stop at attempt n.'+str(opt['a'])+')')
          sleep(1)
          attempts+=1

      if not fetch_successful:
        printerr("ERROR couldn't fetch results. Number of attempts: "+str(attempts)+'. (Last) exception was: ', 1)
        raise

      ###### parsing fetch results
      for record in parsed_results:
        ### note: batch_id_list[batch_id_list_index]   is the taxid we asked for
        out_taxid= record['TaxId']   #string

        while batch_id_list_index < len(batch_id_list) and batch_id_list[batch_id_list_index] != out_taxid :  #### entries should be in the same order in which we asked for. thus we check sequentially if taxid correspond. we enter here if it does not
          if 'AkaTaxIds' in record  and  batch_id_list[batch_id_list_index] in record['AkaTaxIds']:
            #### the id requested was merged into something else, this is why it doesn't match
            if opt['L']: 
              if not opt['silent']:    printerr("WARNING! Reporting entry " +batch_id_list[batch_id_list_index]+" instead of newer: " + out_taxid+' since option -L is active (species: '+record['ScientificName']+')', 1)
              out_taxid = batch_id_list[batch_id_list_index]  ### legacy mode, using the taxid you requested even if it's not the most recent one at ncbi
            else:  
              if not opt['silent']:    printerr("WARNING! Entry " +batch_id_list[batch_id_list_index]+" was merged into this: " +  out_taxid+' | '+record['ScientificName'], 1)
              batch_id_list[batch_id_list_index]=out_taxid
          else:
            ##### id was not found!!
            if not opt['silent']:  printerr("WARNING! Entry " +batch_id_list[batch_id_list_index]+" was not found and was ignored", 1)
            #printerr(    "size batch_id_list ", len(batch_id_list),  'batch_id_list_index', batch_id_list_index,   1)
            output_hash [batch_id_list[batch_id_list_index]]=''
            batch_id_list_index+=1

        
        if batch_id_list[batch_id_list_index]==out_taxid:
          lineage=record['Lineage']
          if opt['s']:   lineage=lineage_string_to_abstract(lineage)
          out=out_taxid+' : '+lineage+' -> '+record['ScientificName']

          #adding additional fields
          if opt['r']:           out+= ' # Rank: '+str(record['Rank'])
          if opt['c']:           out+= ' # GeneticCode: '+str(record['GeneticCode']['GCId'])+'='+record['GeneticCode']['GCName']
          if opt['m']:           out+= ' # MitoGeneticCode: '+str(record['MitoGeneticCode']['MGCId'])+'='+record['MitoGeneticCode']['MGCName']
          if opt['n']:           out+= ' # MaskedName: '+ mask_characters(record['ScientificName']).replace( ' ', '_' )

        
          if not opt['silent']: print(out)
          output_hash [batch_id_list[batch_id_list_index]]=out

        #if batch_id_list_index == len(batch_id_list): return output_hash
        batch_id_list_index+=1
        number_of_records+=1

      if batch_id_list_index<len(batch_id_list): 
        for i in range(batch_id_list_index,len(batch_id_list)):  
          if not opt['silent']: printerr("WARNING! entry " +batch_id_list[i]+" was not found and was ignored", 1)
          output_hash [batch_id_list[batch_id_list_index]]=''
          
      id_retrieved_until_now+=len(batch_id_list)
      batch_index+=1
    except IOError as e:
      if 'Requested URL too long' in e.args[0]:
        batch_size=int(batch_size*0.85)
        if not opt['silent']: printerr('requested URL too long... trying to decrease batch size to:'+str(batch_size), 1)


  return output_hash
#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:
    pass


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
    raise 
