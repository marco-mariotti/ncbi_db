#!/usr/bin/env python3
from string import *
import sys
from subprocess import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
from Bio import Entrez
from time import sleep

Entrez.email = "dpto.bioinformatics@crg.es"

help_msg="""Fetches entries from one of the NCBI databases (nucleotide, protein, taxonomy) through a "all names" search. It wraps functions from Bio->Entrez. Useful to get to the gi code or taxonomy id from any other code (or species name).

search_ncbi_entries.py "search query string" [options]

## Define mode (ncbi database):
-m  +      Possibilities: N -> nucleotide, P -> protein, T -> taxonomy

## Define output:  by default, the id of results are printed, one per line. if any of these options is active, this output is not produced.
-c         print a count of results

### Options:
-a     +                  max attempts numbers. Sometimes unstability of network results in failed attempts to reach ncbi. The program tries by default this number of times before giving up. Set it to 0 to never give up. Careful: if this number is high (or 0), it means it takes a long time before it realizes there's a problem, like: no internet connection.
-print_opt                print currently active options
-h OR --help              print this help and exit

### used as function (imported from another python program), the main function must be run using the opt hash which you would use to get the information, and using the opt['silent']=1.
An output_list is returned, containing the ids of entries matching the search.
For tools more specialized than this, see fetch_ncbi_sequences.py and ncbi_taxonomy.py
"""

command_line_synonyms={}

def_opt= {'temp':'/users-d3/mmariotti/temp', 
's':'',
'm':'P',
'v':0, 'Q':0, 'a':100,
'c':0, 'silent':0,
'M':10000
}

mode_synonyms={'N':'nucleotide', 'P':'Protein', 'T':'Taxonomy'}

#########################################################
###### start main program function

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 's', synonyms=command_line_synonyms )
  else:  
    opt=args
    for k in def_opt: 
      if k not in opt: opt[k]=def_opt[k]
  #global temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  output_list=[]
  
  #mode
  try:    db=mode_synonyms[opt['m'].upper()] 
  except: raise Exception("ERROR mode not recognized: "+opt['m'])
  #determining input ids

  if not opt['s']: raise Exception("ERROR must define a double-quoted search string with opt -s or as first argument. Run with -h for help")

  max_items=opt['M']  
  fetch_successful=False
  attempts=0
  while not fetch_successful and (  attempts < opt['a'] or opt['a']==0):
    try:
      search_handle=Entrez.esearch(db=db, term=opt['s'], retmode='xml', retmax=max_items)    
      record=Entrez.read(search_handle)      
      fetch_successful=True
    except KeyboardInterrupt:          raise
    except:
      service('Something failed (network problem?). Trying again... attempt number: '+str(attempts)+ ' (will keep trying'*int(opt['a']==0)+int(opt['a']!=0)*(' (will stop at attempt n.'+str(opt['a']))+')')
      sleep(1)
      attempts+=1

  if not fetch_successful:
    printerr("ERROR couldn't fetch results. Number of attempts: "+str(attempts)+'. (Last) exception was: ', 1)
    raise

  count= record["Count"]
  id_list= record["IdList"]
  #print record
  if len(id_list)==max_items and not opt['c']:   printerr('WARNING The maximum limit for the number of returned ids was hit ('+str(max_items)+')! you may want to rerun with a higher limit (option -M)', 1)

  if opt['c']:
    if not opt['silent']:    print(count)
  else:
    for id_t in id_list:
      if not opt['silent']: print(id_t)
      output_list.append(id_t)
  return output_list
  

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
