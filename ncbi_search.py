#! /usr/bin/python -u
from ncbi_lib import *
help_msg="""Utility to search/fetch any NCBI db and print info about results.

Usage #1 (search&fetch):  ncbi_search.py [dbname] -KEY1 "VALUE1" [-KEY2 "VALUE2"] [options]
Usage #2a (fetch):        ncbi_search.py [dbname] -i file_list_of_ids  [options]
Usage #2b (fetch):        ncbi_search.py [dbname] -I id1,id2,id3       [options]

Example:   ncbi_search.py  sra  -titl RRBS -orgn "Mus musculus"

Quotes can be omitted for single-word values.
For a list of databases, run:                     ncbi_db_info.py
For a list of keywords for a certain dbname, run: ncbi_db_info.py dbname

### options:
-o  f1,f2,f3  limits the output to these fields only
-x            do not attempt to expand XML strings for pretty printing
-t            tabular output; implies option -x
-r            reverse table; each entry becomes a line, each field a column
-print_opt    print currently active options
-h OR --help  print this help and exit"""

command_line_synonyms={}

def_opt= { #'temp':'/home/mmariotti/temp', 
'i':0, 'I':0, 'd':0, 'o':0,
'v':0, 'x':0, 't':0, 
'retmax':250, 'max_attempts':10, 'sleep_time':5, 
}

#########################################################
###### start main program function

from lxml import etree   #to parse xml and print it pretty

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'd', synonyms=command_line_synonyms, nowarning=1)
  else:  opt=args
  set_MMlib_var('opt', opt)
  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input

  dbname=opt['d']
  
  if not opt['i'] or opt['I']:
    ### building query for esearch
    query=''
    for k in opt:
      if not k in def_opt and not k.startswith('__'):
        query+='{v} [{k}] '.format(k=k, v=opt[k])
    if opt['v']:  printerr( 'search query: {}'.format(query), 1 )
    ids=esearch(db=dbname, term=query)
  else:
    ## ids provided
    ids=opt['I'].split(',') if not opt['i'] else  [line.strip() for line in open(opt['i']) if line.strip()]

  ### fetching entries from ncbi
  entries=   efetch(db=dbname, id=ids)
  if len(ids)!=len(entries): printerr('WARNING esearch returned {} ids, but efetch returned {} entries!'.format(len(ids), len(entries)), 1)
  if not entries: 
    printerr('No entries found!', 1)
    sys.exit(9)

  if opt['t']: opt['x']=1
  prefield='#'             if not opt['t'] else ''
  field_template='{:<12}'  if not opt['t'] else '{}'
  separator=':'            if not opt['t'] else '\t'

  possible_fields=None if not opt['o'] else    set( map(lower, opt['o'].split(','))  ) 
  if opt['r']:
    all_fields=set()
    for e in entries: all_fields.update(e.keys())
    all_fields=sorted(all_fields, cmp=lambda x,y:-1 if x=='Id' else   (+1 if y=='Id' else cmp(x,y))) 

  if opt['r']: write('\t'.join([f for f in all_fields if possible_fields is None or f.lower() in possible_fields]), 1)  #header

  for ei, e in enumerate(entries):
    if opt['r']:   target_fields=all_fields
    else:          target_fields=sorted(e.keys(), cmp=lambda x,y:-1 if x=='Id' else   (+1 if y=='Id' else cmp(x,y))) 
    for i, k in enumerate(target_fields):      
      if possible_fields is None or k.lower() in possible_fields:
        if opt['r'] and not k in e: 
          write('{}'.format('\t' if i else ''))   #filling a blank
          continue
        kstr='{}{}'.format(prefield, k)
        try:       vstr='{}'.format(e[k])
        except UnicodeEncodeError:
          vstr=e[k].encode('ascii', 'ignore')
        if not opt['x'] and vstr.startswith('<') and vstr.endswith('>'):
          try:
            x =  etree.fromstring('<ROOT>'+vstr+'</ROOT>')
            vstr='XML=\n' +  etree.tostring(x, pretty_print=True)  [7:-9]  #removing root
          except XMLSyntaxError: pass
        
        if opt['r']: 
          write('{}{}'.format('\t' if i else '',   vstr))
        else:
          write(field_template.format(kstr), how=  ['blue','green'][ei%2]  )
          write('{} {}'.format(separator, vstr), 1)
    write('',1)

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
