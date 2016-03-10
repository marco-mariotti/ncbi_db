#! /usr/bin/python -u
from ncbi_lib import *
help_msg="""Utility to search pubmed entries online and get formatted reference-style entries 
Usage:
$ ncbi_pubmed.py   [search]   [output]

###  Search options:
-P   pubmed id(s), comma-separated  |  -p   pubmed ids file
-kp  citationKey -tab- pubmed id file; output will be in style citationKey -tab- citation
##  instead of providing pubmed ids you can combine any of these:
-T   title search
-A   author search
-S   default search; this is like standard pubmed search, you can also include field names e.g. "Nature [JOUR]"

###  Output options:
#    Reference style:
-s  [style]   reference style is in the specified format. Available: pnas, info
#    Format options: 
-f  [format]  Use html and open with a browser to have formatted text. Available: html, txt

### other options:
-print_opt      print currently active options
-h OR --help    print this help and exit"""

command_line_synonyms={}

def_opt= { #'temp':'/home/mmariotti/temp', 
'P':0, 'p':0, 'kp':0,
's':'info',
'f':'txt',
'S':0, 'T':0, 'A':0, 
'retmax':250, 'max_attempts':10, 'sleep_time':5, 
'o':0,
'v':0,
}

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

  available_output_formats={'html':1, 'txt':1}
  output_format=opt['f']
  if not output_format in available_output_formats: raise Exception, 'ERROR invalid argument to option -f !  See -h'

  available_ref_style={'pnas':1, 'info':1}
  ref_style=opt['s']
  if not ref_style in available_ref_style: raise Exception, 'ERROR invalid argument to option -s !  See -h'
  
  ###   determining pubmed ids
  if opt['p'] or opt['P']:
    pubmed_ids= str(opt['P']).split(',') if opt['P'] else  [line.strip() for line in open(opt['p']) if line.strip() ]
  elif opt['kp']:
    pubmed_id2citation_key={}
    for line in open(opt['kp']):
      if not line.strip(): continue
      splt=line.strip().split('\t')
      citation_key, pubmed_id= splt
      pubmed_id2citation_key[pubmed_id]=citation_key
    pubmed_ids=pubmed_id2citation_key.keys()
    #print pubmed_id2citation_key
  else: 
    term=''
    if opt['T']: 
      for i in opt['T'].split(): term+= ''+i +'[TITL] '
    if opt['A']: 
      for i in opt['A'].split(): term+= ''+i +'[AUTH] '
    if opt['S']: term+= ''+ opt['S'] +''
    if not term:
      raise Exception, "ERROR no article search defined. See -h"
    
    if opt['v']:  printerr( 'search term: ' +term, 1)
    pubmed_ids=esearch(db='pubmed', term=term)

  ### fetching entries from ncbi
  entries=   efetch(db='pubmed', id=pubmed_ids)
  if not entries: 
    printerr('No entries found!', 1)
    sys.exit(9)

  ## building output
  textout=''
  if output_format=='html':  textout+='<html><head><title>ncbi_pubmed output</title></head><body>'
  for entry in entries:
    #print entry
    pubmed_id=entry['Id'];   doi=entry['DOI'] if 'DOI' in entry else 'None'
    journal_brief=entry['Source']; journal_full=entry['FullJournalName']
    pages=entry['Pages']
    ref_count=entry['PmcRefCount'];   
    all_authors=entry['AuthorList']
    so_bit=entry['SO'];     
    year=entry['PubDate'].split()[0]  #so_bit.split()[0];     
    where_bit=join(so_bit.split(';')[1:], ';')
    title=entry['Title']
    pubstatus=entry['PubStatus'] 
    ## common ?
    if pubstatus=='aheadofprint' and not where_bit:   where_bit='Published online {date} doi:{doi}'.format(date=entry['EPubDate'], doi=doi)

    if ref_style == 'info':
      textout+='\n### PUBMED ID: '+pubmed_id+'\n'
      for k in sorted(entry.keys()):
        try: val=str(entry[k])
        except: val=entry[k].encode('ascii', errors='replace')
        textout+=  k.ljust(20)+' : '+val  +'\n'

    elif ref_style =='pnas':
      
      ##author
      author_summary=''
      for author in all_authors:
        splt=author.split()
        author_summary+=' '+splt[0]+(' '+splt[1] if len(splt)>1 else '')+(', '+join(splt[2:], ', ') if len(splt)>2 else '')+','
        if len(all_authors)>5: author_summary+=' et al. '; break
      author_summary=author_summary.strip().strip(',')
        
      #### out
      citk=pubmed_id2citation_key[pubmed_id]+' ' if opt['kp'] else ''
        
      if output_format=='txt':
        textout+='{citk}{aut} ({year}) {title} {jour} {where}\n'.format(aut=author_summary, year=year, title=title, jour=journal_brief, where=where_bit, citk=citk)
      elif output_format=='html':
        textout+=u'<p>{citk}{aut} ({year}) {title} <i>{jour}</i> {where}</p>'.format(aut=author_summary, year=year, title=title, jour=journal_brief, where=where_bit, citk=citk)

   
    #all_authors_str=join
    

  if output_format=='html':  textout+='</body>'

  if output_format=='txt':
    print textout.encode('utf8', errors='replace')
  elif output_format=='html':
    print textout.encode('ascii', 'xmlcharrefreplace')

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
