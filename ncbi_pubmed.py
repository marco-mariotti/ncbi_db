#! /usr/bin/env python3
from ncbi_lib import *
help_msg="""Utility to search pubmed entries online and get formatted reference-style entries 
Usage:
$ ncbi_pubmed.py   [search]   [output]

###  Search options:
-P   pubmed id(s), comma-separated  |  -p   pubmed ids file
-kp  citationKey -tab- pubmed id file; output will be in style citationKey -tab- citation

##  instead of providing pubmed ids you can combine any of these:
-S   default search; this is like standard pubmed search, you can also include field names e.g. "Nature [JOUR]"
-T   title search
-A   author; use double quotes to include name or initial; comma-separate multiple authors
-D   publication date (year only is accepted)

###  Output options:
#    Reference style:
-s  [style]   reference style is in the specified format. Available: pnas, elife, info
-ma max authors listed (if different than default for style)
#    Format options: 
-f  [format]  Use html and open with a browser to have formatted text. Available: html, txt
#    Misc:
-sa  sort per citation_in_text
-ct  text only output: id|citationInText    
     where citationInText is like "Mariotti et al., 2015" or "Chapple and Guigo, 2008" and id is the pubmed id or the citationKey if -kp

### other options:
-print_opt      print currently active options
-h OR --help    print this help and exit"""

command_line_synonyms={}

def_opt= { #'temp':'/users-d3/mmariotti/temp', 
'P':0, 'p':0, 'kp':0,
's':'info',
'f':'txt',
'ma':0,
'S':0, 'T':0, 'A':0, 'D':0,
'retmax':250, 'max_attempts':10, 'sleep_time':5, 
'o':0,
'ct':0, 'sa':0,
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
  if not output_format in available_output_formats: raise Exception('ERROR invalid argument to option -f !  See -h')

  available_ref_style={'pnas':1, 'info':1, 'elife':1}
  ref_style=opt['s']
  if not ref_style in available_ref_style: raise Exception('ERROR invalid argument to option -s !  See -h')

  style2max_authors={'pnas':5, 'elife':30, 'info':0} #not used for info
  max_authors= opt['ma'] if opt['ma'] else style2max_authors[ref_style]
  
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
    pubmed_ids=list(pubmed_id2citation_key.keys())
    #print pubmed_id2citation_key
  else: 
    term=''
    if opt['T']: 
      for i in opt['T'].split(): term+= ''+i +'[TITL] '
    if opt['A']: 
      for i in opt['A'].split(','): term+= ''+i.strip() +'[AUTH] '
    if opt['D']: 
      term+= ' '+str(opt['D']) +'[PDAT] '
    if opt['S']: term+= ''+ opt['S'] +' ' #[ALL] '
    #write(term, 1, how='red')
    if not term:
      raise Exception("ERROR no article search defined. See -h")
    
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

  def get_misc(entry): 
      all_authors=entry['AuthorList']
      year=entry['PubDate'].split()[0]  #so_bit.split()[0];     
      author_summary=''
      citation_in_text=''
      for auth_index, author in enumerate(all_authors):
        splt=author.split()
        surname=''; name_bit=''
        for splt_index, splt_bit in enumerate(splt):
          if splt_bit.upper() == splt_bit: name_bit+=splt_bit+', '
          else:                           surname+=splt_bit+' '
        
        surname=surname.strip()
        name_bit=name_bit.strip();         name_bit=name_bit.strip(','); 
        author_summary+=   surname + ' '+name_bit  + ', '

        if    auth_index==0: citation_in_text+=surname
        elif  auth_index==1:
          if len(all_authors)==2: citation_in_text+=' and '+surname
          else:                   
            #if citation_in_text[-2:]==', ':             citation_in_text=citation_in_text[:-2]
            citation_in_text+= ' et al.'

        if len(all_authors)>max_authors: 
          if author_summary[-2:]==', ':             author_summary=author_summary[:-2]
          author_summary+=' et al. '; break
      author_summary=author_summary.strip().strip(',')
      citation_in_text+=', '+str(year)
      return author_summary, citation_in_text

  if opt['sa']:    entries.sort(key=lambda x:get_misc(x)[1])

  for entry in entries:
    #print entry
    pubmed_id=entry['Id'];   doi=entry['DOI'] if 'DOI' in entry else 'None'
    journal_brief=entry['Source']; journal_full=entry['FullJournalName']
    pages=entry['Pages']
    ref_count=entry['PmcRefCount'];   
    all_authors=entry['AuthorList']
    so_bit=entry['SO'];     
    year=entry['PubDate'].split()[0]  #so_bit.split()[0];     
    where_bit=';'.join(so_bit.split(';')[1:])
    title=entry['Title']
    pubstatus=entry['PubStatus'] 
    ## common ?
    if pubstatus=='aheadofprint' and not where_bit:   where_bit='Published online {date} doi:{doi}'.format(date=entry['EPubDate'], doi=doi)
    author_summary, citation_in_text =  get_misc(entry) 

    #### opt -kp
    citk=pubmed_id2citation_key[pubmed_id]+' ' if opt['kp'] else ''
      
    if opt['ct']: 
      if opt['kp']: textout+=pubmed_id2citation_key[pubmed_id]+'|'+citation_in_text+'\n'
      else:         textout+=pubmed_id+'|'+citation_in_text+'\n'
      continue ### no other output!


    if ref_style == 'info':
      textout+='\n### PUBMED ID: '+pubmed_id+'\n'
      for k in sorted(entry.keys()):
        try: val=str(entry[k])
        except: val=entry[k].encode('ascii', errors='replace')
        textout+=  k.ljust(20)+' : '+val  +'\n'

    else:

      if ref_style =='pnas':

        if output_format=='txt':
          textout+='{citk}{aut} ({year}) {title} {jour} {where}\n'.format(aut=author_summary, year=year, title=title, jour=journal_brief, where=where_bit, citk=citk)
        elif output_format=='html':
          textout+='<p>{citk}{aut} ({year}) {title} <i>{jour}</i> {where}</p>'.format(aut=author_summary, year=year, title=title, jour=journal_brief, where=where_bit, citk=citk)

      elif ref_style =='elife':
        splt=author_summary.split(', ')
        first_author=splt[0];     rest_authors=', '.join(splt[1:])
        if rest_authors: rest_authors=', '+rest_authors
        splt=where_bit.split(':')
        first_where=splt[0];     rest_where=':'.join(splt[1:])
        if '(' in first_where: first_where=first_where.split('(')[0]
        if rest_where: rest_where=':'+rest_where
        doi_bit=' doi: '+doi+'.' if doi !='None' else ''

        if output_format=='txt': pass
          #textout+='{citk}{aut} ({year}) {title} {jour} {where}\n'.format(aut=author_summary, year=year, title=title, jour=journal_brief, where=where_bit, citk=citk)

        elif output_format=='html':
          ### formatted
          #textout+=u'<p>{citk}<b>{faut}</b>{raut}. {year}. {title} <i>{jour}</i> <b>{fwhere}</b>{rwhere}. <font color="blue">{doi}</font>.</p>'.format(faut=first_author, raut=rest_authors, year=year, title=title, jour=journal_brief, fwhere=first_where, rwhere=rest_where, doi=doi_bit,  citk=citk)
          textout+='<p>{citk}{faut}{raut}. {year}. {title} <i>{jour}</i> <b>{fwhere}</b>{rwhere}.{doi}</p>'.format(faut=first_author, raut=rest_authors, year=year, title=title, jour=journal_brief, fwhere=first_where, rwhere=rest_where, doi=doi_bit,  citk=citk)
        

    

  if output_format=='html':  textout+='</body>'

  if output_format=='txt':
    print(textout.encode('utf8', errors='replace'))
  elif output_format=='html':
    print(textout.encode('ascii', 'xmlcharrefreplace'))

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
