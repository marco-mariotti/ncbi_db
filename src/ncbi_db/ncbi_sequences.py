#!/usr/bin/env python3
from string import *
import sys
from subprocess import *
from .MMlib3 import *
from .ncbi_taxonomy import main as run_ncbi_taxonomy
from .ncbi_lib import email_setup
from Bio import Entrez, SeqIO
import time

help_msg="""Program to fetch entries from the ncbi sequence databases through the internet and print information such as its sequence or the source organism. It wraps Bio->Entrez and Bio->SeqIO.

ncbi_sequences  [mode] [ids] [options]

## Define mode (which ncbi database is queried):
-m  +      Possibilities: N -> nucleotide, P -> protein

## Define IDs; anything that ncbi understand is accepted (e.g. AAB88790.1):
-i  +      define inputfile  with IDs, one per line
-I  +[,+]  define IDs on the command line, comma separated

## Define output: (if multiple output is defined, they are all printed)
-M         minimal header: only the primary accession id is printed
-L         one entry per line: no newlines are printed, tab separators are printed between output information
-f         fasta sequence output
-o         species output 
-t         taxonomy+species output (this is slower than just -o, since an additional web query is performed)
-x         output all available fields, but do not output sequence

The output contains the desired fields (fasta, species, taxonomy) in this order. So if you asked for fasta and taxonomy for example ( -f -t ) you will have: 

>AAB88790 selenophosphate synthetase [Drosophila melanogaster]
7227 : cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Coelomata; Protostomia; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Hexapoda; Insecta; Dicondylia; Pterygota; Neoptera; Endopterygota; Diptera; Brachycera; Muscomorpha; Eremoneura; Cyclorrhapha; Schizophora; Acalyptratae; Ephydroidea; Drosophilidae; Drosophilinae; Drosophilini; Drosophilina; Drosophiliti; Drosophila; Sophophora; melanogaster group; melanogaster subgroup -> Drosophila melanogaster
MSYAADVLNSAHLELHGGGDAELRRPFDPTAHDLDASFRLTRFADLKGRGCKVPQDVLSK
LVSALQQDYSAQDQEPQFLNVAIPRIGIGLDCSVIPLRHGGLCLVQTTDFFYPIV

### Options:
-email +         email address to use with NCBI Entrez. Required for first usage of ncbi online tools
-a     +         max attempts numbers. Sometimes unstability of network results in failed attempts to reach ncbi. The program tries by default this number of times before giving up. There's a 1 sec wait between attempts. Set -a to 0 to never give up. Careful: if this number is high (or 0), it means it takes a long time before it realizes there's a problem, like: no internet connection.
-w     +         compulsory waiting time between NCBI connections
-b               batch size. default is 50, but it is lowered by 15% if the following error is found: IOError("Requested URL too long (try using EPost?)") 
-print_opt       print currently active options
-h OR --help     print this help and exit

### When used as function imported in another python program, the "main" function must be run using the opt hash which you would use to get the information, and using the opt['silent']=1.  Example:
from ncbi_db import ncbi_sequences; out_hash=ncbi_sequences.main({'I':'AAB88790', 'f':1, 'silent':1})
An out_hash is returned,   with as key the id provided and as value a list [full_id, seq_if_requested, taxonomy_+_species_if_requested]. Full_id is minimal if option M is active
"""

command_line_synonyms={}

def_opt= {'temp':'/tmp/', 
          'm':'', 'a':100,
          'i':0, 'I':'',
          'M':0, 'L':0,
          'f':0, 't':0, 'o':0,
          'v':0, 'Q':0, 'x':0,
          'w':0.5,
          'b':50, 'silent':0,
          'email':'',
}

mode_synonyms={'N':'nucleotide', 'P':'Protein'}

#########################################################
###### start main program function

def main(args={}):
#########################################################
  def id_is_included_in_list(main_id, id_list):
    for a_id in id_list:
      if main_id in a_id:
        return True
    return False


############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms )
  else:  
    opt=args
    for k in def_opt: 
      if k not in opt: opt[k]=def_opt[k]
      
  def verbose(msg):
    if opt['v']:       printerr( msg, 1)
  #global temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 


  #mode
  if not opt['m']:  raise Exception("ERROR you must specify a mode with option -m (N|P). See --help")
  try:
    db=mode_synonyms[opt['m'].upper()] 
  except: raise Exception("ERROR mode not recognized: "+opt['m'])
  opt['m']=opt['m'].upper()
  #determining input ids
  id_list=[]
  if opt['i']:
    check_file_presence(opt['i'], 'input_file_with_ids')
    for line in open(opt['i']):       id_list.append(line[:-1])
  elif opt['I']:    id_list=str(opt['I']).split(',')
  else: raise Exception("ERROR must provide ids either with option -i FILE or option -I id1,id2,...,idN. Run with option -h for help")
  #output preparing
  line_sep='\n'
  if opt['L']: line_sep='\t'
  #batch preparation
  batch_size=opt['b']
  batch_index=0
  id_retrieved_until_now=0
  headers_checklist={} ### to deal with duplicate fetching

  output_hash={}

  email_setup(opt['email'])
  
  while id_retrieved_until_now!=len(id_list):
    try:
      batch_id_list=id_list[id_retrieved_until_now:id_retrieved_until_now+batch_size]
      if len(batch_id_list)!=len(id_list) and not opt['silent']:
        verbose('Batch n'+str(batch_index+1)+' --- entries:'+str(len(batch_id_list)))
      
      #retrieving! Seq_accver': 'CM000066.5', u'TSeq_sequence
      rettype='gb'
      if opt['m']=='N' and not opt['x'] : rettype='fasta'

      number_of_records=0
      batch_id_list_index=0

      fetch_successful=False
      error_occurred=None
      attempts=0
      while not fetch_successful and (  attempts < opt['a'] or opt['a']==0):
        if opt['w'] and attempts>0:  time.sleep(opt['w'])
        try:
          #print( f"db: {db} id: {','.join(batch_id_list)} retmode=xml  rettype={rettype}")
          fetch_handle=Entrez.efetch(db=db, id=','.join(batch_id_list), retmode='xml', rettype=rettype)
          parsed_results=[r for r in Entrez.parse(fetch_handle)]
          fetch_successful=True
        except KeyboardInterrupt:          raise
        except Exception as error: ##debug
          #print('Something failed (network problem?). Trying again... attempt number: '+str(attempts)+ ' (will keep trying'*int(opt['a']==0)+int(opt['a']!=0)*(' (will stop at attempt n.'+str(opt['a']))+')')            
          if not opt['silent']:
            service('Something failed (network problem?). Trying again... attempt number: '+str(attempts)+ ' (will keep trying'*int(opt['a']==0)+int(opt['a']!=0)*(' (will stop at attempt n.'+str(opt['a']))+')')
            time.sleep(1)
          error_occurred=error
          attempts+=1

      if not fetch_successful:
        printerr("ERROR couldn't fetch results. Number of attempts: "+str(attempts)+'. (Last) exception was: ', 1)
        raise error_occurred
      
      for record in parsed_results:
        if opt['m']=='N':
          #[u'TSeq_accver', u'TSeq_sequence', 'GBSeq_other-seqids', u'TSeq_length', u'TSeq_taxid', u'TSeq_orgname', u'TSeq_sid', u'TSeq_gi', u'TSeq_seqtype', u'TSeq_defline']
          #interactive_mode
          if opt['x']:
            write('#'*100, 1, how='blue')
            for k in sorted(record.keys()):
              write( str(k.ljust(25)), 0, how='green')
              write(' '+str(record[k]), 1)
              header=''
          else: 
            seq=record['TSeq_sequence'].upper()
            organism=record['TSeq_orgname']
            if 'TSeq_accver' in record:
              header=record['TSeq_accver']
            else: 
              if 'TSeq_sid' not in record: raise Exception("ERROR record format not recognized for this entry: "+str( '\n'. join([str(k)+': '+str(record[k]) for k in record])))
              header=record['TSeq_sid']

            if not opt['M']:
              ## COMPUTING NICE STRING FOR IDS
              #header=header+" gi|"+record['TSeq_gi']+'|'+header+' '+record['TSeq_defline']
              header=header+' '+record['TSeq_defline']
            if header.split()[0] in headers_checklist: continue
            headers_checklist[header.split()[0]]=True

          #print header
          
          record['GBSeq_other-seqids']=[[batch_id_list[batch_id_list_index]]]   ## to keep homogenous procedures below...

        elif opt['m']=='P':        

          if opt['x']:
            write('#'*100, 1, how='blue')
            for k in sorted(record.keys()):
              write( str(k.ljust(20)) + ' '+str(record[k]), 1)
              header=''
          if True: #
            seq=record['GBSeq_sequence'].upper()          
            organism=record['GBSeq_organism']
            #taxonomy=record['GBSeq_taxonomy']+' -> '+organism   ### ---> this is approximate and does not contain all the sublevels. also, it is not available for nucleotides, which makes the whole thing undesirable.

            #for k in record:
            #  print k+' '+str(record[k])[:100]
            #raw_input('..')
            header=  record['GBSeq_accession-version'] if 'GBSeq_accession-version' in record else    record['GBSeq_primary-accession']
            if not opt['M']:       
            ## COMPUTING NICE STRING FOR IDS
              gi_id=''; ref_id=''; other_ids=[]
              for a_id in record['GBSeq_other-seqids']:
                if    a_id.split('|')[0]=='gi': gi_id=a_id
                elif  a_id.split('|')[0]=='ref': 
                  ref_id=a_id
                else: other_ids.append(a_id)
              id_string=''
              if gi_id:   id_string+=gi_id
              if ref_id:  
                if id_string and id_string[-1]!='|': id_string+='|' #adding a | if is not already there
              id_string+=  ref_id
              for a_id in other_ids:
               if id_string and id_string[-1]!='|': id_string+='|' 
               id_string+=  a_id 
              header+=' '+id_string+' '+record['GBSeq_definition']        

          if header.split():
            if header.split()[0] in headers_checklist: continue        
            headers_checklist[header.split()[0]]=True
              
          ##########
            

         
        while batch_id_list_index < len(batch_id_list) and not id_is_included_in_list( batch_id_list[batch_id_list_index], record['GBSeq_other-seqids'] ):
          if not opt['silent']:  printerr("WARNING! entry " +batch_id_list[batch_id_list_index]+"  was not found and was ignored", 1)
          output_hash[batch_id_list[batch_id_list_index]]=[]
          batch_id_list_index+=1

        if batch_id_list_index == len(batch_id_list): sys.exit()
        output_hash[batch_id_list[batch_id_list_index]]= [header]
          
        out= '>'+header+line_sep

        if opt['o']:
          output_hash[batch_id_list[batch_id_list_index]].append(organism)        
          out+=organism+line_sep
          
        if opt['t']:
          ncbi_tax_opt={'silent':True, 'S':organism, 'a':opt['a']}
          results=run_ncbi_taxonomy( ncbi_tax_opt )
          if len(results)==0:
            if not opt['silent']: printerr("WARNING! organism " +organism+" returned no results in taxonomy! keeping this empty", 1)
            taxonomy = 'unknown'
          else:
            if len(results)>1: 
              if not opt['silent']: printerr("WARNING! organism " +organism+" returned more than a result in taxonomy! taking a random one", 1)
#            raise Exception, "ERROR searching species in ncbi_taxonomy: "+organism+" ; number of results: "+str(len(results))
            taxonomy = results[ list(results.keys())[0]  ]
  
          output_hash[batch_id_list[batch_id_list_index]].append(taxonomy)
          out+=taxonomy+line_sep
        
        if opt['f']: #fasta
          output_hash[batch_id_list[batch_id_list_index]].append(seq)
          if opt['L']:         out+=seq+line_sep
          else:                out+= fasta(seq)+line_sep

        if not opt['silent'] and not opt['x']: write(out[:-1], 1)

        batch_id_list_index+=1
        number_of_records+=1
      id_retrieved_until_now+=len(batch_id_list)

      if batch_id_list_index<len(batch_id_list): 
        for i in range(batch_id_list_index,len(batch_id_list)):  
          output_hash[batch_id_list[batch_id_list_index]]=[]
          if not opt['silent']: printerr("WARNING! entry " +batch_id_list[i]+"  was not found and was ignored", 1)

      batch_index+=1
      

    except IOError as e:
      if e.args and 'Requested URL too long' in e.args[0]:
        batch_size=int(batch_size*0.85)
        if not opt['silent']: verbose('requested URL too long... trying to decrease batch size to:'+str(batch_size))
      else: raise

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
