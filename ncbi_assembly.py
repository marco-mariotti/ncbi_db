#! /usr/bin/python -u
from string import *
import sys
from commands import *
sys.path.append('/home/mmariotti/scripts')
from MMlib import *
from Bio import Entrez 
from urllib2 import URLError
from time import strptime, sleep

Entrez.email = "marco.mariotti@crg.es"  #for Entrez
help_msg="""Script that interrogates ncbi online through the Bio.Entrez module and retrieves information about the assemblies for a certain organism or lineage. Normally the 'genome' db is queried, and the corresponding 'assembly' entries are fetched. If you wish, you can query the 'assembly' db instead with -A, which ends up in more results (potentially more than 1 per species).
## Usage:   $ ncbi_assembly.py   -S species_or_lineage   [options]

### Options:
-a              shows detailed information for each assembly entry. As arguments:     
     0           -> shows nothing (quiet mode)
     1 [default] -> shows a few built-in selected fields
     2           -> shows all the non-null fields in this object
     f1,f2,..    -> shows this list of comma-separated fields 
     +f1,f2,..   -> add these fields to the default fields
 Note: for each assembly the attributes available are those found in the xml object returned by the search; additionally, these are also derived and made available as fields:
   props      all strings in PropertyList joined together
   date       NCBIReleaseDate made shorter
   ftp        ftp folder for this assembly, derived from the string in the Meta attribute
   ftp:XXX    derived path of a specific file. Other options can accept these "XXX" file classes. Here's the possible ones:
 | dna   genome file   (*_genomic.fna)  | gff   annotation file    (*_genomic.gff)
 | pep   proteome file (*_protein.faa)  | fea   feature table file (*_feature_table.txt)
 | gbf   genbank format annoation file  | md5  md5checksum file                 

-G        query the 'genome' db, instead of 'assembly' directly; this gets fewer entries but better quality
 ## if -G is active:
 -g       show detailed  information for each genome entry. Syntax just like -a (no additional fields available though)
 -z       display also the entries with AssemblyID=0; normally they are skipped
 ## if -G is NOT active:
 -e       do not merge different versions of the same assembly (normally only the latest version is displayed)
 -n       do not keep only the newest assembly per species
 -tax     provide a ncbi taxid; in this way you don't have to provide option -S. Note: only exact matches are found, the descendants are not reported (unlike with -S)

## check presence of web files and/or download them
-c  [x1,x2..]   test if certain files are found in the ncbi ftp site. Possible x values: see XXX file classes above. With no argument, they are all tested
 ## if option -c  is active:
 -s             show file size (compressed) instead of just presence
-d  [x1,x2..]   download certain files; syntax just as -c. Md5sum is checked (HAVE TO IMPLEMENT MD5SUM)
 ## if option -d  is active:
 -f             download master folder. The files are downloaded in a subfolder named after the species. By default they are also uncompressed and linked
 -dl            do not link files with a standard name (e.g. genome.fa for dna, proteome.fa for pep); if used as "-dl 2", files are not even uncompressed with gunzip
 -F             force download, even if files are found in the local master folder 

## technical
-retmax         max items requested at each connection to ncbi
-max_attempts   max attempts querying a database
-sleep_time     time between attempts in seconds

## other options
-tab          tab separated output. The main messages then go to stderr, to allow redirection of the tables to an output file
-v            verbose output, for debugging
-print_opt    print currently active options
-h OR --help  print this help and exit"""

command_line_synonyms={}

def_opt= { 'temp':'/home/mmariotti/temp', 
'S':'',   'tax':False,
'z':0, 'e':0, 'n':0,
'tab':0, 
'c':0, 's':0, 'd':0, 
'f':'', 'F':0, 'dl':0, 
'G':0,   
'g':1, 'a':1,
'retmax':250, 'max_attempts':10, 'sleep_time':5, 
'v':0, 
}

ftp_file_description={'dna':'genome.fa',  'pep':'proteome.fa', 
'gff':'annotation.gff', 'fea':'feature_table.txt', 'gbf':'genbank.gbff',
'md5':'md5sum.txt'}
ftp_file_types=ftp_file_description.keys()


def populate_with_file_paths(assembly_e):
  """ Given an assembly entry, it parses its fields and creates some fields with the path to files that can be downloaded"""
  assembly_e['ftp']= assembly_e['Meta'].split('<FtpSites> ')[1].split('</FtpSites>')[0].split('<FtpPath type="GenBank">')[1].split('</FtpPath>')[0].strip()
  last_bit=assembly_e['ftp'].split('/')[-1]
  assembly_e['ftp:dna']= '{0}/{1}_genomic.fna.gz'.format(assembly_e['ftp'], last_bit)
  assembly_e['ftp:gff']= '{0}/{1}_genomic.gff.gz'.format(assembly_e['ftp'], last_bit)
  assembly_e['ftp:gbf']= '{0}/{1}_genomic.gbff.gz'.format(assembly_e['ftp'], last_bit)
  assembly_e['ftp:pep']= '{0}/{1}_protein.faa.gz'.format(assembly_e['ftp'], last_bit)
  assembly_e['ftp:fea']= '{0}/{1}_feature_table.txt.gz'.format(assembly_e['ftp'], last_bit)
  assembly_e['ftp:md5']= '{0}/md5checksums.txt'.format(assembly_e['ftp'])

default_fields_displayed_genome  =['Id', 'Organism_Name', 'DefLine', 'Assembly_Accession', 'AssemblyID']
default_fields_displayed_assembly=['ChainId', 'AssemblyAccession','RsUid', 'SpeciesName', 'AssemblyStatus', 'date', 'AssemblyName'] #, 'props'] 
def message(msg):
  """Function to handle printing the main messages of the program, e.g. start, how many hits, how many removed and so on """
  if opt['tab']: printerr('#'+str(msg).center(90, '=')+'#', 1, how='green')
  else:          write(   '#'+str(msg).center(90, '=')+'#', 1, how='green')

#########################################################
###### start main program function

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'S', synonyms=command_line_synonyms )
  else:  opt=args
  set_MMlib_var('opt', opt)
  global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input

  ### checking for illegal options
  species=opt['S']
  if not species and not opt['tax']:  raise Exception, "ERROR you must specify a species or lineage name with option -S (or alternatively a taxid with -tax). See --help"
  if species and opt['tax']:          raise Exception, "ERROR option -S and option -tax are mutually exclusive. See --help"
  if opt['tax'] and opt['G']:         raise Exception, "ERROR option -tax cannot be used with -G. See --help"
  if opt['n'] and opt['G']:           raise Exception, "ERROR option -n makes no sense with -G. See --help"
  if opt['e'] and opt['G']:           raise Exception, "ERROR option -e makes no sense with -G. See --help"
  if opt['s'] and not opt['c']:       raise Exception, "ERROR option -s makes sense only if -c is active. See --help"
  if opt['d']: 
    if opt['d']==1 or any([not c in ftp_file_types for c in opt['d'].split(',')]): raise Exception, "ERROR illegal argument for option -d ! See -help"
    if not opt['f']:                  raise Exception, "ERROR to download (option -d) you must provide a local folder with option -f. See --help"
    local_master_folder=Folder(opt['f']); test_writeable_folder(local_master_folder)
  if opt['c'] and opt['c']!=1 and any(  [not c in ftp_file_types for c in opt['c'].split(',') ] ): raise Exception, "ERROR illegal argument for option -c ! See -help"
  #########

  ####### program start
  ## NOTE message is used to display any non-data message, e.g. how many results, how many filtered etc
  message(' Searching for "{0}" in ncbi database "{1}" '.format(species, {False:'assembly',True:'genome'}[bool(opt['G'])]))

  if opt['G']:
    assembly_uid_list=[]  ## we have to populate this, either by querying genome db first, or by querying directly assembly db (-A)
    discarded_entries=0
    ## searching genomes for this lineage (or species)
    genome_uid_list= esearch(db='genome', term=species, field='Organism')  
    if not genome_uid_list:  message(' No genome entries were found -- Exiting... ' .center(100, '=')); return 
    genome_entries = efetch (db='genome', id=genome_uid_list)
    ## if genome_entries had a different number of entries than genome_uid_list, efetch would have crashed
    message(' {0} genome entries found '.format( len(genome_entries) ))
    genome_entries.sort(key=lambda x:x['Organism_Name'])
    ### if we're printing info, let's prepare the stage
    if opt['g']:  #and genome_entries    implicit
      fields_displayed=default_fields_displayed_genome
      if opt['g']==1:    pass
      elif opt['g']==2:  fields_displayed=genome_entries[0].keys()  ## assuming all entries have the same fields
      else:              
        if opt['g'].startswith('+'): fields_displayed.extend( opt['g'][1:].split(',') )
        else:                        fields_displayed=opt['g'].split(',')
      max_length_per_fields={}
      for field in fields_displayed:  max_length_per_fields[field]= max([len(str(e[field])) for e in genome_entries] + [len(field)]) #getting max string length for pretty tabular like output
      if opt['tab']:  fields_displayed.sort()   ##sorting alphabetical
      else:           fields_displayed.sort(key=lambda x:max_length_per_fields[x])
      if opt['tab']: header_line = join(fields_displayed, '\t')
      else:          header_line = join([field.ljust(max_length_per_fields[field]) for field in fields_displayed ], ' ')
      write(header_line, 1, how='reverse')

    #### cycle for printing information, and also storing assembly ids
    for genome_e in genome_entries:
      if not genome_e['AssemblyID']=='0': assembly_uid_list.append( genome_e['AssemblyID'] )
      else: 
        discarded_entries+=1
        if not opt['z'] and genome_e['AssemblyID']=='0': continue  ## skipping strange entries
      if opt['g']:
        if opt['tab']: summary_line=join([ str(genome_e[field]) for field in fields_displayed ], '\t') 
        else:          summary_line=join([ str(genome_e[field]).ljust(max_length_per_fields[field]) for field in fields_displayed ], ' ')
        write(summary_line, 1)
    if discarded_entries:   message(' {0} entries were discarded, since missing an AssemblyID '.format( discarded_entries ) ) 
 
  elif not opt['tax']:  ### querying assembly db directly   
    ######## DEFAULT BEHAVIOUR
    assembly_uid_list= esearch(db='assembly', term=species, field='Organism')  
  else:
    taxid=str(opt['tax'])
    assembly_uid_list= esearch(db='assembly', term=taxid, field='Taxonomy ID')  

  if not assembly_uid_list:  message(' No assembly entries were found -- Exiting... '); return       
  #############
  ######  now getting assembly entries
  assembly_entries = efetch (db='assembly', id=assembly_uid_list)
  ## if assembly_entries had a different number of entries than assembly_uid_list, efetch would have crashed

  ## deriving ftp folder and other useful attributes for each assembly
  for assembly_e in assembly_entries:
    assembly_e['ftp']='Failed to parse!'
    for ft in ftp_file_types: assembly_e['ftp:'+ft]='None'
    try:        populate_with_file_paths(assembly_e)
    except: pass
    assembly_e['props']= join(assembly_e['PropertyList'], ' ')
    assembly_e['date']= assembly_e['NCBIReleaseDate'].split()[0]

  ####################
  ####### filtering the assembly entries
  if not opt['G'] and not opt['e']:
    # removing duplicates of different versions, e.g. GCA_000188675.1  GCA_000188675.2 
    assembly_accession_dict={} # k: root_accession -> value:  [best_version_accession, best_version_index, [indexes_to_remove...]]
    for index, assembly_e in enumerate(assembly_entries):
      try:      root_accession, version_accession =assembly_e['AssemblyAccession'].split('.');   version_accession=int(version_accession)
      except:   write( "ERROR trying to dot-split the AssemblyAccession id in this entry: "+str(assembly_e), 1); raise
      if not root_accession in assembly_accession_dict:    assembly_accession_dict[root_accession]=[version_accession, index, []] 
      else: 
        if assembly_accession_dict[root_accession][0] > version_accession:    assembly_accession_dict[root_accession][2].append(index)  #stored is already better
        else:   assembly_accession_dict[root_accession]= [ version_accession, index, assembly_accession_dict[root_accession][2]+[assembly_accession_dict[root_accession][1]] ]  
    indexes_to_remove=[]
    for root_accession in assembly_accession_dict: indexes_to_remove.extend( assembly_accession_dict[root_accession][2] )
    for index in sorted(indexes_to_remove, reverse=True): assembly_entries.pop(index)
    if indexes_to_remove: message(' Removed {0} assembly entries with a newer version available '.format( len(indexes_to_remove) ) )
    del assembly_accession_dict; del indexes_to_remove

  if not opt['G'] and not opt['n']:
    # removing duplicates of different species, keeping the most recent one
    species_dict={} # k: species_name -> value:  [most_recent_date, most_recent_index, [indexes_to_remove...]]
    for index, assembly_e in enumerate(assembly_entries):
      try:      
        this_species=assembly_e['SpeciesName']
        this_date   =strptime(assembly_e['date'], "%Y/%m/%d")        
      except:   write("ERROR trying to get SpeciesName and Date in this entry: "+str(assembly_e), 1); raise
      if not this_species in species_dict:    species_dict[this_species]=[this_date, index, []] 
      else: 
        if species_dict[this_species][0] > this_date:    species_dict[this_species][2].append(index)  #stored is already better
        else:   species_dict[this_species]= [ this_date, index, species_dict[this_species][2]+[species_dict[this_species][1]] ]  
    indexes_to_remove=[]
    for this_species in species_dict: indexes_to_remove.extend( species_dict[this_species][2] )
    for index in sorted(indexes_to_remove, reverse=True): assembly_entries.pop(index)
    if indexes_to_remove: message(' Removed {0} assemblies in favor of a newer version for the same species '.format( len(indexes_to_remove) ) ) 
    del species_dict; del indexes_to_remove
  ##########
  ##############

  message(' {0} assembly entries found '.format( len(assembly_entries) ))
  assembly_entries.sort(key=lambda x:x['SpeciesName']+'&'+x['AssemblyAccession'])

  if opt['a']:  #and assembly_entries    implicit
    fields_displayed=default_fields_displayed_assembly
    if opt['a']==1: pass
    elif opt['a']==2:  fields_displayed=assembly_entries[0].keys()  ## assuming all entries have the same fields
    else:              
      if opt['a'].startswith('+'): fields_displayed.extend( opt['a'][1:].split(',') )
      else:                        fields_displayed=opt['a'].split(',')
    max_length_per_fields={}
    for field in fields_displayed:  max_length_per_fields[field]= max([len(str(e[field])) for e in assembly_entries] + [len(field)]) #getting max string length for pretty tabular like output
    if opt['tab']: fields_displayed.sort()
    else:          fields_displayed.sort(key=lambda x:max_length_per_fields[x])

  ### preparing to check files
  if opt['c'] or opt['d']:
    if   opt['c']==1: fields_to_test=  map(lambda x:'ftp:'+x, ftp_file_types)
    elif opt['c']:    fields_to_test=  map(lambda x:'ftp:'+x, opt['c'].split(','))
    else:             fields_to_test=[]  #for -d switch below to work
    if opt['d']:      fields_to_test.extend( [ftt for ftt in map(lambda x:'ftp:'+x, opt['d'].split(',')) if not ftt in fields_to_test] ) 
    for field in fields_to_test: 
      new_field_name= '?'+field.split(':')[1]
      if opt['s']:  new_field_name = '#'+new_field_name[1:]
      max_length_per_fields[new_field_name]=4
      fields_displayed.append(new_field_name)
  ### preparing to download files


  ####### writing header line
  if opt['a']:
    if opt['tab']: header_line = join(fields_displayed, '\t')
    else:          header_line = join([field.ljust(max_length_per_fields[field]) for field in fields_displayed ], ' ')
    write(header_line, 1, how='reverse')

  for assembly_e in assembly_entries:
    if opt['c'] or opt['d']:
      ### checking file presence, setting attribute of assembly_e
      for field in fields_to_test:
        ffile= assembly_e[field]
        new_field_name='?'+field.split(':')[1]
        if opt['s']:  
          new_field_name = '#'+new_field_name[1:]
          ssize=wget_spider(ffile, return_size=1)        
          if ssize!='---': ssize=human_readable_size(ssize)
          assembly_e[new_field_name]= ssize
        else:      assembly_e[new_field_name]= wget_spider(ffile)
        
    if opt['a']:
      ### printing to stdout detailed information for each entry
      if opt['tab']:  summary_line=join([ str(assembly_e[field]) for field in fields_displayed ], '\t')
      else:           summary_line=join([ str(assembly_e[field]).ljust(max_length_per_fields[field]) for field in fields_displayed ], ' ')
      write(summary_line, 1)
    
    if opt['d']: 
      ###### download files!   ################ TODO: add md5sum control
      species_name=assembly_e['SpeciesName']
      types_to_download=  opt['d'].split(',')
      for file_type in types_to_download:
        ftp_path=       assembly_e['ftp:'+file_type]
        check_ftp_path= assembly_e['?'+file_type]

        if ftp_path=='None':            continue
        if check_ftp_path=='---':       continue
        file_destination = ftp_file_to_local_path(ftp_path, species_name, local_master_folder);    file_base_name_destination= base_filename(file_destination)
        species_folder= Folder( abspath( directory_name(file_destination) ) )
        gunzipped_file= file_destination.split('.gz')[0];       file_base_name_gunzipped= base_filename(gunzipped_file)
        file_base_name_link= ftp_file_description[file_type];   file_link= species_folder + file_base_name_link
        lets_unzip= not opt['dl']==2;  lets_link= not opt['dl'];  lets_download=True

        if is_file(file_link):     
          #printerr(file_link, 1)
          if lets_link  and opt['F']: message('REMOVE existing link to replace it: "{0}" '.format(file_link)); bbash('rm '+file_link)
          else:                       lets_download=False; lets_unzip=False; lets_link=False
        if is_file(gunzipped_file):     
          if lets_unzip and opt['F']: message('REMOVE existing file to replace it: "{0}" '.format(gunzipped_file)); bbash('rm '+gunzipped_file)
          else:                       lets_download=False; lets_unzip=False
        if is_file (file_destination):
          if opt['F']:  message('REMOVE existing file.gz to replace it: "{0}" '.format(file_destination)); bbash('rm '+file_destination)          
          else:         lets_download=False


        #printerr( str( (lets_download, lets_link, lets_unzip)), 1)
        actions= [ {True:'DOWNLOAD', False:''}[lets_download], {True:'GUNZIP', False:''}[lets_unzip],  {True:'LINK', False:''}[lets_link] ] 
        if any(actions):  message('{0} {1} to "{4}" '.format(join(actions, ' & ') , file_base_name_link.split('.')[0], file_base_name_destination, species_name, species_folder))  #removed   "{2}"   and also   for species: "{3}"
        if lets_download: bbash('cd {0} && wget   {1}'.format(species_folder, ftp_path))  ## now DOWNLOADING with wget
        if lets_unzip:    bbash('cd {0} && gunzip {1}'.format(species_folder, file_base_name_destination ))  ## now GUNZIPPING with gunzip
        if lets_link:     bbash('cd {0} && ln -s "{1}" {2}'.format(species_folder, file_base_name_gunzipped, file_base_name_link))  ## now GUNZIPPING with gunzip


#### later move these functions somewhere else:          
""" Wrapping default Entrez methods to connect to ncbi to allow network problems and batch requests"""
def esearch(**keyargs):
  """ Generic wrap. You should use term=..  and db=..  and field=..  
  Returns a list of uids (strings) """
  opt=get_MMlib_var('opt');  max_attempts=opt['max_attempts']; sleep_time=opt['sleep_time']; retmax=opt['retmax']
  uid_list=[]; retstart=0
  while not uid_list or (int(searched_g['Count']) > len(uid_list) ):  ## to manage when more results are returned than the allowed number
    if  uid_list:    retstart+=int(searched_g['RetMax'])  ## same as setting it to len(genome_uid_list)
    ### add stuff for network problems
    n_attempt=0 #-1 means success
    while n_attempt < max_attempts:
      try:                
        verbose('ESEARCH.n{3} retstart:{0} retmax:{1} {2}'.format(retstart, retmax, join([k+':'+str(keyargs[k]) for k in keyargs], ' '), n_attempt+1), 1)
        searched_g = Entrez.read( Entrez.esearch(retstart=retstart, retmax=retmax,  **keyargs), validate=False ) 
        break
      except URLError:    sleep(sleep_time); n_attempt+=1;  service('esearch FAILED attempt n{0} trying again in {1}s'.format(n_attempt, sleep_time) )
    if n_attempt == max_attempts: write('esearch ERROR network problem??', 1); raise
    ###### 
    #write(searched_g, 1, how='green')
    if not int(searched_g['Count']):  break
    verbose('ESEARCH-OUTIDS: {0}'.format( join(searched_g['IdList'], ' ')), 1 )
    uid_list.extend( searched_g['IdList'] )
  assert int(searched_g['Count']) == len(uid_list)  ## control to see that code is written properly
  return uid_list

def efetch(**keyargs):
  """ Generic wrap. You should use db=..  and id=..list_of_ids..  
  in here, id MUST Be a list (Entrez.efetch accepts a single id as well)
  The fetch is performed by with rettype='docsum'.
  Returns a list of the objects parsed with Entrez.read; the fields of the objects depend on the db queried
  Normally we expect the same quantity of objects in input and in output, a control is performed.
 """
  opt=get_MMlib_var('opt');  max_attempts=opt['max_attempts']; sleep_time=opt['sleep_time']; retmax=opt['retmax']
  batch_size=int(retmax)
  ### add stuff
  list_out=[]
  if not 'id' in keyargs   or  type(keyargs['id'])!=list:  raise Exception, "efetch ERROR id must be provided and must be a list! keyargs: {0} ".format(keyargs)
  if not keyargs['id']: return []
  for batch_index in range( 1+  (len(keyargs['id'])-1) / batch_size ):
    batch_this_list = keyargs['id'][ batch_size*batch_index : batch_size*(batch_index+1) ]
    this_keyargs= dict(keyargs);   del this_keyargs['id']
    ### add stuff for network problems
    n_attempt=0
    while n_attempt < max_attempts:
      try:              
        verbose('EFETCH.n{2} {1} idlist: {0} '.format(join([b for b in batch_this_list], ' '), join([k+':'+str(keyargs[k]) for k in keyargs if k!='id'], ' '), n_attempt+1 ), 1)
        parsed  =Entrez.read( Entrez.efetch(id=batch_this_list, rettype='docsum', validate=False, **this_keyargs), validate=False )
        break
      except URLError:    sleep(sleep_time); n_attempt+=1; service('efetch FAILED attempt n{0} trying again in {1}s'.format(n_attempt, sleep_time) )
    if n_attempt == max_attempts: write('efetch ERROR network problem??', 1); raise

    if    issubclass(type(parsed),list):    fetched_obj_list = parsed
    elif  'DocumentSummarySet' in parsed:   fetched_obj_list = parsed['DocumentSummarySet']['DocumentSummary']
    else:   raise Exception, 'efetch ERROR I do not know how to handle this object to get a list of results! {0}'.format(parsed)
    if len(fetched_obj_list) != len(batch_this_list): raise Exception, "efetch ERROR some ids could not be found! idlist: {0} \n{1}".format(join([b for b in batch_this_list], ' '), join([str(o) for o in fetched_obj_list], '\n'))
    list_out.extend(fetched_obj_list)
  return list_out
    
def wget_spider(ffile, return_size=False):
  """ given a file path in a ftp server, checks with wget if it exists """
  while True:
    if ffile=='None': return '   '
    wget_command='wget --spider {0} '.format(ffile); verbose(wget_command, 1); log = bash(wget_command)[1]
    splt=log.split()
    if   splt and splt[-1]=='exists.' :                                 
      if not return_size:      return 'web'
      else: 
        for index, sstring in enumerate(splt):
          if sstring=='SIZE': return int(splt[index+3])
        raise Exception, 'couldnt parse size here! \n'+log

    elif len(splt)>3 and join(splt[-4:-1], ' ')=='No such file' :       return '---'    
    else:  
      raise Exception,  'wget_spider ERROR with file: "{0}" \nThis is the wget log:{1}'.format(ffile, log)

def load_md5sum(md5file):
  """ Returns a dictionary with the md5sum codes loaded from the remote ncbi ftp site. 
  E.g.
"""
  

  ###############

### copy-pasted from hurry.filesize:
size_names = [    (1024 ** 5, 'P'),    (1024 ** 4, 'T'),     (1024 ** 3, 'G'),     (1024 ** 2, 'M'),     (1024 ** 1, 'K'),    (1024 ** 0, 'B')    ]
def human_readable_size(bytes, system=size_names):
    """Human-readable file size.
    Using the traditional system, where a factor of 1024 is used::     >>> size(10)    '10B'
    >>> size(20000)    '19K'       >>> size(100000)    '97K'           >>> size(200000)    '195K'        """
    for factor, suffix in system:
        if bytes >= factor:            break
    amount = int(bytes/factor)
    if isinstance(suffix, tuple):
        singular, multiple = suffix
        if amount == 1:    suffix = singular
        else:              suffix = multiple
    return str(amount) + suffix

## temp, to improve
def ftp_file_to_local_path(ftp_path, species_name, local_master_folder):
  """ """
  return local_master_folder.rstrip('/') +'/'+  replace(mask_characters(species_name), ' ', '_') + '/' + base_filename(ftp_path)


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
