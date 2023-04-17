#! /usr/bin/env python3
from .ncbi_lib import *


help_msg="""Interrogates NCBI online through Bio.Entrez, retrieves info (and files) related to the genomic assemblies for an organism or lineage. 

## Usage:   $ ncbi_assembly  -S species_or_lineage   [options]

############## Options ############## 
### Which entries
-G        query the 'genome' db, instead of 'assembly' directly; this gets fewer entries but better quality
 # if -G is active:
 -g       show detailed  information for each genome entry. Syntax just like -a (no additional fields available though)
 -z       display also the entries with AssemblyID=0; normally they are skipped
 # if -G is NOT active:
 -e       do not merge different versions of the same assembly (normally only the latest version is displayed)
 -n       do not keep only the newest assembly per species

### Alternative query methods
 -tax     provide a ncbi taxid; in this way you don't have to provide option -S. Note: only exact matches are found, the descendants are not reported (unlike with -S)
 -aa      provide assembly accessions (style: GCA_001940725.1); optionally multiple, comma-separated. Overrides -S and -tax

### Information displayed
-a              shows detailed information for each assembly entry. As arguments:
     0           -> shows nothing (quiet mode)
     1 [default] -> shows a few built-in selected fields
     2           -> shows all the non-null fields in this object
     f1,f2,..    -> shows this list of comma-separated fields 
     +f1,f2,..   -> add these fields to the default fields
 Note: available fields are those found in the NCBI xml; additionally, these are made available here:
   props      all strings in PropertyList joined together
   date       NCBIReleaseDate made shorter
   ftp        ftp folder for this assembly, derived from the string in the Meta attribute
   ftp:XXX    derived path of a specific file.  Here's the possible XXX file classes:
 | dna   genome file   (*_genomic.fna)  | gff   annotation file    (*_genomic.gff)
 | pep   proteome file (*_protein.faa)  | fea   feature table file (*_feature_table.txt)
 | gbf   genbank format annoation file  | md5  md5checksum file

### Check or download files
-c  [x1,x2..]   test if files are found locally or at NCBI ftp. Possible x values: see XXX file classes above (Default: all)
 ## if option -c  is active:
 -s             show file size (compressed) instead of just presence
-d  [x1,x2..]   download certain files; syntax just as -c. Md5sum is checked
 ## if option -d  is active:
 -f             download master folder. Files are downloaded in a subfolder named after the species, then uncompressed and linked
 -dl            do not link files with a standard name (e.g. genome.fa for dna); if used as "-dl 2", do not even uncompress with gunzip
 -F             force download, even if files are found in the local master folder 

### Technical
-retmax         max items requested at each connection to ncbi
-max_attempts   max attempts querying a database
-sleep_time     time between attempts in seconds
-max_down       max attempts at (re-)downloading a file if md5sum does not match
-k              skip and keep running if error occurs for a certain assembly / species, instead of crashing

### Other options
-email        email address to use with NCBI Entrez. Required for first usage of ncbi online tools
-tab          tab separated output. The main messages then go to stderr, to allow redirection of the tables to an output file
-markup       comma separated terminal markup colors to show results in output. Default: 'yellow,,magenta,red'
-v            verbose output, for debugging
-print_opt    print currently active options
-h OR --help  print this help and exit"""

command_line_synonyms={}

def_opt= { 'temp':'/tmp/', 
           'S':'',   'tax':False, 'aa':0,
           'z':0, 'e':0, 'n':0,
           'tab':0, 
           'c':0, 's':0, 'd':0, 
           'f':'', 'F':0, 'dl':0, 
           'G':0,   
           'g':1, 'a':1,
           'k':0, 'max_down':4,
           'retmax':250, 'max_attempts':10, 'sleep_time':5, 
           'v':0,
           'markup':'yellow,,magenta,red',
           'email': '',
}

ftp_file_description={'dna':'genome.fa',  'pep':'proteome.fa', 
'gff':'annotation.gff', 'fea':'feature_table.txt', 'gbf':'genbank.gbff',
'md5':'md5sum.txt'}
ftp_file_types=list(ftp_file_description.keys())

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
  if not args:
    opt=command_line(def_opt, help_msg, 'S', synonyms=command_line_synonyms )
  else:
    opt=args
    for k in def_opt:
      if k not in opt: opt[k]=def_opt[k]


  set_MMlib_var('opt', opt)
  global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input

  ### checking for illegal options
  species=opt['S']
  n_query_options=sum([int(bool(opt[x])) for x in ['tax','S','aa']])
  if not n_query_options:  raise Exception("""ERROR you must specify one input query among:
 -S    species or lineage name 
 -tax  NCBI taxonomy taxid
 -aa   assembly accession
See --help""")
  if  n_query_options>1:
    raise Exception("ERROR options -S, -tax and -aa are mutually exclusive. See --help")
  if (opt['tax'] or opt['aa']) and opt['G']:         raise Exception("ERROR option -G cannot be used with -aa or -tax. See --help")
  if opt['n'] and opt['G']:           raise Exception("ERROR option -n makes no sense with -G. See --help")
  if opt['e'] and opt['G']:           raise Exception("ERROR option -e makes no sense with -G. See --help")
  if opt['s'] and not opt['c']:       raise Exception("ERROR option -s makes sense only if -c is active. See --help")
  if opt['f']:
    local_master_folder=Folder(opt['f'])
    test_writeable_folder(local_master_folder)
  else: local_master_folder='./'
  if opt['d']: 
    if opt['d']==1 or any([not c in ftp_file_types for c in opt['d'].split(',')]): raise Exception("ERROR illegal argument for option -d ! See -help")
    if not opt['f']:                  raise Exception("ERROR to download (option -d) you must provide a local folder with option -f. See --help")    
  if opt['c'] and opt['c']!=1 and any(  [not c in ftp_file_types for c in opt['c'].split(',') ] ): raise Exception("ERROR illegal argument for option -c ! See -help")
  markup=opt['markup'].split(',')
  
  ######### config: email
  email_setup(opt['email'])
  
  # ncbi_config=load_ncbi_config()
  # if not opt['email'] and not ncbi_config['ncbi_email']:
  #   raise Exception('ERROR on your first usage, you must provide an email address to use NCBI services with -email')

  # if opt['email'] and opt['email']!=ncbi_config['ncbi_email']:
  #   printerr(f"Setting email = {opt['email']} into the user configuration --> ~/.ncbi_config", 1)
  #   ncbi_config['ncbi_email']=opt['email']    
  #   try:
  #     save_ncbi_config()
  #   except:
  #     printerr('WARNING: could not save configuration to ~/.ncbi_config', 1)


  
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
      elif opt['g']==2:  fields_displayed=list(genome_entries[0].keys())  ## assuming all entries have the same fields
      else:              
        if opt['g'].startswith('+'): fields_displayed.extend( opt['g'][1:].split(',') )
        else:                        fields_displayed=opt['g'].split(',')
      max_length_per_fields={}
      for field in fields_displayed:  max_length_per_fields[field]= max([len(str(e[field])) for e in genome_entries] + [len(field)]) #getting max string length for pretty tabular like output
      if opt['tab']:  fields_displayed.sort()   ##sorting alphabetical
      else:           fields_displayed.sort(key=lambda x:max_length_per_fields[x])
      if opt['tab']:
        header_line = '\t'.join(fields_displayed)
        write(header_line, 1, how='reverse')
      else:
        for i,field in enumerate(fields_displayed):
          write(field.ljust(max_length_per_fields[field])+' ', 0,  how=  markup[i%len(markup)] )
        write('', 1)

    #### cycle for printing information, and also storing assembly ids
    for genome_e in genome_entries:
      if not genome_e['AssemblyID']=='0': assembly_uid_list.append( genome_e['AssemblyID'] )
      else: 
        discarded_entries+=1
        if not opt['z'] and genome_e['AssemblyID']=='0': continue  ## skipping strange entries
      if opt['g']:
        if opt['tab']:
          summary_line='\t'.join([ str(genome_e[field]) for field in fields_displayed ]) 
          write(summary_line, 1)
        else:
          for i,field in enumerate(fields_displayed):
            write(str(genome_e[field]).ljust(max_length_per_fields[field])+' ', 0,  how=  markup[i%len(markup)] )
          write('', 1)
            #summary_line=' '.join([ str(genome_e[field]).ljust(max_length_per_fields[field]) for field in fields_displayed ])

    if discarded_entries:   message(' {0} entries were discarded, since missing an AssemblyID '.format( discarded_entries ) ) 
 
  elif opt['aa']:
    assembly_uid_list=[]
    for assembly_accession in opt['aa'].split(','):
      assembly_uid_list+=esearch(db='assembly', term=assembly_accession, field='AssemblyAccession')      
  elif opt['tax']:
    taxid=str(opt['tax'])
    assembly_uid_list= esearch(db='assembly', term=taxid, field='Taxonomy ID')  
  else:  ### querying assembly db directly   
    ######## DEFAULT BEHAVIOUR, with -S
    assembly_uid_list= esearch(db='assembly', term=species, field='Organism')  


  #write(assembly_uid_list, 1, how='yellow')
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
    assembly_e['props']= ' '.join(assembly_e['PropertyList'])
    # print '\n'.join(sorted(assembly_e.keys()))
    assembly_e['date']= assembly_e['SeqReleaseDate'].split()[0] #was NCBIReleaseDate

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
    elif opt['a']==2:  fields_displayed=list(assembly_entries[0].keys())  ## assuming all entries have the same fields
    else:              
      if opt['a'].startswith('+'): fields_displayed.extend( opt['a'][1:].split(',') )
      else:                        fields_displayed=opt['a'].split(',')
    max_length_per_fields={}
    for field in fields_displayed:  max_length_per_fields[field]= max([len(str(e[field])) for e in assembly_entries] + [len(field)]) #getting max string length for pretty tabular like output
    if opt['tab']: fields_displayed.sort()
    else:          fields_displayed.sort(key=lambda x:max_length_per_fields[x])

  ### preparing to check files
  if opt['c'] or opt['d']:
    if   opt['c']==1: fields_to_test=  ['ftp:'+x for x in ftp_file_types]
    elif opt['c']:    fields_to_test=  ['ftp:'+x for x in opt['c'].split(',')]
    else:             fields_to_test=[]  #for -d switch below to work
    if opt['d']:
      fields_to_test.extend( [ftt for ftt in ['ftp:'+x for x in opt['d'].split(',')] if not ftt in fields_to_test] )
      if 'ftp:md5' in fields_to_test: fields_to_test.remove('ftp:md5')
      fields_to_test.insert(0, 'ftp:md5')
    for field in fields_to_test: 
      new_field_name= '?'+field.split(':')[1]
      if opt['s']:  new_field_name = '#'+new_field_name[1:]
      max_length_per_fields[new_field_name]=4
      fields_displayed.append(new_field_name)
  ### preparing to download files


  ####### writing header line
  if opt['a']:
    if opt['tab']:
      header_line = '\t'.join(fields_displayed)
      write(header_line, 1, how='reverse')      
    else:
      for i,field in enumerate(fields_displayed):
        write(field.ljust(max_length_per_fields[field])+' ', 0,  how=  markup[i%len(markup)] )
      write('', 1)

  for assembly_e in assembly_entries:
    species_name=assembly_e['SpeciesName']
    try:
    
      if opt['c'] or opt['d']:
        ### checking file presence, setting attribute of assembly_e
        for field in fields_to_test:
          ffile= assembly_e[field]
          ffile_cut=ffile[ 6+ (ffile[6:].find('/')+1): ]
          new_field_name='?'+field.split(':')[1]
          ftp_handler=get_ftp_handler()

          if not ffile =='None':
            try:
              ssize=ftp_handler.size(ffile_cut)
            except Exception as e:
              if 'No such file or directory' in str(e):
                ssize='---'
              else:
                printerr(str(e), 1, how='red')
                ssize='?!?'
                #raise e

          else:
            ssize='---'

          # except:
          #   print ffile_cut
          #   for k in sorted(assembly_e):
          #     write('{} : {}'.format(k, [assembly_e[k]]), 1)


          if opt['s']:  
            new_field_name = '#'+new_field_name[1:]  # changing ?gff to #gff
            #ssize=wget_spider(ffile, return_size=1)        
            if not ssize in (None, '---', '?!?'): ssize=human_readable_size(ssize)
            assembly_e[new_field_name]= ssize
          else:          
            #assembly_e[new_field_name]= wget_spider(ffile)
            assembly_e[new_field_name]= 'web' if not ssize in [None, '---']  else '---'
            file_destination = ftp_file_to_local_path(ffile, species_name, local_master_folder)
            gunzipped_file= file_destination.split('.gz')[0]
            if os.path.isfile(file_destination) or os.path.isfile(gunzipped_file):            assembly_e[new_field_name]='loc'

      if opt['a']:
        ### printing to stdout detailed information for each entry
        if opt['tab']:
          summary_line='\t'.join([ str(assembly_e[field]) for field in fields_displayed ])
          write(summary_line, 1)          
        else:
          summary_line=' '.join([ str(assembly_e[field]).ljust(max_length_per_fields[field]) for field in fields_displayed ])
          for i,field in enumerate(fields_displayed):
            write(str(assembly_e[field]).ljust(max_length_per_fields[field])+' ', 0,  how=  markup[i%len(markup)] )
          write('', 1)

      if opt['d']: 
        ###### download files!   
        types_to_download=  opt['d'].split(',')
        if 'md5' in types_to_download: types_to_download.remove('md5')
        types_to_download.insert(0, 'md5')
        md5sum_hash=None

        for file_type in types_to_download:
          ftp_path=       assembly_e['ftp:'+file_type]            
          check_ftp_path= (assembly_e['?'+file_type]    if ('?'+file_type) in assembly_e else
                           assembly_e['#'+file_type])          

          ftp_path_cut=   ftp_path[ 6+ (ftp_path[6:].find('/')+1): ]

          if ftp_path=='None':            continue
          if check_ftp_path=='---':       continue

          file_destination = ftp_file_to_local_path(ftp_path, species_name, local_master_folder)
          file_base_name_destination= base_filename(file_destination)
          species_folder= Folder( abspath( directory_name(file_destination) ) )
          gunzipped_file= file_destination.split('.gz')[0]
          file_base_name_gunzipped= base_filename(gunzipped_file)
          file_base_name_link= ftp_file_description[file_type]
          file_link= species_folder + file_base_name_link

          lets_unzip= not opt['dl']==2
          lets_link= not opt['dl']
          lets_download=True

          if file_type=='md5':
            lets_unzip=False
            lets_link=False
            local_md5_file=file_destination   #useful for all other files

          if os.path.isfile(file_link):     
            #printerr(file_link, 1)
            if lets_link  and opt['F']: message('REMOVE existing link to replace it: "{0}" '.format(file_link));      os.remove(file_link)
            else:                       lets_download=False; lets_unzip=False; lets_link=False
          if os.path.isfile(gunzipped_file):     
            if lets_unzip and opt['F']: message('REMOVE existing file to replace it: "{0}" '.format(gunzipped_file)); os.remove(gunzipped_file)
            else:                       lets_download=False; lets_unzip=False
          if os.path.isfile(file_destination):
            if opt['F']:  message('REMOVE existing file.gz to replace it: "{0}" '.format(file_destination));          os.remove(file_destination)          
            else:         lets_download=False

          actions= [ {True:'DOWNLOAD,', False:''}[lets_download], {True:'GUNZIP,', False:''}[lets_unzip],  {True:'LINK,', False:''}[lets_link] ] 
          if any(actions):  message('{0} {1} to "{4}" '.format(' '.join(actions).strip().strip(',') , file_base_name_link.split('.')[0], file_base_name_destination, species_name, species_folder))  #removed   "{2}"   and also   for species: "{3}"
          if lets_download:                     
            ftp_handler=get_ftp_handler()
            ftp_handler.set_pasv(False)

            if file_type=='md5':
              file_destination_tmp=file_destination +'.downloading'
              fh = open( '{}'.format(file_destination_tmp), "wb")            #fh = open( '{}'.format(file_destination_tmp), "w")
              ftp_handler.retrbinary("RETR " + ftp_path_cut, fh.write)       #ftp_handler.retrlines("RETR " + ftp_path_cut, fh.write) #, 8*1024)
              fh.close()
              os.rename(file_destination_tmp,  file_destination)
              continue  ## breaking, passing to next file
              
            if md5sum_hash is None:
              ## loading md5sum file
              md5sum_hash={}
              for line in open(local_md5_file):
                md5sum, the_file=line[:-1].split()
                md5sum_hash[the_file.strip('.').strip('/')]=md5sum
            
            target_md5sum=md5sum_hash[file_base_name_destination]

            for download_attempt in range(opt['max_down']):
              file_destination_tmp=file_destination +'.downloading'
              fh = open( '{}'.format(file_destination_tmp), "wb")            #fh = open( '{}'.format(file_destination_tmp), "w")
              ftp_handler.retrbinary("RETR " + ftp_path_cut, fh.write)       #ftp_handler.retrlines("RETR " + ftp_path_cut, fh.write) #, 8*1024)
              fh.close()

              downloaded_md5sum=bbash('md5sum "{}"'.format(file_destination_tmp)).split()[0]
              if downloaded_md5sum!=target_md5sum:
                message('md5sum did not match for {} at download attempt {}'.format(file_base_name_destination, download_attempt+1) )                
              else:
                break  #ok, no more download attempts

            if downloaded_md5sum!=target_md5sum:
              m='ERROR md5sum not matching for {} : local {} vs web {}'.format(file_base_name_destination, downloaded_md5sum, target_md5sum)
              raise Exception(m)            

            os.rename(file_destination_tmp,  file_destination)
            #bbash('cd {0} && wget   {1}'.format(species_folder, ftp_path))  ## now DOWNLOADING with wget

          if lets_unzip:
            gunzipped_file_tmp=gunzipped_file+'.gunzipping'
            with gzip.open(file_destination, 'rb') as f_in:
              with open(gunzipped_file_tmp, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                os.rename(gunzipped_file_tmp,  gunzipped_file)
                os.remove(file_destination)
                #bbash('cd {0} && gunzip {1}'.format(species_folder, file_base_name_destination ))  ## now GUNZIPPING with gunzip

          if lets_link:
            os.symlink( file_base_name_gunzipped,  file_link) 
            #bbash('cd {0} && ln -s "{1}" {2}'.format(species_folder, file_base_name_gunzipped, file_base_name_link))  

    except Exception as e:
      if opt['k']:
        write('ERROR {}'.format(species_name), 1)
        message('ERROR {} {}'.format(species_name,  traceback.format_exc( sys.exc_info()[2]) ))
      else:
        raise
      ftp_handler=get_ftp_handler(force=True)

  ###############
  return assembly_entries

## temp, to improve
def ftp_file_to_local_path(ftp_path, species_name, local_master_folder):
  """ """
  return local_master_folder.rstrip('/') +'/'+  mask_characters(species_name).replace(' ', '_') + '/' + base_filename(ftp_path)


#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('ls '+temp_folder)
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
