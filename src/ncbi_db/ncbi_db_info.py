#! /usr/bin/env python3
from .ncbi_lib import *
from .MMlib3 import options
help_msg=""" Displays information about ncbi db that can be queried through the Entrez module. Usage:

$  ncbi_db_info    db_name

To get a list of possible db_name, run with no options

### Options:
-email         email address to use with NCBI Entrez. Required for first usage of ncbi online tools
-h OR --help   print this help and exit"""

command_line_synonyms={}

def_opt= { 'email':'',
'i':0, }


#########################################################
###### start main program function

def main(args=None):
#########################################################
############ loading options
  global opt

  if args is None:
    opt=command_line(def_opt, help_msg, 'i', synonyms=command_line_synonyms )
  else:
    opt=options(def_opt)
    opt.update(args)

  email_setup(opt['email'])
  
  set_MMlib_var('opt', opt)
  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input
  input_name= opt['i']

  if not input_name:
    ## no input specified, let's see the list of dbs
    write('These is the list of possible DBs available at ncbi:', 1)
    write(' '+   '\n '.join (  sorted(Entrez.read( Entrez.einfo() )['DbList'])), 1  )

  else:
    write( 'Name'.ljust(5)+' '+ 'FullName'.ljust(35)  +' '+  'Description', 1)  
    for k in sorted(Entrez.read( Entrez.einfo(db=input_name) )['DbInfo']['FieldList'], key=lambda x:x['FullName']):
      write( k['Name'].ljust(5)+' '+ k['FullName'].ljust(35)  +' '+  k['Description'], 1)


  ###############



#######################################################################################################################################

def close_program():
  # if 'temp_folder' in globals() and is_directory(temp_folder):
  #   bbash('rm -r '+temp_folder)
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
