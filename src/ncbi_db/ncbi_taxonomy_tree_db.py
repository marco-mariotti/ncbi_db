#!/usr/bin/env python2.7
import sys,os
import sqlite3
import argparse
from .MMlib3 import bash, printerr
from .ncbi_lib import load_ncbi_config, save_ncbi_config
from pathlib import Path

def main():

  userfile=str(Path.home())+'/.ncbi_taxdb'

  parser = argparse.ArgumentParser(description='Generate a tax.db for ncbi_taxonomy_tree.\nWith no arguments, the program downloads the current "taxdump" from ncbi taxonomy and generates the tax.db in the current directory. It also writes a ~/.ncbi_taxdb file so that later ncbi_taxonomy_tree will automatically find tax.db', prog='ncbi_taxonomy_tree -makedb')
  parser.add_argument("-f","--folder",  help="folder containing names.dmp and nodes.dmp files from ncbi. tax.db will be created there")
  parser.add_argument("-o","--output",  help="folder where tax.db will be created. Also, if --folder is not provided, files will be downloaded here")  
  parser.add_argument("-s","--skip",  help="skip saving the DB location as default (in file ~/.ncbi_config)")
  #parser.add_argument("-w","--overwrite",  help="overwrite without prompting the .ncbi_taxdb file in your home folder, if existing")

  args = parser.parse_args()

  download_files=False
  if not args.folder:
    args.folder='./'
    download_files=True
  if not args.output:
    args.output=args.folder

  print(f'Working directory: {args.output}')
  os.chdir(args.output)
    
  if download_files:
      print( 'Downloading taxdump.tar.gz...', end='', flush=True)
      b = bash('wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz')
      if b[0]: raise Exception(b[1])
      b = bash('wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5')
      if b[0]: raise Exception(b[1])
      md5 = [ x.rstrip() for x in open('taxdump.tar.gz.md5') ][0]
      b = bash('md5sum taxdump.tar.gz')
      if b[0]: raise Exception(b[1])
      if not md5 == b[1]:
          raise Exception("ERROR: something went wrong with wget, md5sums do not match: %s != %s" % (md5,b[1]))
      b=bash('gunzip -c taxdump.tar.gz | tar xf -')
      if b[0]: raise Exception(b[1])
      args.folder = '.'
      print( 'DONE!')
      


  args.folder = args.folder.rstrip('/') +'/'

  names_dict = {}
  for i in open(args.folder + 'names.dmp'):
      tax_id, name_txt, unique_name, name_class = i.rstrip('\t|\n').split('\t|\t')
      names_dict.setdefault(tax_id, {}).setdefault(name_class, name_txt)

  nodes_dict = {}
  for i in open(args.folder +'nodes.dmp'):
      tax_id, parent_tax_id, rank= i.split('\t|\t')[:3]
      nodes_dict.setdefault(tax_id, {}).setdefault('parent_tax_id', parent_tax_id)
      nodes_dict.setdefault(tax_id, {}).setdefault('rank', rank)

  print(f"Preparing to write file: {args.output.rstrip('/')+'/tax.db'}")

  print("I'm going to insert %s taxids into tax.db, this may take a while\n" % len(nodes_dict))
  con = sqlite3.connect(args.output + 'tax.db')
  with con:
      cur = con.cursor()    
      cur.execute("DROP TABLE IF EXISTS species")
      cur.execute("CREATE TABLE species (taxid INT PRIMARY KEY, parent INT, rank VARCHAR(50), name VARCHAR(50))")

      print( 'INSERTING...', sep=' ', flush=True )
      for taxid in nodes_dict:
          parent = nodes_dict[ taxid ][ 'parent_tax_id' ]
          rank   = nodes_dict[ taxid ][ 'rank' ]
          name   = names_dict[ taxid ][ 'scientific name' ]
          parent,rank,name = [x.replace('"','') for x in [parent,rank,name]]
          cur.execute('INSERT INTO species VALUES('+ taxid +', '+ parent +', "'+ rank +'", "'+ name +'")')

      print( 'DONE!\n' )


  print('File tax.db created!')
      
  if not args.skip:
    ncbi_config=load_ncbi_config()
    taxdb_path=os.path.abspath(str(args.output) + 'tax.db')
    ncbi_config['ncbi_taxdb']= taxdb_path
    printerr(f"Setting ncbi_taxdb = {taxdb_path} into the user configuration --> ~/.ncbi_config", 1)
    save_ncbi_config()
    
    # if os.path.isfile(userfile) and not args.overwrite:
    #     old_content=open(userfile).readline().strip()
    #     if old_content!=content:
    #         a=None
    #         while not a in ['', 'Y', 'N']:
    #             a=input(f'\nExisting: {userfile} with content: {old_content}\nOverwriting with new content: {content}\nConfirm? (Y or Enter to proceed, N to cancel) > ')
    #         if a=='N':
    #             sys.exit(1)
    # with open(userfile, 'w') as fw:
    #     print(f'Writing user config file --> {userfile}')
    #     fw.write(content+'\n')

if __name__ == "__main__":
    main()
