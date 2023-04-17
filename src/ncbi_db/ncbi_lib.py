from string import *
import sys
from .MMlib3 import *
from Bio import Entrez
from urllib.error import URLError
from ftplib import FTP
import gzip
import shutil
from time import strptime, sleep
import traceback
import easyterm
from pathlib import Path
Entrez.email = None

ncbi_config=easyterm.commandlineopt.CommandLineOptions({
  'ncbi_taxdb':'',
  'email':'',
  })

def load_ncbi_config():
  #nonlocal ncbi_config
  ncbi_config_file=Path( Path.home()+'/.ncbi_db' )
  if ncbi_config_file.is_file():
    ncbi_config.update(
      easyterm.read_config_file(str(ncbi_config_file))
      )

def save_ncbi_config():
  #nonlocal ncbi_config
  ncbi_config_file=Path( Path.home()+'/.ncbi_db' )
  ncbi_config.write_config_file(str(ncbi_config_file))

#### NOTE some values must be defined in MMlib.opt:
# sleep_time
# retmax
# max_attempts
set_MMlib_var('opt', options( {'retmax':250, 'max_attempts':10, 'sleep_time':5, 'v':0}))  ## this is to make this lib work without MMlib command line calls

""" Wrapping default Entrez methods to connect to ncbi to allow network problems and batch requests"""
def esearch(**keyargs):
  """ Generic wrap. You should use term=..  and db=..  and field=..
  Note: these will be used:  opt['max_attempts']   opt['sleep_time']; opt['retmax']
   where opt is get_MMlib_var('opt')   [you can set it with  set_MMlib_var('opt', opt) ]
  Returns a list of uids (strings) """
  opt=get_MMlib_var('opt');  max_attempts=opt['max_attempts']; sleep_time=opt['sleep_time']; retmax=opt['retmax']
  uid_list=[]; retstart=0
  while not uid_list or (int(searched_g['Count']) > len(uid_list) ):  ## to manage when more results are returned than the allowed number
    if  uid_list:    retstart+=int(searched_g['RetMax'])  ## same as setting it to len(genome_uid_list)
    ### add stuff for network problems
    n_attempt=0 #-1 means success
    while n_attempt < max_attempts:
      try:                
        verbose('ESEARCH.n{3} retstart:{0} retmax:{1} {2}'.format(retstart, retmax, ' '.join([k+':'+str(keyargs[k]) for k in keyargs]), n_attempt+1), 1)
        searched_g = Entrez.read( Entrez.esearch(retstart=retstart, retmax=retmax,  **keyargs), validate=False ) 
        break
      except URLError:    sleep(sleep_time); n_attempt+=1;  service('esearch FAILED attempt n{0} trying again in {1}s'.format(n_attempt, sleep_time) )
    if n_attempt == max_attempts: write('esearch ERROR network problem??', 1); raise
    ###### 
    #write(searched_g, 1, how='green')
    if not int(searched_g['Count']):  break
    verbose('ESEARCH-OUTIDS: {0}'.format( ' '.join(searched_g['IdList'])), 1 )
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
  if not 'id' in keyargs   or  type(keyargs['id'])!=list:  raise Exception("efetch ERROR id must be provided and must be a list! keyargs: {0} ".format(keyargs))
  if not keyargs['id']: return []
  for batch_index in range( 1+  (len(keyargs['id'])-1) // batch_size ):
    batch_this_list = keyargs['id'][ batch_size*batch_index : batch_size*(batch_index+1) ]
    if not batch_this_list: continue # bad design, inherited. sorry!
    #write( (batch_index, len(batch_this_list)), 1, how='blue' )
    this_keyargs= dict(keyargs);   del this_keyargs['id']
    ### add stuff for network problems
    n_attempt=0
    while n_attempt < max_attempts:
      try:              
        verbose('EFETCH.n{2} {1} idlist: {0} '.format(' '.join([b for b in batch_this_list]), ' '.join([k+':'+str(keyargs[k]) for k in keyargs if k!='id']), n_attempt+1 ), 1)
        parsed  =Entrez.read( Entrez.efetch(id=batch_this_list, rettype='docsum', validate=False, **this_keyargs), validate=False )
        break
      except URLError:    sleep(sleep_time); n_attempt+=1; service('efetch FAILED attempt n{0} trying again in {1}s'.format(n_attempt, sleep_time) )
    if n_attempt == max_attempts: write('efetch ERROR network problem??', 1); raise

    if    issubclass(type(parsed),list):    fetched_obj_list = parsed
    elif  'DocumentSummarySet' in parsed:   fetched_obj_list = parsed['DocumentSummarySet']['DocumentSummary']
    else:   raise Exception('efetch ERROR I do not know how to handle this object to get a list of results! {0}'.format(parsed))
    if len(fetched_obj_list) != len(batch_this_list): raise Exception("efetch ERROR some ids could not be found! idlist: {0} \n{1}".format(' '.join([b for b in batch_this_list]), '\n'.join([str(o) for o in fetched_obj_list])))
    list_out.extend(fetched_obj_list)
  return list_out

ftp_handler=None
def get_ftp_handler(url='ftp.ncbi.nlm.nih.gov', force=False):
  global ftp_handler
  if ftp_handler is None or force:
    verbose('FTP: establishing connection with {}'.format(url), 1)
    ftp_handler=FTP(host=url)
    ftp_handler.login()
  return(ftp_handler)

### copy-pasted from hurry.filesize:
size_names = [    (1024 ** 5, 'P'),    (1024 ** 4, 'T'),     (1024 ** 3, 'G'),     (1024 ** 2, 'M'),     (1024 ** 1, 'K'),    (1024 ** 0, 'B')    ]
def human_readable_size(bytes, system=None):   
  """Human-readable file size.
  Using the traditional system, where a factor of 1024 is used::     >>> size(10)    '10B'
  >>> size(20000)    '19K'       >>> size(100000)    '97K'           >>> size(200000)    '195K'        """
  the_size_names=size_names if system is None else system
  for i, (factor, suffix) in enumerate(the_size_names):
    if bytes >= factor:            break
  amount = int(bytes/factor)
  if system is None and amount>=1000:
    amount=1
    suffix=the_size_names[i+1][1]
  return str(amount) + suffix


#def wget_spider(ffile, return_size=False):
#  """ given a file path in a ftp server, checks with wget if it exists """
#  while True:
#    if ffile=='None': return '   '
#    wget_command='wget --spider {0} '.format(ffile); verbose(wget_command, 1); log = bash(wget_command)[1]
#    splt=log.split()
#    if   splt and splt[-1]=='exists.' :                                 
#      if not return_size:      return 'web'
#      else: 
#        for index, sstring in enumerate(splt):
#          if sstring=='SIZE': return int(splt[index+3])
#        raise Exception, 'couldnt parse size here! \n'+log
#
#    elif len(splt)>3 and join(splt[-4:-1], ' ')=='No such file' :       return '---'    
#    else:  
#      raise Exception,  'wget_spider ERROR with file: "{0}" \nThis is the wget log:{1}'.format(ffile, log)
