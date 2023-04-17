#!/usr/bin/env python3
import sys,os
import json
import sqlite3
import argparse
from collections import defaultdict
from .MMlib3 import bbash, unmask_characters, mask_characters, get_taxids_from_ncbi_db, printerr, random_folder
import shutil, os
from easyterm import NoTracebackError
from .ncbi_taxonomy_tree_db import main as run_ncbi_taxonomy_tree_db
from .ncbi_lib import load_ncbi_config, save_ncbi_config
from pathlib import Path

import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)
from ete3 import Tree

def get_parser():
    """returns the argparse parser object with the argument option requirements"""

    parser     = argparse.ArgumentParser(description='Get the phylogenetic tree from a list of species (based on NCBI taxonomy)')
    taxid_file = parser.add_argument_group('input')
    tree       = parser.add_argument_group('output')
    render     = parser.add_argument_group('render')
    database   = parser.add_argument_group('database')

    parser.add_argument('-temp', help='temporal folder', default='/tmp/')
    parser.add_argument('-q','--quiet', help='do not show warnings', action="store_true", default=False)

#    taxid_file.add_mutually_exclusive_group(required = True)
    taxid_file.add_argument("-n", "--name",
                        help="list of scientific names, one name per line. Only scientific names from NCBI Taxonomy (not synonyms!)\
                              are recognized")
    taxid_file.add_argument("-t", "--taxid",
                        help="list of taxonomy ids, one taxid per line. The scientific name in NCBI associated to the taxid is used\
                              in the tree")
    taxid_file.add_argument("-u", "--unmask",
                        help="allow masked characters in the input file", action="store_true")
    taxid_file.add_argument("-m", "--mask",
                        help="include masked characters in the tree", action="store_true")
    taxid_file.add_argument("-tn","--taxid_name",
                        help="two columns input file: taxid<TAB>species_name. species_name is used in the tree")
    taxid_file.add_argument("-ta","--taxid_annotation",
                        help="include any attribute(s) as a feature in any of the nodes of the tree. The input file format is\
                              'taxid=value[;key=value]' ( taxid key is required). The node corresponding to the specified taxid\
                               (including internal nodes) will include each key=value as a feature")
    taxid_file.add_argument("-na","--name_annotation",
                        help="same as --taxid_annotation, but including name instead of taxid: 'name=value[;key=value]'")
    tree.add_argument("-f","--fix_internal",
                      help="internal nodes present in the input list will be added a sister leaf node with the same name and taxid",
                      action='store_true')
    tree.add_argument("-N","--normalize",
                      help="normalize the number of ranks.", action='store_true')
    tree.add_argument("-l","--ladderize",
                      help="sort the branches of a given tree (swapping children nodes) according to the size of each partition.", action='store_true')
    tree.add_argument("-um","--ultrametric",
                      help="converts a tree to ultrametric topology (all leaves must have the same distance to root).", action='store_true')
    tree.add_argument("-nw","--newick",
                      help="output tree in newick format")
    tree.add_argument("-js","--json",
                      help="output tree in json format")
    tree.add_argument("-x",
                      help="output file with taxids list")
    tree.add_argument("--lineage",
                      help="get the full lineage for each species. Each terminal node includes a 'lineage' attribute", action="store_true")
    tree.add_argument("-nn","--no_name",
                      help="do not include species name in the leaves, use taxid instead", action="store_true")
    tree.add_argument("-d", "--non_dicotomic",
                      help="keep non_dicotomic nodes", action="store_false")
    tree.add_argument("-s", "--show_tree",
                      help="show interactive tree", action="store_true")
    tree.add_argument("-p", "--print_stdout",
                      help="print tree to stdout", action="store_true")
    
    render.add_argument("-r","--render",
                        help="render tree as image. File extension will determine the image format (.PDF, .SVG or .PNG)")
    render.add_argument("-H","--height",
                        help="height of the image in units",
                        type=float,
                        default=None)
    render.add_argument("-W","--width",
                        help="width of the image in units",
                        type=float,
                        default=None)

    render.add_argument("-U","--units",
                        help="px: pixels (default), mm: millimeters, in: inches",
                        default='px')
    render.add_argument("-DPI",
                        help="dots per inches (300)",
                        type=float,
                        default=300)
    database.add_argument('-db', '--database',
                          help="tax.db file. See --makedb to obtain it. If not specified, reads path from ~/.ncbi_taxdb if existing",
                          default=None)
    database.add_argument('-mdb', '--makedb',
                        help="Create a tax.db file. For a dedicated help page, run: ncbi_taxonomy_tree --makedb -h",
                          default=None)

    return parser


def parse_opts(x):
    """to be used when the program is imported as a module:
import ncbi_taxonomy_tree
args = ncbi_taxonomy_tree.parse_opts('-y your -o options -f for -t the -p program')
ncbi_taxonomy_tree.main(args)
"""
    parser = get_parser()                         
    args = parser.parse_args(x.split())
    return args


def get_taxid_from_names(input_names, temp_dir, unmask, database, is_list=False, quiet=False):
    """ gets a file with the list of names, and returns the list of taxids and a dict {taxid:name}"""
    if not is_list:
        names_list = set()
        for line in open(input_names):
            if unmask:
                names_list.add( unmask_characters(line.strip().replace('_',' ')) )
            else:
                names_list.add(line.strip())
    else:
        if not input_names: raise Exception('ERROR if is_list=True, arg is mandatory')
        names_list = set(input_names)
    print_stderr( '%s unique names in input\n' % (len(names_list)), quiet)
    names_db = '/'.join( database.split('/')[:-1] + ['names.dmp'] )
    if not os.path.exists(names_db):
        raise Exception('ERROR file not found: '+ names_db)
    out_hash, ambigous_hash, not_found_hash = get_taxids_from_ncbi_db(names_list, ncbi_db=names_db, temp_dir=temp_dir, full=True)
    if ambigous_hash: print_stderr(  'WARNING: %s ambiguous names: %s\n' % (len(ambigous_hash),' '.join( ambigous_hash )) , quiet)
    if not_found_hash: print_stderr( 'WARNING: %s names not found in ncbi: %s\n' % (len(not_found_hash), ','.join(list(not_found_hash.keys()))) , quiet)
    # check if there is a correspondence 1:1 between names and taxids
    taxids_from_input_dict={}
    for k,v in list(out_hash.items()):
        taxids_from_input_dict.setdefault(v,[]).append(k)
    for x,v in list(taxids_from_input_dict.items()):
        if len(v)>1:
            print_stderr( 'WARNING: %s elements in your input returned the same taxid, only the first one will be included in the tree: %s\n' % (len(v),', '.join(v)) , quiet)
    return out_hash,taxids_from_input_dict

def get_taxid_and_annotation():
    """ """
    taxid_list = set()
    if not os.path.exists(args.taxid_annotation): raise Exception('file not found: %s\n' % args.taxid_annotation)
    infile = args.taxid_annotation
    taxid_attrs = {}
    for line in open(infile):
        if not 'taxid=' in line:
            raise Exception('ERROR: taxid key is missing in '+ line.strip() +'. Please, use "taxid=value[;key=value]" in the input file.') 
        attrs_dict = {}
        for attr in line.strip().split(';'):
            k,v = attr.split('=')
            if 'taxid' in k:
                taxid = int(v)
                continue
            attrs_dict[ k ] = v
        taxid_attrs[ taxid ] = attrs_dict
        taxid_list.add(taxid)
    return list(taxid_list), taxid_attrs

def get_name_and_annotation(arg, temp_dir):
    """ """
    if not os.path.exists(arg): raise Exception('file not found: %s\n' % arg)
    name_attrs = {}
    input_names = []
    for line in open(arg):
        if not 'name=' in line:
            raise Exception('ERROR: taxid key is missing in "'+ line.strip() +'". Please, use "name=value[;key=value]" in the input file.')
        attrs_dict = {}
        for attr in line.strip().split(';'):
            k,v = attr.split('=')
            if 'name' in k:
                if args.unmask:
                    name = unmask_characters(v.replace('_',' '))
                else:
                    name = v
                continue
            attrs_dict[ k ] = v
        name_attrs[ name ] = attrs_dict
        input_names.append(name)
    out_hash, taxids_from_input_dict = get_taxid_from_names(input_names, temp_dir, args.unmask, args.database, is_list=True, quiet=args.quiet)
    taxid_list = list(out_hash.values())
    # map taxid and attributes
    taxid_attrs = {}
    for name,taxid in list(out_hash.items()):
            taxid_attrs[ taxid ] = name_attrs[ name ]
    return taxid_list, taxids_from_input_dict, taxid_attrs

def node_lineage(node, t):
    """DEPRECATED! returns the lineage of node as a list with the ancestral nodes. each element in the list corresponds to (distance, rank, name), where distance is the number of nodes between
 node and the root of the tree."""
    dist = int(node.get_distance(t, topology_only=True))
    lineage = [ (dist,node.rank,node.name,node.taxid) ]
    parent = node.up
    while parent and not parent.name == 'Root':
        dist = int(parent.get_distance(t, topology_only=True))
        lineage.append( (dist,parent.rank,parent.name,parent.taxid) )
        parent = parent.up
    return sorted( lineage, key = lambda x:x[0] )

def include_lineage(t):
    """ """
    lineages_dict={}
    for node in t.get_descendants('preorder'):
        if not node.is_leaf():
            for name in node.get_leaf_names():
#                lineages_dict.setdefault(name,[]).append( (node.rank,node.name,node.taxid) )
                lineages_dict.setdefault(name,[]).append( node.name )
        else:
            node.add_feature('lineage', '||'.join( lineages_dict[node.name] + [node.name] )   )


def print_stderr(msg, quiet=False):
    if not quiet:
        sys.stderr.write(msg)


def main(args=None):

    if args is None:

        if '--makedb' in sys.argv or '-mdb' in sys.argv:
            sys.argv = [x for x in sys.argv if not x in ['--makedb', '-mdb']]
            run_ncbi_taxonomy_tree_db()
            sys.exit()

        parser = get_parser()
        args = parser.parse_args()

    ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


    # ncbi database file
    ncbi_config=load_ncbi_config()    
    if args.database is None:
        # userfile=str(Path.home())+'/.ncbi_taxdb'
        # if os.path.exists(userfile):
        #     with open(userfile) as fr:
        #         content=fr.readline().strip()
        #args.database=content.split('ncbi_taxdb')[1].split('=')[1].strip()
        args.database=ncbi_config['ncbi_taxdb']
        if not args.database:
            args.database=None #just to ensure type is consistent
    elif not args.database is None and not ncbi_config['ncbi_taxdb']:
        taxdb_path=os.path.abspath(args.database)
        ncbi_config['ncbi_taxdb']=taxdb_path
        printerr(f"Setting ncbi_taxdb = {taxdb_path} into the user configuration --> ~/.ncbi_config", 1)        
        save_ncbi_config()
        
    if args.database is None or not os.path.exists(args.database):
        raise NoTracebackError('ERROR: --database file not found: %s\nPerhaps you never built a NCBI tax.db? For help, run: ncbi_taxonomy_tree --makedb -h' % args.database)

    temp_folder=random_folder(args.temp).rstrip('/')+'/'
    
    #### INPUT
    # if input file contains taxonomy ids
    if args.taxid:
        if not os.path.exists(args.taxid): raise Exception('file not found: %s\n' % args.taxid)
        taxid_list = list(set([x.rstrip() for x in open(args.taxid)]))
        try:
            [int(x) for x in taxid_list ]
        except:
            raise Exception('ERROR: this does not look like a taxonomy id: %s; try --name instead\n' % x)

    # if input file contains names
    elif args.name:
        if not os.path.exists(args.name): raise Exception('file not found: %s\n' % args.name)
        out_hash,taxids_from_input_dict = get_taxid_from_names(args.name, temp_folder, args.unmask, args.database, quiet=args.quiet)
        taxid_list = list(out_hash.values())

    # if both taxid and name in input
    elif args.taxid_name:
        if not os.path.exists(args.taxid_name): raise Exception('file not found: %s\n' % args.taxid_name)
        infile = args.taxid_name
        taxid_name_dict = {}
        taxid_list = set()
        for line in open(infile):
            taxid, name = line.strip().split("\t")
            taxid_name_dict[ taxid ] = name
            taxid_list.add(taxid)

    # if an annotation is provided for the taxids
    elif args.taxid_annotation:
        taxid_list, taxid_attrs = get_taxid_and_annotation()

    # if an annotation is provided for the names
    elif args.name_annotation:
        taxid_list, taxids_from_input_dict, taxid_attrs = get_name_and_annotation(args.name_annotation, temp_folder)

    # list of taxids is ready for building the tree
    if not taxid_list:
        sys.stderr.write('no taxids to build the tree, probably your input does not match with ncbi names,\
 please try --unmask option, in case your input names are masked\n')
        sys.exit()

    set_taxid=set(taxid_list)
    taxid_list = list(set_taxid)


    #### BUILD
    print_stderr( 'input for tree: ' + str(len(taxid_list)) + ' taxids\n' , args.quiet)
    print_stderr( 'building tree...', args.quiet)
    nodes_dict = {}
    taxids_remove = []
    cursor = sqlite3.connect(args.database).cursor()
    t, nodes_dict, taxids_remove = buildTree(taxid_list, nodes_dict, taxids_remove, cursor)
    t.add_feature('rank', 'Root')
    print_stderr( 'DONE!\n\n', args.quiet)
    if taxids_remove:
        print_stderr('WARNING: the following ids were not found in the taxonomy database used: '+ ','.join(map(str,taxids_remove)) +'\n', args.quiet)

    for n in t.traverse():
        n._requested= n.taxid in set_taxid
        
    ### add control: if requested, do not trim 
        
    #### OUTPUT
    # get the full lineage of species, this is done before removing non-dicotomic nodes, so all nodes are presenent (no need for -d option)
    if args.lineage:
        include_lineage(t)


    # taxonomy ids file
    if args.x:
        open(args.x, 'w').write('\n'.join( map(str, taxid_list) )+'\n')
        print_stderr( str(len(taxid_list)) +' taxids stored in file: '+ args.x +'\n', args.quiet)            

    # normalize ranks
    if args.normalize:
        t=normalize_ranks(t, ranks)

    # remove non_dicotomic nodes
    if args.non_dicotomic and not args.normalize:
        if len(t.get_children()) == 1:
            t.get_children()[0].delete(prevent_nondicotomic=False)
        for node in t.get_descendants():
#            if len(node.get_children()) == 1 and not node.taxid in taxid_list:
            if len(node.get_children()) == 1 and not node.taxid in set_taxid:            
#                node.delete(prevent_nondicotomic=True)
                node.delete(prevent_nondicotomic=False)                

                
    # include node annotation
    if args.taxid_annotation or args.name_annotation:
        for taxid, attrs_dict in list(taxid_attrs.items()):
            try:
                node = t.search_nodes(taxid=taxid)[0]
                node.add_features(**attrs_dict)
            except:
                # problems detected so far:
                # two taxids in the list point to the same node, then only one is present in the tree and the other one is missing (although it was not include in the taxids_remove list because it returned a result)
                print_stderr("WARNING: "+str(taxid)+' not found in the tree, its attributes could not be included\n', args.quiet)

    # json outfile (node.name is still the taxid)
    if args.json:
        if args.ladderize:
            t.ladderize()
        jsoned = ete2json(t)
        with open(args.json, 'w') as json_outfile:
            json.dump(jsoned, json_outfile)
        print_stderr('Tree saved in json format: '+ args.json +'\n', args.quiet)


    # include species name
    internal_nodes_in_input=[]
    if not args.no_name:
        if args.name:
            for taxid,names_list in list(taxids_from_input_dict.items()):
                nodes = t.search_nodes(taxid = taxid)
                if len(nodes)>1: raise Exception('ERROR: more than one node with taxid: %s' % taxid)
                if not nodes:  raise Exception('ERROR: no node has taxid: %s' % taxid)
                nodes[0].add_features(name=names_list[0])
                if not nodes[0].is_leaf(): internal_nodes_in_input.append( nodes[0] )
        elif args.taxid_name:
            for leaf in t.get_leaves():
                leaf.name = taxid_name_dict[ leaf.taxid ]
        else:
            for leaf in t.get_leaves():
                leaf.name = leaf.name


    if args.fix_internal:# and internal_nodes_in_input:
        for taxid in taxid_list:
            
            nodes = t.search_nodes(taxid = taxid)
            if len(nodes)>1: raise Exception('ERROR: more than one node with taxid: %s' % taxid)
            if not nodes:
                sys.stderr.write( 'WARNING: no node has taxid: %s\n' % taxid )
                continue
            if nodes[0].is_leaf(): continue
            node = nodes[0]
            print_stderr('MESSAGE: internal node duplicated as a leaf: %s\n' % node.name , args.quiet)
            new_node = node.up.add_child(name='artificial_node')
            new_node.add_feature('rank','no rank')
            d = node.detach()
            new_node.add_child(d)
            new_leaf = new_node.add_child(name=node.name)
            new_leaf.add_features( taxid=node.taxid, rank='no rank' ) #scientific_name=node.scientific_name,
            if len(node.get_children()) == 1: node.delete(prevent_nondicotomic=False)



            
    print_stderr( '\n'+ str(len(taxid_list)) +' taxids in input\n'+ str(len(t)) +' leaves in tree\n', args.quiet)


    if args.mask:
        for leaf in t.get_leaves():
            leaf.name = mask_characters( leaf.name ).replace(' ','_')
        
    if args.ladderize:
        t.ladderize()

    if args.ultrametric:
        t.convert_to_ultrametric(t.get_farthest_node()[1])

    # newick outfile
    if args.newick:
        t.write(outfile=args.newick, features=[])
        print_stderr('Tree saved in newick format: '+ args.newick +'\n')
    
    if args.render:
        t.render(args.render, w=args.width, h=args.height, units=args.units, dpi = args.DPI)
    
    if args.show_tree:
        t.show()

    if args.print_stdout:
        print(t)

    if os.path.isdir(temp_folder):
        shutil.rmtree(temp_folder)
    return t

# END OF PROGRAM
# FUNCTIONS
def normalize_ranks(t,ranks):
    normtree = Tree()
    for l in t.get_leaves():
        out = []
        lineage={}
        taxids={}
        parent = l.up
        while parent:
            if parent.rank in ranks:
                lineage[parent.rank]=parent.scientific_name
                taxids[parent.rank]=parent.taxid
            parent = parent.up

        for rank in ranks:
            if not rank in lineage:
                lineage[rank] = 'No '+rank
                taxids[rank]= 'No taxid'

        for i in range(len(ranks)):
            rank=ranks[i]
            if i==0: # superkingdom
                if not normtree.search_nodes(name=lineage[rank], rank=rank):
                    child = normtree.add_child(name=lineage[rank])
#                    child.add_features(rank=rank,lineage=lineage[rank],taxid=taxids[rank])
                    child.add_feature('rank', rank )
                    child.add_feature('lineage', lineage[rank])
                    child.add_feature('taxid', taxids[rank])
            else:
                if rank == 'species' and lineage[rank] == 'No species':
                    lineage['species'] = l.name
                parent_node_lineage =':'.join([ lineage[x] for x in ranks[:i]])
                node_lineage =':'.join([ lineage[x] for x in ranks[:i] + [rank]])
                parent = normtree.search_nodes(lineage=parent_node_lineage)
                if not parent or len(parent)>1:
                    raise
                if not normtree.search_nodes(name=lineage[rank], lineage=node_lineage):
                    child = parent[0].add_child(name=lineage[rank])
#                    child.add_features(rank=rank,lineage=lineage[rank],taxid=taxids[rank])
                    child.add_feature('rank', rank )
                    child.add_feature('lineage', node_lineage)
                    child.add_feature('taxid', taxids[rank])
                    if rank == 'species':
                        child.scientific_name = l.scientific_name
                        child.taxid = l.taxid
    return normtree


#def node_lineage(node):
#    """returns lineage of node"""
#    lineage = [ node.rank ]
#    parent = node.up
#    while parent and not parent.name == '1':
#        lineage.insert(0,parent.rank)
#        parent = parent.up
#    return lineage

def collapse_species(t, taxid_list, ranks):
    """collapse internal nodes if they are present in input list.
this is meant to solve the situation in which a taxid has children taxids also present in the list.
Once solved, the resulting tree has each taxid from the list as a leaf.
It is inconsistennt with their lineages, but this step is necessary for its representation in the sunburst plot.
"""
    stop = False
    for node in t.get_descendants():
        if not node.is_leaf() and node.name in taxid_list:
            
            if node.rank in ranks:
                print_stderr('ERROR: taxid "'+ node.name +'" is present in infile and has children also present in infile.\nIt can not be collapsed because its rank "'+ node.rank +'" is among the items in ranks list '+ str(ranks) +'. Please, remove it from the infile.\n')
                print_stderr( str(len(list(set(node.get_leaf_names()).intersection(taxid_list)))) +' leaves names under this node: ' + str(list(set(node.get_leaf_names()).intersection(taxid_list))) +'\n', args.quiet)
                stop = True
                
            for child in node.get_descendants('postorder'):
                if child.name in taxid_list:
                    detached = child.detach()
                    node.up.add_child(detached)
                else:
                    child.detach()

    if stop: # this tag allows to identify all the problematic nodes in one run
        sys.exit()

#def normalize_ranks(t, input_ranks):
#    """WORKING ON IT...
#gets a tree and normalizes the number of ranks.
#the returned tree will contain only the ranks in the ranks list,
#the same number of ranks for all the species.
#"""
#    sys.stderr.write('MESSAGE: normalizing ranks in tree to ' + ','.join(input_ranks) +'\n')
#    ranks = ['root'] + input_ranks[:-1]
#    rank_node_down_dict = {} # {rank:rank_down} 
#    for r in range(len(ranks)-1):
#        rank_node_down_dict[ ranks[r] ] = ranks[r+1]
#
#    t.rank = 'root'
#
#    for node in t.get_descendants('postorder'):
#        if node.is_leaf():
#            ## Normalize leaves
#            ## I want the parent of the leaf to be the last element in the ranks list
#            ## so I need to collapse everything under that rank
#            
#            lineage = node_lineage(node)
#
#            ## check if duplicate ranks in lineage
#            lineage_copy = lineage[:]
#            while 'no rank' in lineage_copy:
#                lineage_copy.remove('no rank')
#            if not len(lineage_copy) == len(set(lineage_copy)):
#                sys.stderr.write('ERROR: lineage with duplicated ranks. Please, remove taxid: '+ node.name + ' from input file.\nInconsistent lineage: '+ str(lineage) +'\n')
#                sys.exit()
#
#            if not node.up.rank == ranks[-1]:
#                ## parent rank is not last item in ranks list
#                parent = node.up
#                if ranks[-1] in lineage:
#                    ## last item in ranks list is present in lineage
#                    ## delete nodes until rank of parent is ranks[-1]
#                    while not parent.rank == ranks[-1]:
#                        parent.delete(prevent_nondicotomic=False)
#                        parent = node.up
#
#                else:
#                    ## ranks[-1] is not present in lineage
#                    ## I will reconstruct the lineage
#
#                    ## nodes with rank not present in ranks will be removed
#                    ## until the rank of is present in ranks
#                    while not parent.rank in ranks:
#                        child = parent
#                        parent = parent.up
#                        child.delete(prevent_nondicotomic=False)
#                        
#                    while not parent.rank == ranks[-1]:
#                        new = parent.add_child(name='und')
#                        new.rank = rank_node_down_dict[ parent.rank ]
#                        new.scientific_name = ''
#                        child = node.detach()
#                        new.add_child(child)
#                        parent = new
#        else:
#            ## node is not leaf
#            ## I will simply delete the node if its rank is not in ranks
#            if not node.rank in ranks:
#                node.delete(prevent_nondicotomic=False)
#                
#    ## include missing ranks
#    for node in t.get_descendants('levelorder'):
#        if not node.is_leaf():
#            while not node.rank == rank_node_down_dict[ node.up.rank ]:
#                new = node.up.add_child(name='und')
#                new.rank = rank_node_down_dict[ node.up.rank ]
#                new.scientific_name = ''
#                child = node.detach()
#                new.add_child(child)
#                
#    for node in t.get_leaves():
#        node.rank = 'species'
#        if not node_lineage(node) == input_ranks:
#            sys.stderr.write('WARNING: There is a problem with taxid: '+ node.name +'\n' )
#            sys.stderr.write('         lineage != input_ranks...\n')
#            sys.stderr.write('         lineage: '+ str(node_lineage(node)) +'\n')
#            sys.stderr.write('         input_ranks: '+ str(input_ranks) +'\n')


def annotate_tree(t,infile):
    """Annotates tree from infile (key=value[;key=value]). taxid key required. node.name must be the taxid."""
    taxid_attrs = defaultdict()
    for line in open(infile):
        if not 'taxid=' in line:
            sys.stderr.write('taxid key is missing in '+ line.strip() +'. Please, use "taxid=value;" in the input file.\n')
            sys.exit()
        else:
            attrs_dict = {}
            for attr in line.strip().split(';'):
                k,v = attr.split('=')
                attrs_dict[ k ] = v
            taxid_attrs[ attrs_dict['taxid'] ] = attrs_dict
            
    for taxid, attrs_dict in list(taxid_attrs.items()):
        try:        node = t&taxid
        except:
            sys.stderr.write("ERROR "+str(taxid)+' not found! ')
            raise 
        node.add_features(**attrs_dict)

    
    
def ete2json(t, jsoned = dict()):
    """Recursive function to convert Tree into serializable json object"""
    if not t.is_leaf():
        jsoned[ 'children' ] = []
        for node in t.get_children():
            children = {}
            children['taxid'] = node.name
            children['rank']  = node.rank
            children['name']  = node.scientific_name
            children = ete2json(node, children)
            jsoned['children'].append(children)
    else:
        jsoned['size'] = '1'
    return jsoned

def query_a_list(taxid_list, cursor):
    """gets a list of taxids, splits it in batches of size 999,
and does the query to the ncbi database local file.
Returns a list where each item has (taxid, parent_taxid, rank, name)"""
    n = 999
    results = []
    for i in range(0, len(taxid_list), n):
        batch = taxid_list[i:i+n]
        placeholders = ', '.join(['?']*len(batch))
        cur_query = 'SELECT * FROM species WHERE taxid IN (%s)' % placeholders
        cursor.execute(cur_query, batch)
        results += cursor.fetchall()
    return results


def buildTree(taxid_list, nodes_dict, taxids_remove, cursor):
    """Recursive function, returns a ete tree object from a list of taxids.
Requires a cursor connected to a sqlite db build using the script /users/rg/didac/NCBI/Taxonomy/update_sqlite_DB.py
nodes_dict is an empty dict
taxids_remove is an empty list """
    results = query_a_list(taxid_list, cursor)
    
    # check if all taxids returned a result
    if len(set(taxid_list)) != len(results):
        taxids_with_result = set([ x[0] for x in results])
        taxids_remove += list(set(map(int, taxid_list)) - taxids_with_result )

    parent_taxid_list = []
    for result in results:
        taxid, parent_taxid, rank, name = result
        parent_taxid_list.append(parent_taxid)

        if not taxid in nodes_dict:
            c = Tree()
            c.add_feature('name', name)
            nodes_dict[ taxid ] = c

        # I don't have scientific name and rank for parent_taxid yet, but next iteration it will be the taxid
        nodes_dict[ taxid ].add_features(name=name, taxid=taxid, rank=rank)
        # add child to node parent_taxid
        if not parent_taxid in nodes_dict:
            p = Tree()
            p.add_feature('taxid', parent_taxid)
            p.add_child( nodes_dict[ taxid ] )
            nodes_dict[ parent_taxid ] = p

        else:
            # check if taxid is a child of parent_taxid (already in nodes_dict), otherwise adding it
            for descendant in nodes_dict[ parent_taxid ].iter_descendants():
                if taxid == descendant.taxid:
                    break
            else:
                nodes_dict[ parent_taxid ].add_child( nodes_dict[ taxid ] )

    parent_taxid_list = list(set(parent_taxid_list))

    try:
        # "1" is the root of the NCBI tree, if "1" is in parent_taxid_list, and it will become an empty list inside this try
        parent_taxid_list.remove(1)
    except:
        pass

    if parent_taxid_list:
        t,nodes_dict,taxids_remove = buildTree(parent_taxid_list, nodes_dict, taxids_remove, cursor)
    else:
        nodes_dict[ 1 ].add_features(name='Root', rank='Root')
        return nodes_dict[ 1 ], nodes_dict, taxids_remove
    return t, nodes_dict, taxids_remove

      
if __name__ == '__main__':
    main()

