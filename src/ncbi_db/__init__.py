#__all__=['__version__']
import functools
from ._version import __version__

from .ncbi_assembly import main as run_ncbi_assembly
from .ncbi_sequences import main as run_ncbi_sequences
#from .ncbi_db_info import main as run_ncbi_db_info
from .ncbi_pubmed import main as run_ncbi_pubmed
from .ncbi_search import main as run_ncbi_search
from .ncbi_taxonomy import main as run_ncbi_taxonomy
from .ncbi_taxonomy_tree import main as run_ncbi_taxonomy_tree
from .parse_genbank import main as run_parse_genbank

def suppress_return(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        f(*args, **kwargs)
        return None
    return wrapper

run_ncbi_assembly=suppress_return(run_ncbi_assembly)
run_ncbi_sequences=suppress_return(run_ncbi_sequences)
run_ncbi_pubmed=suppress_return(run_ncbi_pubmed)
run_ncbi_search=suppress_return(run_ncbi_search)
run_ncbi_taxonomy=suppress_return(run_ncbi_taxonomy)
run_ncbi_taxonomy_tree=suppress_return(run_ncbi_taxonomy_tree)
run_parse_genbank=suppress_return(run_parse_genbank)

