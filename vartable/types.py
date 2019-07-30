from typing import Union, NamedTuple, List, Dict, Tuple, Callable, TypeVar, Iterator, Optional, Callable, Iterable, Sequence
from typing_extensions import NewType, Literal
from abc import ABCMeta, ABC

T = TypeVar('T') 
Nuc = Literal['A', 'C', 'G', 'T', 'N'] 

AA = Literal['A', 'R', 'N', 'D', 'C', 'E', \
        'Q', 'G', 'H', 'I', 'L', 'K', 'M', \
        'F', 'P', 'S', 'T', 'W', 'Y', 'V']

Ref = NewType('Ref', str) 

Codon = NewType('Codon', str)# Codon = Tuple[Nuc, Nuc, Nuc]

#TODO: check math of one-based indexing

#class CDS(metaclass=ABCMeta):
#class CDS(ABC):
#    pass
class SimpleCDS:
    '''note: uses 1-based indexing.'''
    def __init__(self, is_reverse: bool, start: int, end: int):
        self.isReverse = is_reverse
        self.start = start
        self.end = end 

class CompoundCDS:
    def __init__(self, is_reverse: bool, regions: List[SimpleCDS]):
        self.regions = regions
        self.isReverse = is_reverse


CDS = Union[SimpleCDS, CompoundCDS]

class Result:
    '''note: uses 1-based indexing.'''
    def __init__(self, refCodon: Codon, altCodon: Codon, refProtein: AA, altProtein: AA, isSynoynymous: bool, codonPosition: int, alt: Nuc):
        self.altCodon           = altCodon
        self.refCodon           = refCodon
        self.refProtein      = refProtein
        self.altProtein      = altProtein
        self.isSynonymous   = isSynoynymous
        self.codonPostion    = codonPosition
        self.alt             = alt

def first(seq: Union[Iterator[T], Sequence[T]]) -> Optional[T]:
    '''@Warning: mutates seq'''
    if isinstance(seq, Sequence):
        return None if len(seq) == 0 else seq[0]
    else:
        try:
            return next(seq)
        except StopIteration:
            return None
from typing import NoReturn, Any

def assert_never(x: object) -> NoReturn:
    assert False, "Unhandled type: {}".format(type(x).__name__)

def find(p: Callable[[T], bool], seq: Iterable[T]) -> Optional[T]:
    return first(filter(p, seq))

def codon_pos(pos: int) -> Optional[int]:
    positions = (i*3 for i in range(999_9999_999)) # infinite generator
    return find(lambda x: x >= pos, positions)

def compound_translate(ref: Ref, cds: CompoundCDS, alt: Nuc, pos: int) -> Result:
    '''return the effected codon, the resulting protein for the reference and alternate translation, and whether synonymous'''
    base_region = find(lambda r: r.start <= pos <= r.end, cds.regions) 
    if base_region is None:
        raise ValueError(f"Bad CompoundIndex {cds}") 
    base_index = cds.regions.index(base_region) 
    if (pos - base_region.start < 3):
        assert base_index > 0, f"looking for region after {base_region} for pos {pos} but it was the first in the list in {cds}"
        before_region: Optional[SimpleCDS] = cds.regions[base_index - 1]
    else:
        before_region = None 
    if (base_region.end - pos < 3):
        assert base_index < (len(cds.regions) - 1), f"looking for region after {base_region} for pos {pos} but it was the last in the list in {cds}"
        after_region: Optional[SimpleCDS] = cds.regions[base_index + 1]
    else:
        after_region = None 
    regions = filter(None, [before_region, base_region, after_region]) # this works! (typechecks)
    slices  = [ref[r.start:r.end][::-1] if r.isReverse else \
                                ref[r.start:r.end] for r in regions]
    full_region: str = ''.join(slices)
    return simple_translate(ref, SimpleCDS(cds.isReverse, 1, len(full_region)), alt, pos)

def simple_translate(ref: Ref, cds: SimpleCDS, alt: Nuc, alt_pos: int) -> Result:
    '''return the effected codon, the resulting protein for the reference and alternate translation, and whether synonymous'''
    #NOTE: this function doesn't have to worry about slicing usinsg the CDS, because we assume that we've got the below right
    # possibly we should use pyristent with invariants
    # rather than conflating data checing with processing
    assert cds.start <= alt_pos <= cds.end - 3
    assert len(ref) >= (cds.end - cds.start)
    cpos = codon_pos(alt_pos)
    if cpos is None:
        raise ValueError(f"Codon of position {alt_pos} not found in reference.")
    return codons_and_proteins(ref, cpos, alt_pos, alt, cds.isReverse)

def in_region(pos: int, cds: CDS) -> bool:
    if isinstance(cds, CompoundCDS):
        return find(lambda r: r.start <= pos <= r.end, cds.regions) is not None
    else:
        return cds.start <= pos <= cds.end

def codons_and_proteins(ref: Ref, cpos: int, alt_pos: int, alt: Nuc, reverse: bool) -> Result:
    codon = ref[cpos:cpos+3]
    alcp = cpos - alt_pos
    alt_codon = codon[:alcp] + alt + codon[alcp+1:]
    if reverse:
      ref_aa = table.get(codon[::-1], None)
      alt_aa = table.get(alt_codon[::-1], None)
    else:
      ref_aa = table.get(codon, None)
      alt_aa = table.get(alt_codon, None) 
    if ref_aa is None or alt_aa is None:
        #TODO: Shouldn't error here. This is just a junk translation.
        raise ValueError("codons_and_proteins, {locals()}")
    return Result( Codon(alt_codon), Codon(codon), 
            ref_aa, alt_aa, (alt_aa == ref_aa), cpos, alt)

from typing import NamedTuple
class Failure(NamedTuple):
    exception: Exception

#TODO: design error-handling ahead of time and read. 
from typing_extensions import TypedDict


# some way to do this betterer?
def dispatch_translate(ref: Ref, cds: CDS, pos: int, alt: Nuc) -> Result:
    if isinstance(cds, SimpleCDS):
        return simple_translate(ref, cds, alt, pos)
    elif isinstance(cds, CompoundCDS):
        return compound_translate(ref, cds, alt, pos)

def try2fail(f: Callable[...,T], *args: Any) -> Union[T, Failure]:
    try: 
        return f(*args)
    except Exception as e:
        return Failure(e)

def none2fail(msg: str, x: Optional[T]) -> Union[T, Failure]:
    if x is None:
        return Failure(ValueError(msg))
    else:
        return x

ResDict = TypedDict('ResDict', 
        { 'position' : int, 'alts' : List[Nuc], 'codons' : List[Codon], 
          'proteins' : List[AA], 'synonymous' : bool })


def catMaybes(xs: List[Optional[T]]) -> Optional[List[T]]:
    res: List[T] = []
    for x in xs:
        if x is not None:
            res.append(x)
        else:
            return None
    return res
V = TypeVar('V')
E = TypeVar('E')

def catEither(xs: List[Union[E, V]]) -> Union[E, List[V]]:
    res: List[V] = []
    for x in xs:
        if isinstance(x, V):
            res.append(x)
        else:
            return None
    return res
        
def translate_all(ref: Ref, cdss: List[CDS], variants: Dict[int, List[Nuc]]) -> List[Union[Failure, ResDict]]:
    ''' Dict[int, List[Nuc]] are our alternates organized by their position
        note that we may have multiple alternates at a single position.
        In this case, altProtein becomes the ''' 
    cds_alts: List[Tuple[Optional[CDS], Tuple[int, List[Nuc]]]]  = \
            [( find(lambda cds: in_region(pos, cds), cdss),  \
            (pos, variants[pos]) ) \
            for pos in variants] 
    def dispatch(ref: Ref, cds: Optional[CDS], var: Tuple[int, List[Nuc]]) -> Union[Failure, ResDict]:
        pos, nts = var
        if cds is None:
            return Failure(ValueError(f"No CDS found at position {pos}.")) 
        #results = [dispatch_translate(ref, cds, pos, nt) for nt in nts]
        results_ = [try2fail(dispatch_translate, ref, cds, pos, nt) for nt in nts]
        results = catEither(results_)
        reveal_type(results)
        assert bool(results) # non-empty
        dicts = []
        for r in results: 
            if isinstance(r, Failure):
                dicts.append(r)
            else:
                pass
        alts = [x.alt for x in results]
        refProteins = [x.refProtein for x in results]
        proteins = [ refProteins[0] ] + [x.altProtein for x in results]
        codons = [results[0].refCodon] + [x.altCodon for x in results]
        from functools import reduce
        synonymous = reduce(lambda x,y: x and y, [x.isSynonymous for x in results])
        return { 'position' : pos, 'alts' : nts, 'codons' : codons, \
                'proteins' : proteins, 'synonymous' : synonymous}
    from itertools import starmap
    return [dispatch(ref, x[0], x[1]) for x in cds_alts]
        # return Failure(ValueError("foo")) 

    #filter(lambda x: x is not None, 
#    cds_alts: Dict[CDS, Tuple[int, List[Nuc]]] = \

    pass

table: Dict[str, Optional[AA]] = { 
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
    'TAC':'Y', 'TAT':'Y', 'TAA':None, 'TAG':None, 
    'TGC':'C', 'TGT':'C', 'TGA':None, 'TGG':'W', 
} 

AALong = Literal[ 'Ala', 	'Arg', 	'Asn', 	'Asp', 	'Cys', 	'Glu', 	'Gln', 	\
        'Gly', 	'His', 	'Ile', 	'Leu', 	'Lys', 	'Met', 	'Phe', 	\
        'Pro', 	'Ser', 	'Thr', 	'Trp', 	'Tyr', 	'Val'] 	

#def simple_translate(ref: Ref, cds: SimpleCDS, alts: List[Nuc], alt_pos: int) -> Result:
#def get_mult()    
#    codon = ref[cpos:cpos+3]
#    ref_protein = table.get(codon, None)
#    alcp = cpos - alt_pos
#    def make_codon(nt: Nuc) -> str:
#        return ref[:alcp] + nt + ref[alcp+1:]
#    alt_codons = [make_codon(b) for b in alts]
#    alt_proteins = [table.get(cdn, None) for cdn in alt_codons]
