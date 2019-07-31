from toolz import compose
from itertools import starmap
from functools import partial
from Bio import SeqIO 
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.Seq import Seq
from typing import Union, List, Tuple, Any, Iterator, Sequence, TypeVar, Optional
Location = Union[CompoundLocation, FeatureLocation]
Codon = str
AA = str
Nuc = str
T = TypeVar('T')

def first(seq: Union[Iterator[T], Sequence[T]]) -> Optional[T]:
    '''@Warning: mutates seq'''
    if isinstance(seq, Sequence):
        return None if len(seq) == 0 else seq[0]
    else:
        try:
            return next(seq)
        except StopIteration:
            return None
find = compose(first, filter)

# combine multiple alts later or just
# print a position multiple times.
def translate_one(ref: Seq, cdss: List[SeqFeature], pos: int, alt: str) -> Any:
    cds  = find(partial(get_region, pos), cdss)
    if cds is None:
        return f"No matching CDS"
    loc = cds.location
    alt_seq = ref[:pos] + alt + ref[pos+1:]
     # arguably the absolute_codon_start should be calculated elsewhere 
    absolute_codon_start = pos - (pos % 3)
    if cds is None:
        # should still get the other info, like codon_start 
        ref_codon = ref[absolute_codon_start:absolute_codon_start+3]._data
        alt_codon = alt_seq[absolute_codon_start:absolute_codon_start+3]._data
        #return JunkResult(absolute_codon_start, ref_codon, alt_codon)
        return (absolute_codon_start, ref_codon, alt_codon)
    if isinstance(loc, FeatureLocation):
        region_pos = pos - loc._start
    else:
        preceding_region_size =   sum( (x._end - x._start for x in loc.parts  if pos > x._end ) )
        region_pos = pos - loc.parts[0]._start - preceding_region_size
    region_codon_start = region_pos - (region_pos % 3)
    protein_pos =  region_pos // 3
    ref_codon, ref_aa = codon_protein(cds, ref, region_codon_start, protein_pos)
    alt_codon, alt_aa = codon_protein(cds, alt_seq, region_codon_start, protein_pos)
    return TResult(pos, alt, absolute_codon_start, ref_codon, alt_codon, True, \
                   ref_aa, alt_aa, (ref_aa == alt_aa), False)
    #return absolute_codon_start, ref_codon, ref_aa, alt_codon, alt_aa, (ref_aa == alt_aa)

    position:                     int
    alt:                          Nuc

    # codon_position could be determined wihout the ref.
    codon_position:               int

    ref_codon:                    Codon
    alt_codon:                    Codon
    in_coding_region:             bool

    # the beow are optional because of the case where 
    # the alt is not in the coding region.
    ref_aa:                       Optional[AA]
    alt_aa:                       Optional[AA]
    synonymous:                   Optional[bool]
    alt_is_invalid_stop:          Optional[bool]

#TODO: one-based indexing     
def codon_protein(cds: SeqFeature, seq: Seq, region_codon_start: int, protein_pos: int) -> Tuple[Codon, AA]: 
    #TODO: biopython can raise all kinds of exceptions
    region = cds.extract(seq)
    codon = region[region_codon_start:region_codon_start+3]._data
    translation = cds.translate(seq)
    aa = translation[protein_pos]
    return codon, aa
    
def get_region(pos: int, cds: SeqFeature) -> SeqFeature:
    loc = cds.location
    if isinstance(loc, CompoundLocation):
        return cds if find(lambda r: r._start <= pos <= r._end, loc.parts) else None
    else: 
        return cds if loc._start <= pos <= loc._end else None

def dispatch(rec: SeqFeature, variants): 
    #TODO: returns something that's useful with the other vartable components and can be joined. 
    # a @dataclass hierarchy is appropriate for each component: OG vartable, bam-readcount, and translation
    cdss  = list(filter(lambda x: x.type == 'CDS', rec.features))
    results = list(starmap(partial(translate_one, rec.seq, cdss), variants))
    return results

def main():
    #TODO: fill in
    variants = _
    recs = list(SeqIO.parse(open(gb), 'genbank'))
    assert len(recs) == 1, "expected exactly one genbank record in file {gb}"
    results = dispatch(recs[0], variants)
    #TODO: do something with results. mainly, join with other results



from dataclasses import dataclass
@dataclass
class TResult: 
    '''There are four kinds of results: 
       1. An alternate which is outside a coding region.
       2. An alternate that results in an invalid stop codon.
       3. A normal alternate within a CDS.
       4. Some other error state.  
    '''
    # these first two fields are duplicating input.
    position:                     int
    alt:                          Nuc

    # codon_position could be determined wihout the ref.
    codon_position:               int

    ref_codon:                    Codon
    alt_codon:                    Codon
    in_coding_region:             bool

    # the beow are optional because of the case where 
    # the alt is not in the coding region.
    ref_aa:                       Optional[AA]
    alt_aa:                       Optional[AA]
    synonymous:                   Optional[bool]
    alt_is_invalid_stop:          Optional[bool]


from abc import ABCMeta, abstractmethod

#TODO: it's also possible for the ref to be an invalid stop codon.
# In this case, we should throw an error
# Note: @abstractmethod doesn't play well with mypy, 
@dataclass
class BaseResult(metaclass=ABCMeta):
    position: int
    codon_position: int
    ref_codon: Codon
    alt_codon: Codon

    @property
    def within_cds(self) -> bool: ...

# Note that alt_aa or ref_aa could both 
# (and it would have to be both)
# be valid stop codons
class CodingResult(BaseResult):
    ref_aa:                    AA
    alt_aa:                    AA
    synonymous:                bool
    within_cds = True

class OutsideCDSResult(BaseResult): 
    within_cds = False

class InvalidStopAlt(BaseResult):
    ref_aa: AA
    within_cds = True

class InvalidReference(BaseResult):
    error: str
    within_cds = True

#TODO:  Handle the case of multiple alts.

#TODO:  Handle the following error cases below.
'''
gb = 'testdata/adeno.gb'
# 'z' is an invalid codon
invalid_codon_variants = [ (579, 'z'), (579, 'T')]
# 587 as T gives a stop codon
stop_codon_variants = [ (587, 'T'), (579, 'T')]
# TODO:  TranslationError: Codon 'AZG' is invalid
# TODO:  some stop-codon error: TranslationError: Extra in frame stop codon found.
# KeyError: 'TGA'
# TODO: so some kind of exception regex for the above
junk_regoin_variants = [ (570, 'T'), (579, 'T')]
if __name__ == '__main__':
    variants = [(579, 'T')]
    gb = 'testdata/adeno.gb'
    recs = list(SeqIO.parse(open(gb), 'genbank'))
    assert len(recs) == 1, "expected exactly one genbank record in file {gb}"
    results = dispatch(recs[0], variants)
    print(results)
'''

    
    
#TODO: Consider further error states, i.e.
'''
data Degen =  Insert Codon Index
            | WithN Codon Index
            | StopCodon AA Index Codon [Index]
            | Synonymous' AA Index Codon [Index]
            | NonSynonymous [AA] Index Codon [Index] -- Codon index or AA Index? Should make newtypes
            | NormalCodon
'''
