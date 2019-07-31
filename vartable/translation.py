from toolz import compose
from itertools import starmap
from functools import partial
from Bio import SeqIO 
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.Seq import Seq
from typing import Union, List, Tuple, Any, Iterator, Sequence, TypeVar, Optional, Iterable
from Bio.Data.CodonTable import TranslationError
from dataclasses import dataclass
import re

Location = Union[CompoundLocation, FeatureLocation]
Codon = str
AA = str
Nuc = str
T = TypeVar('T')

#TODO: one-based indexing check
#TODO:  Handle the case of multiple alts.

@dataclass
class TResult: 
    '''Object to represent output
    There are four kinds of results: 
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
    alt_seq = ref[:pos] + alt + ref[pos+1:]
     # arguably the absolute_codon_start should be calculated elsewhere 
    absolute_codon_start = pos - (pos % 3)
    if cds is None:
        # should still get the other info, like codon_start 
        ref_codon = ref[absolute_codon_start:absolute_codon_start+3]._data
        alt_codon = alt_seq[absolute_codon_start:absolute_codon_start+3]._data
        return TResult(pos, alt, absolute_codon_start, ref_codon, alt_codon, False, None, None, None, None)
    if isinstance(cds.location, FeatureLocation):
        region_pos = pos - cds.location._start
    else:
        preceding_region_size =   sum( (x._end - x._start for x in cds.location.parts  if pos > x._end ) )
        region_pos = pos - cds.location.parts[0]._start - preceding_region_size
    protein_pos =  region_pos // 3
    ref_codon = get_codon(cds, ref, region_pos)
    alt_codon = get_codon(cds, alt_seq, region_pos)
    ref_aa = translate(cds, ref, protein_pos) 
    try:
        alt_aa = translate(cds, alt_seq, protein_pos)
    except TranslationError as e:
        if re.compile(r'stop\scodon\sfound', re.I).search(str(e)):
           return TResult(pos, alt, absolute_codon_start, \
                               ref_codon, alt_codon, True, ref_aa, None, None, True)
        else:
            raise e
    return TResult(pos, alt, absolute_codon_start, ref_codon, alt_codon, True, \
                   ref_aa, alt_aa, (ref_aa == alt_aa), False)

def get_codon(cds: SeqFeature, seq: Seq, region_pos: int) -> Codon:
    region_codon_start = region_pos - (region_pos % 3)
    region = cds.extract(seq)
    codon = region[region_codon_start:region_codon_start+3]._data
    return codon

def translate(cds: SeqFeature, seq: Seq, protein_pos: int) -> Tuple[Codon, AA]: 
    #TODO: biopython can raise all kinds of exceptions
    translation = cds.translate(seq)
    aa = translation[protein_pos]
    return aa
    
def get_region(pos: int, cds: SeqFeature) -> SeqFeature:
    loc = cds.location
    if isinstance(loc, CompoundLocation):
        return cds if find(lambda r: r._start <= pos <= r._end, loc.parts) else None
    else: 
        return cds if loc._start <= pos <= loc._end else None

def dispatch(rec: SeqFeature, variants: Iterable[Tuple[int, Nuc]]) -> List[TResult]: 
    #TODO: returns something that's useful with the other vartable components and can be joined. 
    # a @dataclass hierarchy could be appropriate for each component: OG vartable, bam-readcount, and translation
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





#TODO: Remove or use the following
#
#  Unused Classes
#####################################
from abc import ABCMeta, abstractmethod

# Note that alt_aa or ref_aa could both 
# (and it would have to be both)
# be valid stop codons
# Note: @abstractmethod doesn't play well with mypy, so not using it here

@dataclass
class BaseResult(metaclass=ABCMeta):
    position: int
    codon_position: int
    ref_codon: Codon
    alt_codon: Codon

    @property
    def within_cds(self) -> bool: ...

# Good result
class CodingResult(BaseResult):
    ref_aa:                    AA
    alt_aa:                    AA
    synonymous:                bool
    within_cds = True

# Abnormal results
class OutsideCDSResult(BaseResult): 
    within_cds = False

class InvalidStopAlt(BaseResult):
    ref_aa: AA
    within_cds = True

class InvalidReference(BaseResult):
    error: str
    within_cds = True

    
#TODO: Consider further error states, i.e.
'''
data Degen =  Insert Codon Index
            | WithN Codon Index
            | StopCodon AA Index Codon [Index]
            | Synonymous' AA Index Codon [Index]
            | NonSynonymous [AA] Index Codon [Index] -- Codon index or AA Index? Should make newtypes
            | NormalCodon
'''
#TODO: turn these into tests.
'''
if __name__ == '__main__':
    good_variants = [(579, 'T')]
    gb = 'testdata/adeno.gb'
    recs = list(SeqIO.parse(open(gb), 'genbank'))
    assert len(recs) == 1, "expected exactly one genbank record in file {gb}"
    results = dispatch(recs[0], good_variants)
    print(results)
    junk_regoin_variants = [ (570, 'T'), (579, 'T')]
    print(dispatch(recs[0], junk_regoin_variants))
    # 587 as T gives a stop codon
    # TranslationError: Extra in frame stop codon found.
    stop_codon_variants = [ (587, 'T'), (579, 'T')]
    print(dispatch(recs[0], stop_codon_variants))
    # KeyError: 'TGA'
    # TranslationError: Codon 'AZG' is invalid
    # The below should be the only one to raise an exception.
    invalid_codon_variants = [ (579, 'z'), (579, 'T')]
    print(dispatch(recs[0], invalid_codon_variants))
'''
