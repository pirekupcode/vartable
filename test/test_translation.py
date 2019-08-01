#type: ignore
from vartable.translation import TResult
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature, ExactPosition
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import ExactPosition
from typing import Union, List
from typing_extensions import NewType, Literal
from nose.tools import eq_, ok_
from dataclasses import dataclass

#@dataclass 
#class Seq:
#    _data: str 


ExactPosition = NewType('ExactPosition', int)
''' all one-indexed '''


Strand = Literal[-1, 1]


#@dataclass
#class FeatureLocation:
#    _start: ExactPosition
#    _end:   ExactPosition
#    strand: Strand
#
#@dataclass
#class CompoundLocation:
#    parts:  List[FeatureLocation]
#    location_operator: str
#
#@dataclass
#class SeqFeature:
#    location: Union[CompoundLocation, FeatureLocation]

def translate_one(ref: Seq, cdss: List[SeqFeature], pos: int, alt: str) -> TResult:
    from typing import cast
    return cast(None, TResult)
from vartable.translation import  translate_one, dispatch, handle_multi_alts

def make_simple_location(start: int, end: int, strand: Strand) -> FeatureLocation:
    return FeatureLocation(ExactPosition(start), ExactPosition(end), strand)


def translate_one_test1() -> None:
    seq = Seq('ACTGGCG') # ref @ 4 is G
    location = CompoundLocation( [FeatureLocation(ExactPosition(1), ExactPosition(7), 1), \
            FeatureLocation(ExactPosition(8), ExactPosition(11), strand=1)], 'join') 
    cds = SeqFeature(location=location)

    expected = TResult(position=4, alt='A', codon_position=1, ref_codon='GGC', alt_codon='AGC', in_coding_region=True, \
            ref_aa='G', alt_aa='S', synonymous=False, alt_is_invalid_stop=False)

    actual_result = translate_one(seq, [cds], 4, 'A') 
    
    eq_(expected, actual_result, f"\n\nGiven {seq._data}:\n\nexpected: {expected} \n actual: {actual_result}")
    #assert expected == actual_result, f"Expected {expected} but got {actual_result}"

def translate_one_test2() -> None:
    pass
    #CompoundLocation(ExactPosition(
