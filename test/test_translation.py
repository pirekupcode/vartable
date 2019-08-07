#type: ignore
from vartable.translation import TResult, translate_one
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature, ExactPosition
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import ExactPosition
from typing import Union, List
from typing_extensions import NewType, Literal
from nose.tools import eq_, ok_
from dataclasses import dataclass
# Warning: all tests should use 0-based indexing.

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

def make_simple_location(start: int, end: int, strand: Strand) -> FeatureLocation:
    return FeatureLocation(ExactPosition(start), ExactPosition(end), strand)


def test_translate_one1() -> None:
    seq = Seq('ACTGGCG') # ref @ 4 is G
    location = CompoundLocation( [FeatureLocation(ExactPosition(0), ExactPosition(6), 1), \
            FeatureLocation(ExactPosition(8), ExactPosition(11), strand=1)], 'join') 
    cds = SeqFeature(location=location)

    expected = TResult(position=3, alt='A', codon_position=3, ref_codon='GGC', alt_codon='AGC', in_coding_region=True, \
            ref_aa='G', alt_aa='S', synonymous=False, alt_is_invalid_stop=False)

    actual_result = translate_one(seq, [cds], 3, 'A') 
    
    eq_(expected, actual_result, f"\n\nGiven {seq._data}:\n\nexpected: {expected} \n actual: {actual_result}")

def translate_one_test2() -> None:
    pass
    #CompoundLocation(ExactPosition(
