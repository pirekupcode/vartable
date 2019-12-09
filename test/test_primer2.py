from typing import *
from nose.tools import eq_, ok_
from vartable.primers import match_primers, Alt, MappingCoordinate, BamRead, PrimerCount, calc_primer_position

def test_match_primers_2() -> None:
    '''Test something'''
    coordinates = [ MappingCoordinate(referenceStart=10, referenceEnd=40) ]
    alts = [ Alt(24, 'G') ]
    mapperReadsOverlappingPosition =  [BamRead(
                       referenceStart=5,
                       referenceEnd=234,
                       queryStart=255,
                       queryEnd=355,
                       querySequence='ACGACGGCGCACACGGCAGCAGCACCCCCGCTTTTTTTATATATGGCACGAGCGCGCAGCGACGCAGCGCGCTTGAGACGAGGGTAGAGAGTGCGACGCCGCGCAGGAGAGGGGAGAGG') ]
    result = match_primers(primers=coordinates, alts=alts, bamReads=mapperReadsOverlappingPosition, thresh=20)
    expected = PrimerCount(24, 'G', altTotal=1, alt_Primers_overlap=1, altPrimersWithinEnd=1, primer_base='A')                                                                                                                               
  # TODO FIX SO IT FAILS again
   # assert True

  # eq_(result, expected)
    assert (result == expected)
    pass
