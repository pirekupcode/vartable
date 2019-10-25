from nose.tools import eq_, ok_
import random
from vartable.types import CDS, CompoundCDS, SimpleCDS, \
    Result, simple_translate, compound_translate, translate_all
# might need the below
# from pyfakefs.fake_filesystem_unittest import Patcher
# from unittest.mock import patch
class Alt:
    position: int
    base: str

@dataclass
class MappingCoordinate:
    referenceStart: int
    referenceEnd: int

@dataclass
class BamRead:
    referenceEnd: int
    referenceStart: int
    queryStart: int
    queryEnd: int
    querySequence: str # without soft clipping
    querySequenceSoftClipped: str
        
  #  def query_base_at_reference_position(position: int) -> str: ....
  
# drawing from https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment 

@dataclass
class PrimerCount:
    alt: Alt  # position and alternate where this count came from 
    baseATotal: int
    baseAPrimer: int
    baseCTotal: int
    baseCPrimer: int
    baseGTotal: int
    baseGPrimer: int
    baseTTotal: int
    baseTPrimer: int
    other: List[Tuple[str, int]] # e.g. N, insertion, etc
        
        
def match_primers(primers: List[MappingCoordinate], alts: List[Alt], bamReads: Iterator[BamRead]) -> List[PrimerCount]: 

#The test primer function consolidates parsed data from the reference, vcf, bamreads consensus fasta and the primer file to verify ambiguous base calls


#def test_primer_annotate(Alt: [str, int], consensus: , MappingCoordinate:[int, int], BamRead[int, int, int, int, str,str] PrimerCount) -> None:
def test_match_primers(primers: List[MappingCoordinate], alts: List[Alt], bamReads: Iterator[BamRead]) -> List[PrimerCount]:
    ok_(True, "True is true")
#Having difficulty instantiating the PrimerCount (expected instantiating is incomplete)
    actual, expected = test_match_primers( primer: [10, 40], alt: [24, G], bamReads:[referenceEnd: 234,referenceStart:5, queryStart:9,queryEnd:154, querySequence:'ACGACGGCGCACACGGCAGCAGCAGCGACGCGGGGTTTATATATGGCACGAGCGCGCAGCGACGCAGCGCGCTTGAGACGAGGGTAGAGAGTGCGACGCCGCGCAGGAGAGGGGAGAGG' querySequenceSoftClipped:'?')
    eq_(3, 5, "Actual {0} is not Expected {1}".format(actual, expected))
    result = primer_annotate(vcf, primerfile, consensus, bamfile, options)

    # duplicate dirs are impossible anyway, but still using Counter
    ok_(True, "True is true")


