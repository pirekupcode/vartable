from dataclasses import dataclass
from typing import List, Iterator

@dataclass
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
    # querySequenceSoftClipped: str

@dataclass
class PrimerCount:
    position: int
    alt: Alt  # position and alternate where this count came from. the alt call is obtained from the VCF file.
    altTotal: int
    alt_Primers_overlap : int
    altPrimersWithinEnd: int
    primer_base:str
#  altPrimersWithinEnd: means it overlaps with the primer and its is within the treshold of the end
#  altPrimers: All the alternate primers that overlap witht the primer file



def match_primers(primers: List[MappingCoordinate], alts: List[Alt], bamReads: Iterator[BamRead], thresh: int) -> List[PrimerCount]:
    pass
# This is a function to calculate the query position on with reference and query coordinates
def calc_primer_position( alt_position :int, ref_start: int, query_start: int ):

  add_2 = alt_position - ref_start
  query_position = add_2 + query_start

  return query_position
