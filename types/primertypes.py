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
    querySequenceSoftClipped: str
        
    def query_base_at_reference_position(position: int) -> str: ....
  
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
        
        
def match_primers(primers: List[MappingCoordinate], alts: List[Alt], bamReads: Iterator[BamRead]) -> List[PrimerCount]: ...
