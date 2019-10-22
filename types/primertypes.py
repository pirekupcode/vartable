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
