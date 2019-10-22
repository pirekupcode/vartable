@dataclass
class Alt:
    position: int
    base: str

@dataclass
class PrimerCoordinate:
    referenceStart: int
    referenceEnd: int


# see https://pysam.readthedocs.io/en/latest/api.html#pysam.PileupColumn  for below
@dataclass
class PileupColumn:
    n: int
    nsegments: int 
    pileups: List[PileupRead] 
    pos: int 
    reference_id: int 
    reference_name: str 
    reference_pos: int # zero-based 
    tid: int

    def get_mapping_qualities(self) -> List[int]: ...

    def get_num_aligned(self) -> int: ...

    def get_query_names(self) -> List[str]: ...

    def get_query_positions(self) -> List[int]: ...

    def get_query_qualities(self) -> List[int]: ...

    def get_query_sequences(self, mark_matches: bool = False, mark_ends: bool = False, add_indels: bool = False) -> List[str]: ... #TODO check return

    def __len__(self) -> int: ...

