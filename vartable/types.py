from typing import Union, NamedTuple, List, Dict


class CDS:
    pass

class SimpleCDS(CDS):
    '''note: uses 1-based indexing.'''
    def __init__(self, is_reverse: bool, start: int, end: int):
        self.isReverse = is_reverse
        self.start = start
        self.end = end 

class CompoundCDS(CDS):
    def __init__(self, is_reverse: bool, regions: List[SimpleCDS]):
        self.regions = regions
        self.isReverse = is_reverse

Base = str # i.e. A, C, Y, N . . .  just an alias for string now. 

class Result:
    '''note: uses 1-based indexing.'''
    def __init__(self, codon: str, refProtein: Base, altProtein: Base, isSynoynymous: bool, codonPosition: int):
        self.codon           = codon
        self.refProtein      = refProtein
        self.altProtein      = altProtein
        self.isSynonymous   = isSynoynymous
        self.codonPostion    = codonPosition

def compound_translate(ref: str, cds: CompoundCDS, alt: List[Base], pos: int) -> Result:
    '''return the effected codon, the resulting protein for the reference and alternate translation, and whether synonymous'''
    pass

def simple_translate(ref: str, cds: SimpleCDS, alts: List[Base], pos: int) -> Result:
    '''return the effected codon, the resulting protein for the reference and alternate translation, and whether synonymous'''
    pass

def translate_all(ref: str, cdss: List[CDS], variants: Dict[int, List[Base]]) -> List[Result]:
    ''' Dict[int, List[Base]] are our alternates organized by their position
        note that we may have multiple alternates at a single position.
        In this case, altProtein becomes the '''
    pass

