from ensembl.api.core.Metadata import Metadata
from enum import Enum
from typing import Dict, Union, Optional

class RegionCode(Enum):
    CHROMOSOME = 1
    PLASMID = 2
    SCAFFOLD = 3
    CONTIG = 4

class RegionTopology(Enum):
    LINEAR = 1
    CIRCULAR = 2

class Region(object):
    """Representation of a region.
    The Region data type describes the coordinate system that contains Features.
    In biological terms, it represents a naturally occurring or an artificially produced polynucleotide.
    Attributes:
      name: The name of the Region.
      code: It identifies if the region is 'CHROMOSOME', 'PLASMID' or 'SCAFFOLD'.
      topology: 'LINEAR' or 'CIRCULAR'.
      length: Length of the region in nucleotides.
      fasta_file: URL to a FASTA file containing the DNA of the region.
      assembly: The Assembly the Region belongs to. - TO ADD
      metadata: Metadata associated to the Region.
    Raises:
      Exception: ???
    """
 
    def __init__(self, 
                 name: str,
                 code: Union[str, RegionCode],
                 is_top_level: bool,
                 rank: Optional[int],
                 topology: Union[str, RegionTopology],
                 length: int,
                 assembly: str,
                 fasta_file: Optional[str] = None
                ) -> None:
        self.__name = name
        self._code = code
        self._is_top_level = is_top_level,
        self._rank = rank,
        self._topology = topology if isinstance(topology, RegionTopology) else RegionTopology[topology.upper()]
        self._length = length
        self._assembly = assembly
        self._fasta_file_url = fasta_file

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.name},{self._topology.name},{self._code},{self._assembly})'

    @property
    def assembly(self) -> str:
        return self._assembly

    @property
    def code(self) -> RegionCode:
        return self._code
    
    @property
    def name(self) -> str:
        return self.__name

    @property
    def fasta_file_url(self) -> str:
        return self._fasta_file_url

    @property
    def length(self) -> int:
        return self._length
    
    @property
    def is_top_level(self) -> bool:
        return self._is_top_level
    
    @property
    def rank(self) -> int:
        return self._rank

    @property
    def topology(self) -> RegionTopology:
        return self._topology
    
    def get_metadata(self) -> Dict[str, Metadata]:
        return self._metadata
    
    def set_metadata(self, metadata_item: Dict[str, Metadata]) -> None:
        self._metadata = metadata_item

