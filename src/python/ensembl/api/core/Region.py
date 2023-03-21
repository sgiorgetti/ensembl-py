from ensembl.api.core.Metadata import Metadata
from ensembl.api.core.Assembly import CoordSystem
from enum import Enum
from typing import Union, Optional

__all__ = [ 'Region', 'RegionCode', 'RegionTopology' ]

class RegionCode(Enum):
    CHROMOSOME = 1
    PLASMID = 2
    SCAFFOLD = 3
    CONTIG = 4

class RegionTopology(Enum):
    LINEAR = 1
    CIRCULAR = 2

class Region():
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
                 name: str, # e.g. '22'
                 topology: Union[str, RegionTopology],
                 length: int,
                 coord_system: Union[CoordSystem, str],
                 fasta_file: Optional[str] = None
                ) -> None:
        self.__name = name
        self._topology = topology if isinstance(topology, RegionTopology) else RegionTopology[topology.upper()]
        self._length = length
        self._coord_system = coord_system
        self._fasta_file_url = fasta_file

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.name},{self._topology.name},{self._code})'

    @property
    def coord_system(self) -> Union[CoordSystem, str]:
        return self._coord_system

    @property
    def code(self) -> Union[RegionCode, str]:
        return self._coord_system.name
    
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
        return self._coord_system.is_toplevel
    
    @property
    def rank(self) -> int:
        return self._coord_system.rank

    @property
    def topology(self) -> RegionTopology:
        return self._topology
    
    def get_metadata(self) -> dict[str, Metadata]:
        return self._metadata
    
    def set_metadata(self, metadata_item: dict[str, Metadata]) -> None:
        self._metadata = metadata_item

