from Metadata import Metadata
from enum import Enum
from typing import Dict

class RegionCode(Enum):
    CHROMOSOME = 1
    PLASMID = 2
    SCAFFOLD = 3

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
 
    def __init__(self, name: str, code: RegionCode, topology: RegionTopology, length: int, fasta_file: str) -> None:
        self.__name = name
        self._code = code
        self._topology = topology
        self._length = length
        self._fasta_file_url = fasta_file

    def __repr__(self) -> str:
        return f"I am {self.__name}: an Ensembl region!"

    @property
    def name(self) -> str:
        return self.__name

    @property
    def code(self) -> RegionCode:
        return self._code

    @property
    def fasta_file_url(self) -> str:
        return self._fasta_file_url

    @property
    def length(self) -> int:
        return self._length

    @property
    def topology(self) -> RegionTopology:
        return self._topology
    
    def get_metadata(self) -> Dict[str, Metadata]:
        return self._metadata
    
    def set_metadata(self, metadata_item: Dict[str, Metadata]) -> None:
        self._metadata = metadata_item

