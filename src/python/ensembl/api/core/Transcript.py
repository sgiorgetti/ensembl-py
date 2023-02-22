from ensembl.api.core.Metadata import Metadata
from ensembl.api.core.Exon import Exon
from typing import Dict, List, Any

class Transcript(object):
    """Representation of a transcript.

    TO DO!!!

    The location of a transcript is always provided 5' -> 3' on the forward strand.
    Thus the start should always be lower than the end.

    Attributes:
      exons: List of exons constituting the transcript.
      start: Start position of the transcript.
      end: End position of the transcript.
      strand: Strand of the transcript, either + or -.
      location_name: Name of the region the transcript is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the transcript.
      external_id: Name of the transcript.
      introns: List of introns constituting the transcript, should be 0 if exons size is one.
      cds_genomic_start: Genomic start position of the CDS, always be lower than cds_genomic_end.
      cds_genomic_end: Genomic end position of the CDS, always be greater than cds_genomic_start.
      cds_sequence: Translateable cDNA sequence.
      translation_sequence: Protein sequence translated from cds_sequence.
    """

    __type = 'Transcript'
 
    def __init__(self, name: str=None) -> None:
        self._name = name

    def __repr__(self) -> str:
        return f"I am an Ensembl transcript! I am {self._name}!!!"

    @property
    def external_id(self) -> str:
        return self._external_id

    @external_id.setter
    def external_id(self, value: str) -> None:
        self._external_id = value

    @property
    def name(self) -> str:
        return self._name
    
    @property
    def stable_id(self) -> str:
        return f"{self._unversioned_stable_id}.{self._version}"

    @stable_id.setter
    def stable_id(self, value: str) -> None:
        (self._unversioned_stable_id, self._version) = value.split('.')
    
    @property
    def symbol(self) -> str:
        return self._symbol

    @symbol.setter
    def symbol(self, value: str) -> None:
        self._symbol = value
    
    @property
    def type(self) -> str:
        return self.__type

    @property
    def unversioned_stable_id(self) -> str:
        return self._unversioned_stable_id

    @unversioned_stable_id.setter
    def unversioned_stable_id(self, value) -> None:
        self._unversioned_stable_id = value

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, value) -> None:
        self._version = value

    def get_exons(self) -> List[Exon]:
        return self._exons
    
    def set_exons(self, exons: List[Exon]) -> None:
        self._exons = exons

    def get_introns(self) -> List[Any]:
        return self._introns
    
    def set_introns(self, introns: List[Any]) -> None:
        self._introns = introns
    
    def get_metadata(self) -> Dict[str, Metadata]:
        return self._metadata
    
    def set_metadata(self, metadata_item: Dict[str, Metadata]) -> None:
        self._metadata = metadata_item

    def get_slice(self) -> Any:
        return self._slice
    
    def set_slice(self, slice: Any) -> None:
        self._slice = slice
