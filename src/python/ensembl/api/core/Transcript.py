from ensembl.api.core.Metadata import Metadata
from ensembl.api.core.Exon import SplicedExon
from ensembl.api.core.Slice import Slice
from ensembl.api.core.Location import Location
from typing import Dict, List, Any

__all__ = ['Transcript']

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
 
    def __init__(self,
                 stable_id: str,
                 slice: Slice,
                 relative_location: Location = None,
                ) -> None:
        if not stable_id:
            raise ValueError()
        if not Slice:
            raise ValueError()
        if stable_id.find('.') < 0:
            raise ValueError()
        (self._unversioned_stable_id, self._version) = stable_id.split('.')
        self._slice = slice
        self._relative_location = relative_location
        self._metadata = {}

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.stable_id})'
    
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

    def get_exons(self) -> List[SplicedExon]:
        return self._exons
    
    def set_exons(self, exons: List[SplicedExon]) -> None:
        self._exons = exons

    def get_introns(self) -> List[Any]:
        return self._introns
    
    def set_introns(self, introns: List[Any]) -> None:
        self._introns = introns
    
    def get_metadata(self) -> Dict[str, Metadata]:
        return self._metadata
    
    def set_metadata(self, metadata_item: Dict[str, Metadata]) -> None:
        self._metadata = metadata_item

    def add_metadata(self, meta_key: str, meta_value: Metadata) -> None:
        self._metadata[meta_key] = meta_value

    def get_slice(self) -> Slice:
        return self._slice
    
    def set_slice(self, slice: Slice) -> None:
        self._slice = slice

    def get_relative_location(self) -> Location:
        return self._relative_location
    
    def set_relative_location(self, relative_location: Location) -> None:
        self._relative_location = relative_location
