from ensembl.api.core.Metadata import Metadata
from ensembl.api.core.Slice import Slice
from ensembl.api.core.Location import Location
from typing import Dict, Optional

__all__ = ['Exon', 'SplicedExon']

class Exon(object):
    """Representation of an exon.
    The location of an exon is always provided 5' -> 3' on the forward strand.
    Thus the start should always be lower than the end.
    We use phase and end phase to provide information about the coding potential
    and where the next codon start at the start of the exon.
    Attributes:
      start: Start position of the exon.
      end: End position of the exon.
      strand: Strand of the exon, either + or -.
      location_name: Name of the region the exon is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the exon.
      public_identifier: Name of the exon.
      exon_start_phase: Phase of the exon, either -1, 0, 1 or 2.
      exon_end_phase: End phase of the exon, either -1, 0, 1 or 2.
    Raises:
      Exception: when start is greater than end
    """

    __type = 'Exon'
 
    def __init__(self,
                 stable_id: str,
                 slice: Slice,
                 phase: int,
                 end_phase: int
                ) -> None:
        if not stable_id:
            raise ValueError()
        if not Slice:
            raise ValueError()
        if stable_id.find('.') < 0:
            raise ValueError()
        (self._unversioned_stable_id, self._version) = stable_id.split('.')
        self._slice = slice
        self._phase = phase
        self._end_phase = end_phase
        self._metadata = {}
    # (1, 1, 76835377, 76835502, -1, -1, -1, 1, 0, 'ENSE00002089356', 1, datetime.datetime(2022, 7, 5, 10, 44, 45), datetime.datetime(2022, 7, 5, 10, 44, 45))

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.stable_id})'
    
    @property
    def stable_id(self) -> str:
        return f"{self._unversioned_stable_id}.{self._version}"

    @stable_id.setter
    def stable_id(self, value: str) -> None:
        (self._unversioned_stable_id, self._version) = value.split('.')
    
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

    @property
    def phase(self) -> int:
        return self._phase

    @property
    def end_phase(self) -> int:
        return self._end_phase

    @property
    def frame(self):
        location = self._slice.location
        if self._phase == -1:
            return '.' # gff convention for no frame info
        if self._phase == 0:
            return location.start%3
        if self._phase == 1:
            return (location.start + 2)%3
        if self._phase == 2:
            return (location.start + 1)%3
        raise Exception(f"bad phase in exon {self.phase}")
    
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


class SplicedExon(Exon):
    def __init__(self,
                 stable_id: str,
                 slice: Slice,
                 phase: int,
                 end_phase: int,
                 index: int,
                 relative_location: Location = None,
                ) -> None:
        Exon.__init__(self,
                       stable_id,
                       slice,
                       phase,
                       end_phase
        )
        self._index = index
        self._relative_location = relative_location
    
    @property
    def index(self) -> int:
        return self._index
    
    def get_relative_location(self) -> Location:
        return self._relative_location
    
    def set_relative_location(self, relative_location: Location) -> None:
        self._relative_location = relative_location