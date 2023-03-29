# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from .Slice import Slice
from .Strand import Strand
from .Feature import Feature
from typing import Union

__all__ = [ 'Exon', 'SplicedExon' ]

class Exon(Feature):
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

    __type = 'exon'
 
    def __init__(self,
                 stable_id: str,
                 version: int,
                 phase: int,
                 end_phase: int,
                 internal_id: Union[str, int] = None,
                 slice: Slice = None,
                 start: int = None,
                 end: int = None,
                 strand: Strand = None,
                 analysis: str = None,
                 is_constitutive: bool = False,
                 is_current: bool = True,
                 created_date = None,
                 modified_date = None
                ) -> None:
        if not stable_id:
            raise ValueError()
        if not Slice:
            raise ValueError()
        if stable_id.find('.') > 0:
            raise ValueError()
        if not isinstance(phase, int):
            raise ValueError('phase argument must be int')
        if phase not in (-1, 0, 1, 2):
            raise ValueError(f'Bad value {phase} for exon phase: it must be any of (-1, 0, 1, 2)')
        if end_phase is None:
            raise ValueError(f"No end phase set in Exon. You must set it explicitly.")
        self._stable_id = stable_id
        self._version = version
        self._slice = slice
        self._phase = phase
        self._end_phase = end_phase
        self._is_constitutive = is_constitutive
        self._is_current = is_current
        self._internal_id = internal_id

        super().__init__(start, end, strand, slice, analysis, internal_id, created_date, modified_date)
    # (1, 1, 76835377, 76835502, -1, -1, -1, 1, 0, 'ENSE00002089356', 1, datetime.datetime(2022, 7, 5, 10, 44, 45), datetime.datetime(2022, 7, 5, 10, 44, 45))

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.stable_id})'
    
    @property
    def stable_id_version(self) -> str:
        return f"{self._stable_id}.{self._version}"
    
    @property
    def stable_id(self) -> str:
        return self._stable_id

    @stable_id.setter
    def stable_id(self, value: str) -> None:
        if value.find('.') > 0:
            (self._stable_id, self._version) = value.split('.')
        else:
            self._stable_id = value
    
    @property
    def type(self) -> str:
        return self.__type

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
        if self._phase == -1:
            return '.' # gff convention for no frame info
        if self._phase == 0:
            return self._start%3
        if self._phase == 1:
            return (self._start + 2)%3
        if self._phase == 2:
            return (self._start + 1)%3
        raise Exception(f"bad phase in exon {self.phase}")
    
    @property
    def is_constitutive(self):
        return self._is_constitutive

    @is_constitutive.setter
    def is_constitutive(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise ValueError(f'New value must be boolean')
        self._is_constitutive = value

    @property
    def internal_id(self):
        return self._internal_id
    


    def get_summary(self) -> dict[str, str]:
        """
        Example       : exon_summary = exon.get_summary()
        Description   : Retrieves a textual summary of this Feature.
        Returns       : Dict[str, str]
        Status        : Alpha - Intended for internal use
        """
        summary = super().get_summary()
        summary['id'] = self._stable_id
        summary['exon_id'] = self._stable_id
        summary['version'] = self.version if self.version else ''
        summary['start'] = self.start
        summary['end'] = self.end
        summary['strand'] = self.strand.value
        summary['seq_region_name'] = self.seqname
        summary['constitutive'] = str(self._is_constitutive)
        summary['ensembl_phase'] = self._phase
        summary['ensembl_end_phase'] = self._end_phase
        return summary


class SplicedExon(Exon):
    def __init__(self,
                 stable_id: str,
                 version: int,
                 phase: int,
                 end_phase: int,
                 index: int,
                 internal_id: Union[str, int] = None,
                 slice: Slice = None,
                 start: int = None,
                 end: int = None,
                 strand: Strand = None,
                 analysis: str = None,
                 is_constitutive: bool = False,
                 is_current: bool = True,
                 transcript_source = None,
                 created_date = None,
                 modified_date = None
                ) -> None:
        
        self._index = index
        self._transcript_source = transcript_source
        super().__init__(stable_id,
                       version,
                       phase,
                       end_phase,
                       internal_id,
                       slice,
                       start,
                       end,
                       strand,
                       analysis,
                       is_constitutive,
                       is_current,
                       created_date,
                       modified_date
                       )
        
    def __repr__(self) -> str:
        return super().__repr__()
    
    @property
    def index(self) -> int:
        return self._index
    
    @property
    def source(self) -> str:
        return self._transcript_source
    
    @property
    def transcript_source(self) -> str:
        return self._transcript_source


    def get_summary(self) -> dict[str, str]:
        """
        Example       : exon_summary = exon.get_summary()
        Description   : Retrieves a textual summary of this Feature.
        Returns       : Dict[str, str]
        Status        : Alpha - Intended for internal use
        """
        summary = super().get_summary()
        summary['exon_index'] = self._index
        return summary
    
    
    # 1       havana  exon    65419   65433   .       +       .       Parent=transcript:ENST00000641515;Name=ENSE00003812156;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003812156;rank=1;version=1
    def gff3_qualifiers(self) -> dict[str, Union[str, tuple[str]]]:
        qualifiers = {
            'source': self.source,
            'score': ".",
            'Name': self._stable_id,
            'constitutive': str(self._is_constitutive),
            'ensembl_phase': self._phase,
            'ensembl_end_phase': self._end_phase,
            'exon_id': self._stable_id,
            'rank': self._index,
            'version': self.version
        }

        return qualifiers