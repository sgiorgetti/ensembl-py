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

from . import SplicedExon
from .Strand import Strand
from typing import Union


__all__ = [ 'Translation' ]


class Translation():
    """A class representing the translation of a
    transcript

    A Translation object defines the CDS and UTR regions of a Transcript
    through the use of start_Exon/end_Exon, and start/end attributes.

    SYNOPSIS
    translation = ensembl.api.core.Translation(
        start_exon = $exon1,
        end_exon   = $exon2,
        seq_start  = 98,
        seq_end    = 39
    );
    # stable ID setter
    translation.stable_id = 'ENSP00053458'
    # get start and end position in start/end exons
    start = translation.start
    end   = translation.end
    """
    __type = 'translation'

    """
    Arg [start_exon] : The Exon object in which the translation (CDS) starts
    Arg [end_exon]   : The Exon object in which the translation (CDS) ends
    Arg [seq_start]  : The offset in the start_Exon indicating the start
                        position of the CDS.
    Arg [seq_end]    : The offset in the end_Exon indicating the end
                        position of the CDS.
    Arg [stable_id]  : The stable identifier for this Translation
    Arg [version]    : The version of the stable identifier
    Arg [internal_id]: The internal identifier of this Translation
    Arg [adaptor]    : The TranslationAdaptor for this Translation
    Arg [seq]        : Manually sets the peptide sequence of this translation.
                        May be useful if this translation is not stored in
                        a database.
    Arg [created_date]: the date the translation was created
    Arg [modified_date]: the date the translation was modified
    Example    : tl = Translation(start_exon = ex1,
                                end_exon = ex2,
                                seq_start  = 98,
                                seq_end    = 39)
    Description: Constructor.  Creates a new Translation object
    Returntype : ensembl.api.core.Translation
    Exceptions : none
    Caller     : general
    Status     : At Risk
                : under development
    """
    def __init__(self,
        start_exon: SplicedExon,
        end_exon: SplicedExon,
        seq_start: int,
        seq_end: int,
        stable_id: str = None,
        version: int = None,
        internal_id: int = None,
        seq: str = None,
        created_date: str = None,
        modified_date: str = None,
        genomic_start: int = None,
        genomic_end: int = None
    ) -> None:
        self._start_exon = start_exon
        self._end_exon = end_exon
        self._seq_start = seq_start
        self._seq_end = seq_end
        self._stable_id = stable_id
        self._version = version
        self._internal_id = internal_id
        self._seq = seq
        self._created_date = created_date
        self._modified_date = modified_date
        self._genomic_start = genomic_start
        self._genomic_end = genomic_end
        self._attributes = {}

        self.__post_init__()
        
    
    def __post_init__(self):
        if self._genomic_start is None:
            if self._start_exon.strand == Strand.REVERSE:
                self._genomic_start = self._end_exon.end - (self._seq_end - 1)
            else:
                self._genomic_start = self._start_exon.start - (self._seq_start - 1)
        if self._genomic_end is None:
            if self._start_exon.strand == Strand.REVERSE:
                self._genomic_end = self._start_exon.end - (self._seq_start - 1)
            else:
                self._genomic_end = self._end_exon.start - (self._seq_end - 1)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(internal_id{self._internal_id})'

    @property
    def start_exon(self) -> SplicedExon:
        return self._start_exon
    
    @start_exon.setter
    def start_exon(self, value: SplicedExon) -> None:
        if not isinstance(value, SplicedExon):
            raise AttributeError(f'Got to have an Exon object, not a {type(value)}')
        self._start_exon = value

    @property
    def end_exon(self) -> SplicedExon:
        return self._end_exon
    
    @end_exon.setter
    def end_exon(self, value: SplicedExon) -> None:
        if not isinstance(value, SplicedExon):
            raise AttributeError(f'Got to have an Exon object, not a {type(value)}')
        self._end_exon = value

    @property
    def seq_start(self) -> int:
        return self._seq_start
    
    @property
    def seq_end(self) -> int:
        return self._seq_end
    
    @property
    def internal_id(self) -> int:
        return self._internal_id
    
    @property
    def dbID(self) -> int:
        return self._internal_id
    
    @property
    def stable_id(self) -> str:
        return self._stable_id
    
    @property
    def version(self) -> int:
        return self._version
    
    @property
    def seq(self) -> str:
        return self._seq
    
    @property
    def length(self) -> int:
        if not self.seq:
            return 0
        return len(self.seq)
    
    @property
    def created_date(self) -> str:
        return self._created_date

    @created_date.setter
    def created_date(self, value: str) -> None:
        self._created_date = value

    @property
    def modified_date(self) -> str:
        return self._modified_date

    @modified_date.setter
    def modified_date(self, value: str) -> None:
        self._modified_date = value
    
    def set_attribs(self, attribs: dict[str, str]) -> None:
        self._attributes = attribs

    def get_attribs(self, att_code: str = None) -> Union[dict, tuple]:
        if att_code is None:
            return self._attributes
        return tuple(att_code, self._attributes.get(att_code))
    
    def add_attrib(self, code: str, value: str):
        self._attributes[code] = value