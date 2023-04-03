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

from .Strand import Strand
from .Slice import Slice

from functools import cache
from typing import Union
import warnings

__all__ = [ 'Feature' ]

class Feature():
    """
    Ensembl specific sequence feature.
    This is the Base feature class from which all Ensembl features inherit.
    It provides a bare minimum functionality that all features require.  It
    basically describes a location on a sequence in an arbitrary coordinate
    system.
    """
    def __init__(self,
            start: int = 1,
            end: int = None,
            strand: Strand = Strand.FORWARD,
            slice: Slice = None,
            analysis: str = None,
            internal_id: int = None,
            created_date: str = None,
            modified_date: str = None,
            seqname = None) -> None:
                 
        if slice and not isinstance(slice, Slice):
            raise ValueError('Slice argument must be a Bio::EnsEMBL::Slice')
        # if analysis and not isinstance(analysis, Slice):
        #     raise ValueError('Analysis argument must be a Bio::EnsEMBL::Analysis')
        if start and end:
            if not isinstance(start, int) or not isinstance(end, int):
                raise ValueError('Start/End arguments must be int')
            if end+1 < start and slice and not slice.is_circular():
                raise ValueError(f'Start {start} must be less than or equal to end+1 {end+1}')
        
        self._start = start
        self._end = end
        self._strand = strand
        self._slice = slice
        self._analysis = analysis
        self._internal_id = internal_id
        self._seqname = seqname
        self._created_date = created_date
        self._modified_date = modified_date

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(internal_id{self._internal_id}){self._slice.name}:{self._start}:{self._end}'
    
    def __eq__(self, __o: object) -> bool:
        if not isinstance(__o, Feature):
            return False
        
        # If the features has the same dbID, they are equal.
        if self.internal_id and __o.internal_id:
            return True if self.internal_id == __o.internal_id else False
        
        # If the features have the same start, end, strand and seq_region_id,
        # and analysis_id, they are equal.
        if (self.start == __o.start 
            and self.end == __o.end
            and self.strand == __o.strand
            and self._slice.seq_region_id == __o._slice.seq_region_id
            and self.analysis == __o.analysis):
            return True
        else:
            return False
    
    @property
    def start(self) -> int:
        return self._start
    
    @start.setter
    def start(self, value: int) -> None:
        if not isinstance(value, int):
            raise('start argument must be int')
        self._start = value

    @property
    def end(self) -> int:
        return self._end
    
    @end.setter
    def end(self, value: int) -> None:
        if not isinstance(value, int):
            raise('end argument must be int')
        self._end = value

    @property
    def strand(self) -> Strand:
        return self._strand
    
    @strand.setter
    def strand(self, value: Union[Strand, int]) -> None:
        if not value:
            warnings.warn(f"Provided value is empty or None", UserWarning)
            return
        if not isinstance(value, Strand) or not isinstance(value, int):
            raise('strand argument must be Strand or int')
        self._strand = value if isinstance(value, Strand) else Strand(value)

    @property
    def analysis(self):
        return self._analysis
    
    @analysis.setter
    def analysis(self, value) -> None:
        self._analysis = value

    @property
    def internal_id(self) -> int:
        return self._internal_id
    
    @property
    def dbID(self) -> int:
        return self._internal_id
    
    @internal_id.setter
    def internal_id(self, value: int) -> None:
        self._internal_id = value

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

    @property
    def length(self) -> int:
        if self._end < self._start:
            # if circular, we can work out the length of an origin-spanning
            # feature using the size of the underlying region.
            if self._slice.is_circular():
                len = self._slice.seq_region_length - (self._start-self._end) + 1
                return len
            raise Exception(f'Cannot determine length of non-circular feature where start > end')
        return self._end - self._start + 1
    
    @property
    @cache
    def seqname(self, value: str) -> str:
        if value:
            self._seqname = value
        if not self._seqname and self._slice:
            self._seqname = self._slice.name
        return self._seqname

    
    def get_slice(self) -> Slice:
        return self._slice
    
    def set_slice(self, slice: Slice) -> None:
        self._slice = slice


    def transform(self, cs_name, cs_version, to_slice):
        raise NotImplementedError()
    
    def transfer(self, to_slice):
        raise NotImplementedError()
    
    def project_to_slice(self, to_slice):
        raise NotImplementedError()
    
    def project(self, cs_name, cs_version):
        raise NotImplementedError()
    
    def feature_slice(self) -> Slice:
        if not self._slice:
            warnings.warn(f"Cannot obtain Feature_Slice for feature without attached slice", UserWarning)
            return None
        
        s = Slice(
            # aaaa
        )

    @property
    def seq_region_name(self) -> str:
        return self._slice.seq_region_name() if self._slice else ''
    
    #  FROM SLICE
    # sub seq_region_length
    @property
    def seq_region_length(self) -> int:
        return self._slice.seq_region_length
    # sub seq_region_strand 
    @property
    def seq_region_strand(self) -> Strand:
        st = self._slice.strand.value * self._strand.value
        return Strand(st) if self._slice.strand else Strand.UNDEFINED
    # sub seq_region_start
    @property
    def seq_region_start(self) -> int:
        if self._slice.strand == Strand.FORWARD:
            return self._slice.seq_region_start + self._start - 1
        else:
            return self._slice.seq_region_end - self._end + 1
    # sub seq_region_end
    @property
    def seq_region_end(self) -> int:
        if self._slice.strand == Strand.FORWARD:
            return self._slice.seq_region_start + self._end - 1
        else:
            return self._slice.seq_region_end + self._start + 1
    # sub coord_system_name
    @property
    def coord_system_name(self) -> str:
        return self._slice.coord_system_name
    
    @property
    def coord_system_version(self) -> str:
        return self._slice.coord_system_version

    # sub seq
    def seq(self):
        raise NotImplementedError()

    def get_summary(self) -> dict:
        summary = {}
        summary['id'] = self.internal_id
        summary['version'] = self._version if self._version else ''
        summary['start'] = self.seq_region_start()
        summary['end'] = self.seq_region_end()
        summary['strand'] = self._strand
        summary['seq_region_name'] = self.seq_region_name()
        return summary

