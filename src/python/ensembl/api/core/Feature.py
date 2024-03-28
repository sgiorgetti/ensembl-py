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
"""Generic Feature module"""

__all__ = [ "Feature", "EnsemblFeature" ]

from typing import Union
import warnings
from . import Analysis, Location, RegionType, Sequence, Slice, Strand

class Feature():
    """
    Ensembl specific sequence feature.
    This is the Base feature class from which all Ensembl features inherit.
    It provides a bare minimum functionality that all features require.  It
    basically describes a location on a sequence in an arbitrary coordinate
    system.
    """

    __type = 'feature'

    def __init__(self,
            location: Location,
            strand: Strand,
            reg_slice: Slice = None,
            analysis: Analysis = None,
            internal_id: int = None,
            sequence: Sequence = None
        ) -> None:

        if reg_slice and not isinstance(reg_slice, Slice):
            raise ValueError(f"Slice argument must be of type {type(Slice)}")
        # if analysis and not isinstance(analysis, Slice):
        #     raise ValueError('Analysis argument must be a Bio::EnsEMBL::Analysis')
        if not location or not isinstance(location, Location):
            raise ValueError(f"Location argument must be defined and of type {type(Location)}")
        if not strand or not isinstance(strand, Strand):
            raise ValueError(f"Location argument must be defined and of type {type(Strand)}")
        if location.start > location.end:
            raise ValueError(f"Invalid Location: start({location.start}) > end({location.end})")
        try:
            if not reg_slice.is_circular and location not in reg_slice.location:
                raise ValueError(f"Feature {internal_id} location not contained in given Slice")
        except TypeError:
            pass
        if not internal_id:
            internal_id = f"{reg_slice.region.name}:{location.start}:{location.end}:{strand.value}"
        self._location = location
        self._slice = reg_slice
        self._strand = strand
        self._analysis = analysis
        self._internal_id = internal_id
        self._sequence = sequence
        self._attributes = {}

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._internal_id})'

    @property
    def feature_type(self) -> str:
        return self.__type

    @property
    def location(self) -> Location:
        return self._location

    @property
    def start(self) -> int:
        return self._location.start

    @start.setter
    def start(self, value: int) -> None:
        if not isinstance(value, int):
            raise ValueError('start argument must be int')
        self._location.start = value

    @property
    def end(self) -> int:
        return self._location.end

    @end.setter
    def end(self, value: int) -> None:
        if not isinstance(value, int):
            raise ValueError('end argument must be int')
        self._location.end = value

    @property
    def strand(self) -> Strand:
        return self._strand

    @strand.setter
    def strand(self, value: Union[Strand, int]) -> None:
        if not value:
            warnings.warn("Provided value is empty or None", UserWarning)
            return
        if not isinstance(value, Strand) or not isinstance(value, int):
            raise ValueError('strand argument must be Strand or int')
        self._strand = value if isinstance(value, Strand) else Strand(value)

    @property
    def analysis(self) -> Analysis:
        return self._analysis

    @analysis.setter
    def analysis(self, value) -> None:
        self._analysis = value

    @property
    def internal_id(self) -> int:
        return self._internal_id

    @internal_id.setter
    def internal_id(self, value: int) -> None:
        self._internal_id = value

    @property
    def length(self) -> int:
        raw_len = self._location.end - self._location.start
        if raw_len >= 0:
            return raw_len + 1
        if self._slice.is_circular():
            # if circular, we can work out the length of an origin-spanning
            # feature using the size of the underlying region.
            return  raw_len%self._slice.length + 1
        raise ValueError('Cannot determine length of non-circular feature where start > end')

    @property
    def seq(self) -> str:
        raise NotImplementedError()

    @seq.setter
    def seq(self, value) -> None:
        raise NotImplementedError()

    def set_attribs(self, attribs: dict[str, str]) -> None:
        self._attributes = attribs

    def get_attribs(self) -> dict[str, str]:
        return self._attributes

    def add_attrib(self, code: str, value: str):
        self._attributes[code] = value

    def get_attrib(self, code: str) -> str:
        return self._attributes.get(code)

    def get_slice(self) -> Slice:
        return self._slice

    def set_slice(self, gen_slice: Slice) -> None:
        self._slice = gen_slice

    def feature_slice(self) -> Slice:
        return Slice(self._slice.region, self._location)

    @property
    def seq_region_name(self) -> str:
        return self._slice.region.name


    #  FROM SLICE - TO BE CONFIRMED!!
    @property
    def seq_region_length(self) -> int:
        return self._slice.length

    @property
    def seq_region_type(self) -> RegionType:
        return self._slice.region_type

    @property
    def seq_region_cs(self) -> str:
        return self._slice.location.coordinate_system.name

class EnsemblFeature(Feature):

    __type = 'ensembl_feature'

    def __init__(self,
                 location: Location,
                 strand: Strand,
                 reg_slice: Slice,
                 analysis: Analysis = None,
                 internal_id: Union[str, int] = None,
                 sequence: Sequence = None,
                 biotype: str = None,
                 source: str = 'ensembl',
                 stable_id: str = None,
                 version: int = 1,
                 is_current: bool = True) -> None:
        self._stable_id = stable_id
        self._version = version
        self._biotype = biotype
        self._source = source
        self._is_current = is_current
        self._attributes = {}

        super().__init__(location, strand, reg_slice, analysis, internal_id, sequence)


    def __repr__(self) -> str:
        if self._stable_id:
            return f'{self.__class__.__name__}({self._stable_id})'
        return super().__repr__()

    @property
    def stable_id_version(self) -> str:
        return f"{self._stable_id}.{self._version}"

    @property
    def stable_id(self) -> str:
        return self._stable_id

    @stable_id.setter
    def stable_id(self, value: str) -> None:
        self._stable_id = value

    @property
    def feature_type(self) -> str:
        return self.__type

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, value) -> None:
        self._version = value

    @property
    def internal_id(self) -> int:
        return self._internal_id

    @property
    def biotype(self) -> str:
        return self._biotype

    @biotype.setter
    def biotype(self, value: str) -> None:
        self._biotype = value

    @property
    def source(self) -> str:
        return self._source

    @source.setter
    def source(self, value: str) -> None:
        self._source = value

    @property
    def is_current(self) -> bool:
        return self._is_current

    @is_current.setter
    def is_current(self, value: bool) -> None:
        self._is_current = value
