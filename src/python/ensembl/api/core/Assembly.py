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

from typing import Optional, Union

__all__ = [ 'Assembly', 'CoordSystem' ]

class Assembly():
    def __init__(self,
                 id: str,
                 name: str,
                 species: str,
                 accession_id: str = None,
                 accessioning_body: str = None,
                 gb_last_geneset_update: str = None,
                 assembly_date: str = None,
                 is_default: bool = False,
                 tolid: str = None) -> None:
        self._id = id
        self._name = name
        self._species = species
        if is_default:
            if not (accession_id and accessioning_body and gb_last_geneset_update and assembly_date):
                raise ValueError(f'accession_id, accessioning_body, gb_last_geneset_update, assembly_date must be specified for the default assembly')
        self._accession_id = accession_id
        self._accessioning_body = accessioning_body
        self._gb_last_geneset_update = gb_last_geneset_update
        self._assembly_date = assembly_date
        self._default = is_default
        self._tolid = tolid

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._species}-{self._id}-{self._accession_id})'

    @property
    def id(self) -> str:
        return self._id
    
    @property
    def name(self) -> str:
        return self._name
    
    @property
    def accession_id(self) -> str:
        return self._accession_id
    
    @property
    def accessioning_body(self) -> str:
        return self._accessioning_body
    
    @property
    def species(self) -> str:
        return self._species
    
    @property
    def gb_last_geneset_update(self) -> str:
        return self._gb_last_geneset_update
    
    @property
    def assembly_date(self) -> str:
        return self._assembly_date
    
    @property
    def is_default(self) -> bool:
        return self._default
    
    @property
    def tolid(self) -> bool:
        return self._tolid
    


class CoordSystem():
    """
    This is a simple object which contains a few coordinate system attributes:
    name, internal identifier, version.  A coordinate system is uniquely defined
    by its name and version.  A version of a coordinate system applies to all
    sequences within a coordinate system.  This should not be confused with
    individual sequence versions.

    Take for example the Human assembly.  The version 'NCBI33' applies to
    to all chromosomes in the NCBI33 assembly (that is the entire 'chromosome'
    coordinate system).  The 'clone' coordinate system in the same database would
    have no version however.  Although the clone sequences have their own sequence
    versions, there is no version which applies to the entire set of clones.

    Coordinate system objects are immutable. Their name and version, and other
    attributes may not be altered after they are created.
    """
    def __init__(self,
                 name: str,
                 version: Optional[str] = None,
                 rank: int = 0,
                 top_level: bool = False,
                 sequence_level: bool = False,
                 default: bool = True,
                 species_id: Union[str, int] = None,
                 internal_id: int = None
                 ) -> None:

        if top_level:
            if rank != 0:
                raise ValueError(f'RANK argument must be 0 if TOP_LEVEL is True')

            if name != 'toplevel':
                raise ValueError(f'The NAME argument must be "toplevel" if TOP_LEVEL is True')

            if sequence_level:
                raise ValueError(f'SEQUENCE_LEVEL argument must be False if TOP_LEVEL is True')

            default = False
        else:
            if rank == 0:
                raise ValueError(f'RANK argument must be non-zero unless TOP_LEVEL is True')
            if name == 'toplevel':
                raise ValueError(f'The NAME argument cannot be "toplevel" if TOP_LEVEL is False')

        assert isinstance(rank, int)
        if rank < 0:
            raise ValueError(f'Rank must be non-negative integer number')

        self._name = name
        self._version = version
        self._rank = rank
        self._top_level = top_level
        self._sequence_level = sequence_level
        self._default = default
        self._species_id = species_id
        self._internal_id = internal_id


    @classmethod
    def new_fast(cls, data: dict) -> None:
        cls.__init__(
            data.get('name'),
            data.get('rank'),
            data.get('top_level'),
            data.get('sequence_level'),
            data.get('default'),
            data.get('species_id'),
            data.get('internal_id')
        )

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._name}:{self._version} - rank {self._rank})'

    def __eq__(self, __o: object) -> bool:
        if not isinstance(__o, CoordSystem):
            raise ValueError('Argument must be a CoordSystem')
        if self._version == __o.version and self._name == __o.name:
            return True
        return False

    @property
    def name(self) -> str:
        return self._name

    @property
    def version(self) -> str:
        return self._version if self._version is not None else ''

    @property
    def is_toplevel(self) -> bool:
        return self._top_level

    @property
    def is_sequencelevel(self) -> bool:
        return self._sequence_level

    @property
    def is_default(self) -> bool:
        return self._default

    @property
    def rank(self) -> int:
        return self._rank

    @property
    def species_id(self) -> str:
        return self._species_id
    
    @property
    def internal_id(self) -> int:
        return self._internal_id