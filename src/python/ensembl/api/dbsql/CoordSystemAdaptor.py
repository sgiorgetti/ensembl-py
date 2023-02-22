__all__ = ['CoordSystemAdaptor', 'CoordSystem']

from sqlalchemy import select, func
from sqlalchemy.orm.exc import NoResultFound

from ensembl.core.models import CoordSystem as CoordSystemORM
from ensembl.core.models import Meta as MetaORM

from typing import Any, List, Optional

from ensembl.database.dbconnection import DBConnection
from ensembl.api.dbsql.DBAdaptor import ArgumentError

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
                 species: Optional[str] = None
                 ) -> None:

        if top_level:
            if rank != 0:
                raise ArgumentError(f'RANK argument must be 0 if TOP_LEVEL is True')

            if name != 'toplevel':
                raise ArgumentError(f'The NAME argument must be "toplevel" if TOP_LEVEL is True')

            if sequence_level:
                raise ArgumentError(f'SEQUENCE_LEVEL argument must be False if TOP_LEVEL is True')

            default = False
        else:
            if rank == 0:
                raise ArgumentError(f'RANK argument must be non-zero unless TOP_LEVEL is True')
            if name == 'toplevel':
                raise ArgumentError(f'The NAME argument cannot be "toplevel" if TOP_LEVEL is False')

        assert isinstance(rank, int)
        if rank <= 0:
            raise ArgumentError(f'Rank must be non-negative integer number')

        self._name = name
        self._version = version
        self._rank = rank
        self._top_level = top_level
        self._sequence_level = sequence_level
        self._default = default
        self._species_prod_name = species


    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._species_prod_name}:{self._name},{self._rank}-{self._version})'

    def __eq__(self, __o: object) -> bool:
        if not isinstance(__o, CoordSystem):
            raise ArgumentError('Argument must be a CoordSystem')
        if self._version == __o.version and self._name == __o.name:
            return True
        return False

    @property
    def name(self) -> str:
        return self._name

    @property
    def version(self) -> str:
        return self._name if self._name is not None else ''

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
    def species(self) -> str:
        return self._species_prod_name


class CoordSystemAdaptor():
    """Contains all the coordinate system related functions over CoordSystem ORM
    This adaptor allows the querying of information on the coordinate
    system.
    Note that many coordinate systems do not have a concept of a version
    for the entire coordinate system (though they may have a per-sequence
    version).  The 'chromosome' coordinate system usually has a version
    (i.e. the assembly version) but the clonal coordinate system does not
    (despite having individual sequence versions).  In the case where a
    coordinate system does not have a version an empty string ('') is used
    instead.
    """
    
    @classmethod
    def fetch_all(cls, dbconnection: DBConnection) -> List[CoordSystem]:
        """
        Arg [1]    : dbconnection: DBConnection
                     The DB connection object
        Example    : for cs in CoordSystemAdaptor.fetch_all(dbconnection):
                        print(f"{cs.name} {cs.version}";
        Description: Retrieves every coordinate system defined in the DB.
                     These will be returned in ascending order of species_id and rank. I.e.
                     The coordinate system with lower rank would be first in the
                     array.
        Returntype : List[ensembl.dbsql.CoordSystem]
        Exceptions : ensembl.api.dbsql.DBAdaptor.ArgumentError
        Caller     : general
        Status     : Alpha
        """
        if not dbconnection:
            raise ArgumentError()

        stmt = (select(CoordSystemORM, MetaORM.meta_value.label('species_prod_name'))
        .join(MetaORM, CoordSystemORM.species_id == MetaORM.species_id)
        .where(MetaORM.meta_key == 'species.production_name')
        .order_by(CoordSystemORM.species_id, CoordSystemORM.rank)
        )

        cs_list = []
        with dbconnection.connect() as conn:
            rows = conn.execute(stmt)
            for cnt, row in enumerate(rows):
                toplevel = True if row.name == 'top_level' else False
                seqlevel = True if 'sequence_level' in row.attrib else False
                default = True if 'default' in row.attrib else False

                cs = CoordSystem(row.name, row.version, row.rank, toplevel, seqlevel, default, row.species_prod_name)
                cs_list.append(cs)
            if cnt <= 0:
                raise NoResultFound(f'Could not find any coordinate system')
            return cs_list
    
    @classmethod
    def fetch_by_rank(cls, dbconnection: DBConnection, rank: int) -> List[CoordSystem]:
        """
        Arg [1]    : int rank
        Example    : cs_list = CoordSystemAdaptor.fetch_by_rank(1)
        Description: Retrieves a CoordinateSystem via its rank. 0 is a special
                     rank reserved for the pseudo coordinate system 'toplevel'.
                     undef is returned if no coordinate system of the specified rank
                     exists.
        Returntype : List[ensembl.dbsql.CoordSystem]
        Exceptions : ensembl.api.dbsql.DBAdaptor.ArgumentError
        Caller     : general
        Status     : Alpha
        """
        if rank < 0:
            raise ArgumentError('Rank argument must be a non-negative integer.')
        
        if rank == 0:
            pass # return fetch_top_level

        stmt = (select(CoordSystemORM, MetaORM.meta_value.label('species_prod_name'))
        .join(MetaORM, CoordSystemORM.species_id == MetaORM.species_id)
        .where(MetaORM.meta_key == 'species.production_name')
        .where(CoordSystemORM.rank == rank)
        .order_by(CoordSystemORM.species_id, CoordSystemORM.rank)
        )

        cs_list = []
        with dbconnection.connect() as conn:
            rows = conn.execute(stmt)
            for cnt, row in enumerate(rows):
                toplevel = True if row.name == 'top_level' else False
                seqlevel = True if 'sequence_level' in row.attrib else False
                default = True if 'default' in row.attrib else False

                cs = CoordSystem(row.name, row.version, row.rank, toplevel, seqlevel, default, row.species_prod_name)
                cs_list.append(cs)
            if cnt <= 0:
                raise NoResultFound(f'Could not find any coordinate system')
            return cs_list

    @classmethod
    def fetch_by_name(cls, dbconnection: DBConnection, name: str, version: Optional[str] = None) -> List[CoordSystem]:
        """
        Arg [1]    : str name
                     The name of the coordinate system to retrieve.  Alternatively
                     this may be an alias for a real coordinate system.  Valid
                     aliases are 'toplevel' and 'seqlevel'.
        Arg [2]    : str version
                     The version of the coordinate system to retrieve.  If not
                     specified the default version will be used.
        Example    : cs_list = CoordSystemAdaptor.fetch_by_name('contig')
                     cs_list = CoordSystemAdaptor.fetch_by_name('chromosome','GRCh37')
        Description: Retrieves a coordinate system by its name
        Returntype : List[ensembl.dbsql.CoordSystem]
        Exceptions : ensembl.api.dbsql.DBAdaptor.ArgumentError
        Caller     : general
        Status     : Alpha
        """
        
        if not version:
            stmt = (select(CoordSystemORM, MetaORM.meta_value.label('species_prod_name'))
            .join(MetaORM, CoordSystemORM.species_id == MetaORM.species_id)
            .where(MetaORM.meta_key == 'species.production_name')
            .where(func.lower(CoordSystemORM.name) == name.lower())
            .order_by(CoordSystemORM.species_id, CoordSystemORM.rank)
        )
        else:
            stmt = (select(CoordSystemORM, MetaORM.meta_value.label('species_prod_name'))
            .join(MetaORM, CoordSystemORM.species_id == MetaORM.species_id)
            .where(MetaORM.meta_key == 'species.production_name')
            .where(func.lower(CoordSystemORM.name) == name.lower())
            .where(func.lower(CoordSystemORM.version) == version.lower())
            .order_by(CoordSystemORM.species_id, CoordSystemORM.rank)
            )

        cs_list = []
        with dbconnection.connect() as conn:
            rows = conn.execute(stmt)
            for cnt, row in enumerate(rows):
                toplevel = True if row.name == 'top_level' else False
                seqlevel = True if 'sequence_level' in row.attrib else False
                default = True if 'default' in row.attrib else False

                cs = CoordSystem(row.name, row.version, row.rank, toplevel, seqlevel, default, row.species_prod_name)
                cs_list.append(cs)
            if cnt <= 0:
                raise NoResultFound(f'Could not find any coordinate system')
            return cs_list
