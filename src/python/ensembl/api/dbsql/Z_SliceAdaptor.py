__all__ = ['SliceAdaptor']

from sqlalchemy import select
from sqlalchemy.schema import Table
from sqlalchemy.orm.exc import NoResultFound

from typing import Any, Union, List, Dict, Optional

from ensembl.database.dbconnection import DBConnection
from DBAdaptor import DBAdaptor, ArgumentError
from ensembl.api.core import Slice, Region, Location, Strand

class SliceAdaptor(DBAdaptor):
    """Contains all the slice related functions over Slice ORM

    Attributes:
        _dbc: DBConnection
        _slice_t: Table object
        _seq_region_t
        _coord_sys_t
        _slice: Slice object
    """
    def __init__(self, dbconnection:DBConnection=None) -> None:
        if not dbconnection:
            Exception()
        self._dbc = dbconnection
        self._slice_t = Table("seq_region", dbconnection._metadata, autoload_with=dbconnection._engine)
        self._seq_region_t = Table("seq_region", dbconnection._metadata, autoload_with=dbconnection._engine)
        self._seq_region_attrib_t = Table("seq_region_attrib", dbconnection._metadata, autoload_with=dbconnection._engine)
        self._coord_sys_t = Table("coord_system", dbconnection._metadata, autoload_with=dbconnection._engine)

    def __repr__(self) -> str:
        """Returns a string representation of this object."""
        return f'{self.__class__.__name__}(Table({self._slice})@{self._dbc._engine})'

    def get_dbconnection(self):
        return super().get_dbconnection()

    def fetch_by_seq_region(self,
                            coord_system: str,
                            seq_region_name: str,
                            start: Optional[int],
                            end: Optional[int],
                            strand: Optional[int],
                            version: Optional[str]
                            ) -> Any:
        """
        Arg [1]    : str coord_system
            the coordinate system name of the seq_region to retrieve
        Arg [2]    : str seq_region_name
            the seq_region name as specified in the core DB
        Arg [3]    : int start
            the seq_region start position
        Arg [4]    : int end
            the seq_region name as specified in the core DB
        Arg [5]    : int strand
            the seq_region name as specified in the core DB
        Arg [6]    : str version
            the coordinate system version as specified in the core DB - e.g. GRCh38
        Example    : region = region_adaptor->fetch_by_region( 'chromosome', 'X' );
        Description: Retrieves a Region from the database via its name and coordinate system
        Returntype : ensembl.api.core.Slice in native coordinates.!!!!!!!!!!!!!!!!!!!!
        Exceptions : NoResultFound
        Caller     : general
        Status     : Alpha
        """
        slices = []
        assert isinstance(start, int)
        assert isinstance(end, int)
        if not start:
            start = 1
        if not strand:
            strand = 1
        if not seq_region_name:
            raise ArgumentError('seq_region_name argument is required')
        if not coord_system:
            raise ArgumentError('coord_system argument is required')
        if not version:
            stmt = (select(
            (self._seq_region_t.c.name).label("sr_name"),
            (self._seq_region_t.c.length).label("sr_length"),
            (self._coord_sys_t.c.name).label("cs_name"),
            (self._coord_sys_t.c.version).label("cs_version"),
            (self._coord_sys_t.c.rank).label("cs_rank")
            )
            .join(self._coord_sys_t, self._seq_region_t.c.coord_system_id == self._coord_sys_t.c.coord_system_id)
            .where(self._coord_sys_t.c.name == coord_system)
            .where(self._coord_sys_t.c.attrib == 'default_version')
            .where(self._seq_region_t.c.name == seq_region_name)
            )
        else:
            stmt = (select(
            (self._seq_region_t.c.name).label("sr_name"),
            (self._seq_region_t.c.length).label("sr_length"),
            (self._coord_sys_t.c.name).label("cs_name"),
            (self._coord_sys_t.c.version).label("cs_version"),
            (self._coord_sys_t.c.rank).label("cs_rank")
            )
            .join(self._coord_sys_t, self._seq_region_t.c.coord_system_id == self._coord_sys_t.c.coord_system_id)
            .where(self._coord_sys_t.c.name == coord_system)
            .where(self._coord_sys_t.c.version == version)
            .where(self._seq_region_t.c.name == seq_region_name)
            )
        with self._dbc._engine.connect() as conn:
            rows = conn.execute(stmt)
            for cnt, row in enumerate(rows):
                print(row)
                region = Region(row.sr_name, row.cs_name, 'linear', row.sr_length, row.cs_version)
                location = Location(start, end, row.sr_length)
                st = Strand.REVERSE if strand == -1 else Strand.FORWARD
                slices.append(Slice(location, region, st))
            if cnt <= 0:
                raise NoResultFound(f'Could not find any region with name {seq_region_name} and coordinate {coord_system}')
            return slices
    