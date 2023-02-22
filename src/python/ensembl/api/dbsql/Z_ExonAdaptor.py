__all__ = ['ExonAdaptor']

from sqlalchemy import select
from sqlalchemy.schema import Table
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.engine.row import Row

from typing import Any, Union, List, Dict

from ensembl.database.dbconnection import DBConnection
from ensembl.api.core import Exon, Slice

class ExonAdaptor():
    """Contains all the exon related functions over Exon ORM

    Attributes:
        dbc: DBConnection
        exon: Table object
    """
    def __init__(self, dbconnection:DBConnection=None) -> None:
        if not dbconnection:
            Exception()
        self._dbc = dbconnection
        self._exon = Table("exon", dbconnection._metadata, autoload_with=dbconnection._engine)
        

    def __repr__(self) -> str:
        """Returns a string representation of this object."""
        return f'{self.__class__.__name__}(Table({self._exon})@{self._dbc._engine})'

    def get_dbconnection(self):
        return self._dbc

    def fetch_by_stable_id(self, stable_id: str) -> Any:
        """
        Arg [1]    : str stable_id
            the stable id of the exon to retrieve
        Example    : exon = exon_adaptor->fetch_by_stable_id('ENSE0000988221');
        Description: Retrieves an Exon from the database via its stable id
        Returntype : ensembl.api.core.Exon in native coordinates.
        Exceptions : NoResultFound
        Caller     : general
        Status     : Alpha
        """
        stmt = select(self._exon).where(self._exon.c.stable_id == stable_id, self._exon.c.is_current == 1)
        with self._dbc._engine.connect() as conn:
            res = (conn.execute(stmt)
            .first())
            if not res:
                 raise NoResultFound(f'Could not find {stable_id} as current stable_id')
            return res

    def fetch_by_stable_id_version(self, unversioned_stable_id: str, version: Union[str, int]) -> Any:
        """
        Arg [1]    : str stable_id
            the stable id of the exon to retrieve
        Example    : exon = exon_adaptor->fetch_by_stable_id('ENSE0000988221');
        Description: Retrieves an Exon from the database via its stable id
        Returntype : ensembl.api.core.Exon in native coordinates.
        Exceptions : NoResultFound, ValueError
        Caller     : general
        Status     : Alpha
        """
        int(version)
        stmt = (select(self._exon)
        .where(self._exon.c.stable_id == unversioned_stable_id)
        .where(self._exon.c.version == version)
        .where(self._exon.c.is_current == 1)
        )
        with self._dbc._engine.connect() as conn:
            res = (conn.execute(stmt)
            .first())
            if not res:
                 raise NoResultFound(f'Could not find {unversioned_stable_id} v.{version} as current stable_id')
            return res

    def _exonrow_to_exon(self, row: Row) -> Exon:
        pass


    def _fetch_slice(self) -> Slice:
        pass