__all__ = ['ExonAdaptor']

from sqlalchemy import select
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.engine.row import Row

from ensembl.core.models import Exon as ExonORM

from typing import Union, List

from ensembl.database.dbconnection import DBConnection
from ensembl.api.core import Exon, Location, Strand, ValueSetMetadata
from ensembl.api.dbsql.SliceAdaptor import SliceAdaptor

class ExonAdaptor():
    """Contains all the exon related functions over Exon ORM
    """

    @classmethod
    def fetch_by_stable_id(cls, dbconnection: DBConnection, stable_id: str) -> Exon:
        """
        Arg [1]    : str stable_id
            the stable id of the exon to retrieve
        Example    : exon = exon_adaptor->fetch_by_stable_id('ENSE0000988221');
        Description: Retrieves an Exon from the database via its stable id
                     The stable id can be versioned or unversioned
        Returntype : ensembl.api.core.Exon in native coordinates.
        Exceptions : NoResultFound
        Caller     : general
        Status     : Alpha
        """
        if stable_id.find('.') > 0:
            (unversioned_stable_id, version) = stable_id.split('.')
            return cls.fetch_by_stable_id_version(dbconnection, unversioned_stable_id, version)

        stmt = select(ExonORM).where(ExonORM.stable_id == stable_id, ExonORM.is_current == 1)
        with dbconnection.connect() as conn:
            res = (conn.execute(stmt)
            .first())
            if not res:
                 raise NoResultFound(f'Could not find {stable_id} as current stable_id')
            return res

    @classmethod
    def fetch_by_stable_id_version(cls, dbconnection: DBConnection, unversioned_stable_id: str, version: Union[str, int]) -> Exon:
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
        stmt = (select(ExonORM)
        .where(ExonORM.stable_id == unversioned_stable_id)
        .where(ExonORM.version == version)
        .where(ExonORM.is_current == 1)
        )
        with dbconnection.session_scope() as session:
            res = (session.execute(stmt)
            .first())
            if not res:
                 raise NoResultFound(f'Could not find {unversioned_stable_id} v.{version} as current stable_id')
            return ExonAdaptor._exonrow_to_exon(session, res[0])

    @classmethod
    def _exonrow_to_exon(cls, session: Session, row: Row) -> Exon:
        slice = SliceAdaptor.fetch_by_seq_region_id(session, row.seq_region_id)
        sr_len = row.seq_region_end - row.seq_region_start
        slice.location = Location(row.seq_region_start, row.seq_region_end, sr_len)
        slice.strand = Strand.REVERSE if row.seq_region_strand == -1 else Strand.FORWARD
        e = Exon('.'.join((str(row.stable_id), str(row.version))), slice, row.phase, row.end_phase)
        e.add_metadata('is_current', ValueSetMetadata('exon.is_current', row.is_current))
        e.add_metadata('is_constitutive', ValueSetMetadata('exon.is_constitutive', row.is_constitutive))
        e.add_metadata('created_date', ValueSetMetadata('exon.created_date', row.created_date))
        e.add_metadata('modified_date', ValueSetMetadata('exon.modified_date', row.modified_date))
        return e
