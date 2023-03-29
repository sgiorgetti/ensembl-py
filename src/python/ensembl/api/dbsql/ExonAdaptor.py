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

from sqlalchemy import select, and_
from sqlalchemy.orm import Session, Bundle
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.engine.row import Row

from ensembl.core.models import Exon as ExonORM, ExonTranscript as ExonTranscriptORM, Transcript as TranscriptORM

from typing import Union

from ensembl.api.core import Exon, SplicedExon, Slice, Strand, Transcript
from ensembl.api.dbsql.SliceAdaptor import SliceAdaptor

# from ensembl.api.dbsql.Utils import timeme

__all__ = [ 'ExonAdaptor' ]

class ExonAdaptor():
    """Contains all the exon related functions over Exon ORM
    """

    @classmethod
    def fetch_by_stable_id(cls, session: Session, stable_id: str) -> Exon:
        """
        Arg [1]    : sqlalchemy.orm.Session session
        Arg [2]    : str stable_id
                     The stable id of the exon to retrieve
        Example    : exon = ExonAdaptor.fetch_by_stable_id(dbc, 'ENSE00001544499')
        Description: Retrieves an Exon from the database via its stable id
                     The stable id can be versioned or unversioned
        Returntype : ensembl.api.core.Exon in native coordinates.
        Exceptions : NoResultFound
        Caller     : general
        Status     : Alpha
        """
        if stable_id.find('.') > 0:
            (unversioned_stable_id, version) = stable_id.split('.')
            return cls.fetch_by_stable_id_version(session, unversioned_stable_id, version)

        stmt = select(ExonORM).where(ExonORM.stable_id == stable_id, ExonORM.is_current == 1)
        res = (session.execute(stmt)
        .first())
        if not res:
                raise NoResultFound(f'Could not find {stable_id} as current stable_id')
        
        slice = SliceAdaptor.fetch_by_seq_region_id(session, res.Exon.seq_region_id)
        return cls._exonrow_to_exon(slice, res)

    @classmethod
    def fetch_by_stable_id_version(cls, session: Session, unversioned_stable_id: str, version: Union[str, int]) -> Exon:
        """
        Arg [1]    : sqlalchemy.orm.Session session
        Arg [2]    : str stable_id
                     The stable id of the exon to retrieve
        Example    : exon = ExonAdaptor.fetch_by_stable_id(dbc, 'ENSE00001544499')
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

        res = (session.execute(stmt)
        .first())
        if not res:
                raise NoResultFound(f'Could not find {unversioned_stable_id} v.{version} as current stable_id')
        slice = SliceAdaptor.fetch_by_seq_region_id(session, res.Exon.seq_region_id)
        return cls._exonrow_to_exon(slice, res)

    @classmethod
    def fetch_all_by_Transcript(cls, session: Session, transcript: Transcript) -> list[SplicedExon]:
        """
        Arg [1]    : sqlalchemy.orm.Session session
        Arg [2]    : ensembl.api.core.Transcript transcript
        Example    : none
        Description: Retrieves all Exons for the Transcript in 5-3 order
        Returntype : List[ensembl.api.core.SplicedExon] on Transcript slice 
        Exceptions : throws if transcript is not specified or None
        Caller     : Transcript->get_all_Exons()
        Status     : Alpha
        """
        if not transcript:
            raise ValueError("Transcript must be specified!")
        
        exon_rows = (
             session.query(Bundle("Transcript",
                            TranscriptORM.source
                           ),
                           Bundle("ExonTranscript", ExonTranscriptORM.rank),
                           ExonORM)
                .join(TranscriptORM.exons)
                .join(ExonTranscriptORM.exon)
                .filter(and_(TranscriptORM.transcript_id == transcript.internal_id, 
                             TranscriptORM.is_current == 1))
                .order_by(ExonTranscriptORM.rank)
                .all()
        )
        exons: list[SplicedExon] = []
        for er in exon_rows:
            exons.append(cls._exonrow_to_splicedexon(transcript._slice, er))
        return exons
            

    @classmethod
    def fetch_by_internal_id(cls, session: Session, internal_id: int) -> Exon:
        """
        Arg [1]    : sqlalchemy.orm.Session session
        Arg [2]    : internal_id: int
                     The internal id (dbID) of the exon to retrieve
        Example    : exon = ExonAdaptor.fetch_by_internal_id(session, 1234)
        Description: Retrieves an Exon from the database via its stable id
        Returntype : ensembl.api.core.Exon in native coordinates.
        Exceptions : NoResultFound, ValueError
        Caller     : general
        Status     : Alpha
        """
        if internal_id is None or not isinstance(internal_id, int):
             return None
        
        stmt = (select(ExonORM)
        .where(ExonORM.exon_id == internal_id)
        .where(ExonORM.is_current == 1)
        )

        res = (session.execute(stmt)
        .first())
        if not res:
                raise NoResultFound(f'Could not find exon with exon_id {internal_id}.')
        slice = SliceAdaptor.fetch_by_seq_region_id(session, res[0].seq_region_id)
        return cls._exonrow_to_exon(slice, res)
    

    @classmethod
    def _exonrow_to_exon(cls, exon_slice: Slice, row: Row) -> Exon:
        e = Exon(
             row.Exon.stable_id,
             row.Exon.version,
             row.Exon.phase,
             row.Exon.end_phase,
             row.Exon.exon_id,
             exon_slice,
             row.Exon.seq_region_start,
             row.Exon.seq_region_end,
             Strand(row.Exon.seq_region_strand),
             None,
             True if row.Exon.is_constitutive else False,
             True if row.Exon.is_current else False,
             created_date=row.Exon.created_date,
             modified_date=row.Exon.modified_date
            )
        return e
    
    @classmethod
    def _exonrow_to_splicedexon(cls, tr_slice: Slice, row: Row) -> SplicedExon:
        e = SplicedExon(
             row.Exon.stable_id,
             row.Exon.version,
             row.Exon.phase,
             row.Exon.end_phase,
             row.ExonTranscript.rank,
             row.Exon.exon_id,
             tr_slice,
             row.Exon.seq_region_start,
             row.Exon.seq_region_end,
             Strand(row.Exon.seq_region_strand),
             None,
             row.Exon.is_constitutive,
             row.Exon.is_current,
             row.Transcript.source,
             created_date=row.Exon.created_date,
             modified_date=row.Exon.modified_date
            )

        return e
