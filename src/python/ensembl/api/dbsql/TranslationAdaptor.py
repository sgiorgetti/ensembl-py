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

from sqlalchemy import select
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.engine.row import Row

from ensembl.core.models import Translation as TranslationORM, TranslationAttrib as TranslationAttribORM, AttribType as AttribTypeORM, Transcript as TranscriptORM

from ensembl.api.core import Translation, Transcript
from ensembl.api.dbsql.ExonAdaptor import ExonAdaptor


__all__ = [ 'TranslationAdaptor' ]

class TranslationAdaptor():
    """Contains all the Translation related functions over Translation ORM
    """

    @classmethod
    def fetch_by_Transcript(cls, session: Session, transcript: Transcript):
        """
        Arg [1]    : session: Session - the ORM session object to connect to the DB to
        Arg [2]    : transcript: ensembl.api.core.Transcript
        Example    : tl = TranslationAdaptor.fetch_by_Transcript(transcript);
        Description: Retrieves a Translation via its associated transcript.
                     If the Translation is not found, undef is returned.
        Returntype : ensembl.api.core.Translation
        Exceptions : throw on incorrect argument
        Caller     : Transcript
        Status     : At Risk
                   : under development
        """
        if not transcript:
            raise AttributeError(f'You need to specify a transcript.')
        if not isinstance(transcript, Transcript):
            raise AttributeError(f'You need to specify a transcript object, instead of {type(transcript)}.')
        
        tr_id = transcript.canonical_translation_id
        if tr_id is None:
            return None
        
        res = (
            session.query(
                TranslationORM
            )
            .filter(TranslationORM.translation_id == tr_id)
            .first()
        )

        if not res:
            raise NoResultFound(f'Could not find any translation with {tr_id}.')
        
        translation_attribs = (
            session.query(TranslationAttribORM, AttribTypeORM)
            .join(TranslationAttribORM.attribs)
            .filter(TranslationAttribORM.translation_id == tr_id)
            .all()
        )
        tl = cls._translationrow_to_translation_fast(transcript, res, translation_attribs)
        transcript.set_translation(tl)

        return tl
    
    @classmethod
    def fetch_by_internal_id(cls, session: Session, internal_id: int) -> Translation:
        """
        Arg [1]    : Session session - the ORM session object to connect to the DB to
        Arg [2]    : internal_id: int
                     the internal id (dbID) of the transcript to retrieve
        Example    : translation = TranslationAdaptor.fetch_by_internal_id(1234);
        Description: Retrieves a translation via its internal id (dbID).
                     The translation will be retrieved in its native coordinate system (i.e.
                     in the coordinate system it is stored in the database).
        Returntype : ensembl.api.core.Translation in native coordinates.
        Exceptions : NoResultFound
        Caller     : general
        Status     : At Risk
                   : under development
        """
        if internal_id is None:
            return None

        res = (
            session.query(
                TranslationORM
            )
            .filter(TranslationORM.translation_id == internal_id)
            .first()
        )

        if not res:
            raise NoResultFound(f'Could not find any translation with {internal_id}.')
        
        translation_attribs = (
            session.query(TranslationAttribORM, AttribTypeORM)
            .join(TranslationAttribORM.attribs)
            .filter(TranslationAttribORM.translation_id == res.translation_id)
            .all()
        )
        return cls._translationrow_to_translation(session, res, translation_attribs)
    

    @classmethod
    def fetch_by_transcript_id(cls, session: Session, transcript_id: int) -> Translation:
        """
        Arg [1]    : Session session - the ORM session object to connect to the DB to
        Arg [2]    : transcript_id: int
                     the internal id (dbID) of the translation-related transcript to retrieve
        Example    : translation = TranslationAdaptor.fetch_by_transcript_id(1234);
        Description: Retrieves a translation via the transcript's internal id (dbID).
                     The translation will be retrieved in its native coordinate system (i.e.
                     in the coordinate system it is stored in the database).
        Returntype : ensembl.api.core.Translation in native coordinates.
        Exceptions : NoResultFound
        Caller     : general
        Status     : At Risk
                   : under development
        """
        if transcript_id is None:
            return None

        res = (
            session.query(
                TranslationORM
            )
            .filter(TranslationORM.transcript_id == transcript_id)
            .first()
        )

        if not res:
            raise NoResultFound(f'Could not find any translation realted to transcript {transcript_id}.')
        
        stmt = (
            select(TranslationAttribORM, AttribTypeORM)
            .join(TranslationAttribORM.attribs)
            .filter(TranslationAttribORM.translation_id == res.Translation.translation_id)
        )
        translation_attribs = session.scalars(stmt).all()
        return cls._translationrow_to_translation(session, res, translation_attribs)


    @classmethod
    def _translationrow_to_translation(cls, session: Session, tr_row: Row, translation_attribs: list[Row]) -> Translation:
        start_exon = ExonAdaptor.fetch_by_internal_id(session, tr_row.start_exon_id)
        end_exon = ExonAdaptor.fetch_by_internal_id(session, tr_row.end_exon_id)
        if not start_exon or not end_exon:
             raise NoResultFound(
                 f'Could not find either start exon or end exon for translation id {tr_row.translation_id}'
                )
        tr = Translation(
            start_exon,
            end_exon,
            tr_row.seq_start,
            tr_row.seq_end,
            tr_row.stable_id,
            tr_row.version,
            tr_row.translation_id,
            created_date=tr_row.created_date,
            modified_date=tr_row.modified_date)

        for ta in translation_attribs:
            tr.add_attrib(ta.AttribType.code, ta.TranslationAttrib.value)

        return tr
    
    @classmethod
    def _translationrow_to_translation_fast(cls, t: Transcript, tr_row: Row, translation_attribs: list[Row]) -> Translation:
        start_exon = None
        end_exon = None
        for e in t.get_exons():
            if start_exon is None and e.internal_id == tr_row.start_exon_id:
                start_exon = e
            if end_exon is None and e.internal_id == tr_row.end_exon_id:
                end_exon = e
            if start_exon is not None and end_exon is not None:
                break

        if not start_exon or not end_exon:
             raise NoResultFound(
                 f'Could not find either start exon or end exon for translation id {tr_row.translation_id}'
                )
        tr = Translation(
            start_exon,
            end_exon,
            tr_row.seq_start,
            tr_row.seq_end,
            tr_row.stable_id,
            tr_row.version,
            tr_row.translation_id,
            created_date=tr_row.created_date,
            modified_date=tr_row.modified_date)

        for ta in translation_attribs:
            tr.add_attrib(ta.AttribType.code, ta.TranslationAttrib.value)

        return tr
