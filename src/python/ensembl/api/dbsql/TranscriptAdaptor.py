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

from sqlalchemy.orm import Session, Bundle
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.engine.row import Row

from ensembl.core.models import Transcript as TranscriptORM, TranscriptAttrib as TranscriptAttribORM, AttribType as AttribTypeORM, Gene as GeneORM
from ensembl.core.models import Xref, Analysis

from typing import Union

from ensembl.api.core import Transcript, Strand, Gene
from ensembl.api.dbsql.SliceAdaptor import SliceAdaptor
from ensembl.api.dbsql.BiotypeAdaptor import BiotypeAdaptor
from enum import Enum
import warnings

__all__ = [ 'TranscriptAdaptor', 'TSLVERSION' ]

class TSLVERSION(Enum):
    one = "tsl1"
    two = "tsl2"
    three = "tsl3"
    four = "tsl4"
    five = "tsl5"
    na = "tslNA"

class TranscriptAdaptor():
    """Contains all the Transcript related functions over Transcript ORM
    """

    @classmethod
    def fetch_by_stable_id(cls, session: Session, stable_id: str) -> Transcript:
        """
        Arg [1]    : Session session - the ORM session object to connect to the DB to
        Arg [2]    : str stable_id
                     the stable id of the transcript to retrieve
        Example    : transcript = TranscriptAdaptor.fetch_by_stable_id('ENST00000309301');
                     transcript = TranscriptAdaptor.fetch_by_stable_id('ENST00000309301.2');
        Description: Retrieves a transcript via its stable id.
                     The stable id can be versioned or unversioned.
                     The transcript will be retrieved in its native coordinate system (i.e.
                     in the coordinate system it is stored in the database). It may
                     be converted to a different coordinate system through a call to
                     transform() or transfer(). If the transcript is not found
                     None is returned instead.
        Returntype : ensembl.api.core.Transcript in native coordinates.
        Exceptions : NoResultFound
        Caller     : general
        Status     : At Risk
                   : under development
        """
        if stable_id.find('.') > 0:
            (unversioned_stable_id, version) = stable_id.split('.')
        else:
            unversioned_stable_id = stable_id
            version = 0

        int(version)

        res = (
            session.query(
                TranscriptORM,
                Bundle("Xref",
                       Xref.dbprimary_acc,
                       Xref.description
                ),
                Bundle("Analysis",
                       Analysis.logic_name
                )
            )
            .join(TranscriptORM.display_xref)
            .join(TranscriptORM.analysis)
            .filter(TranscriptORM.stable_id == unversioned_stable_id)
            .filter(TranscriptORM.is_current == '1')
            .first()
        )

        if not res:
                raise NoResultFound(f'Could not find {unversioned_stable_id} as current stable_id')
        if version > 0 and res.Transcript.version != version:
            warnings.warn(f"Could not find {stable_id}. Found {unversioned_stable_id}.{res.Transcript.version} instead!", UserWarning)
        
        transcript_attribs = (
            session.query(TranscriptAttribORM, AttribTypeORM)
            .join(TranscriptAttribORM.attribs)
            .filter(TranscriptAttribORM.transcript_id == res.Transcript.transcript_id)
            .all()
        )
        return cls._transcriptrow_to_transcript(session, res, transcript_attribs)


    @classmethod
    def fetch_by_stable_id_version(cls, session: Session, unversioned_stable_id: str, version: Union[str, int]) -> Transcript:
        """
        Arg [1]    : str stable_id
                     the stable id of the transcript to retrieve
        Arg [2]    : str/int version
                     The version of the stable_id to retrieve
        Example    : transcript = TranscriptAdaptor.fetch_by_stable_id_version('ENST00000309301', 3);
        Description: Retrieves a transcript object from the database via its 
                     stable id and version.
                     The transcript will be retrieved in its native coordinate system (i.e.
                     in the coordinate system it is stored in the database). It may
                     be converted to a different coordinate system through a call to
                     transform() or transfer(). If the transcript is not found
                     None is returned instead.
        Returntype : ensembl.api.core.Transcript in native coordinates.
        Exceptions : NoResultFound, ValueError
        Caller     : general
        Status     : Alpha
        """
        int(version)
        return cls.fetch_by_stable_id(session, f"{unversioned_stable_id}.{version}")

    
    @classmethod
    def fetch_all_by_gene_id(cls, session: Session, gene_internal_id: int) -> tuple[Transcript]:
        """
        Arg [1]    : Session session - the ORM session object to connect to the DB to
        Arg [2]    : int gene_internal_id
                     the gene internal id of the transcripts to retrieve
        Example    : transcript = TranscriptAdaptor.fetch_all_by_gene_id(1234);
        Description: Retrieves a transcript list via the related gene internal_id.
                     The transcripts will be retrieved in its native coordinate system (i.e.
                     in the coordinate system it is stored in the database).
                     If no transcript is found, an empty list is returned.
        Returntype : list[ensembl.api.core.Transcript] in native coordinates.
        Exceptions : none
        Caller     : general
        Status     : At Risk
                   : under development
        """
        rows = (
            session.query(
                TranscriptORM,
                Bundle("Xref",
                       Xref.dbprimary_acc,
                       Xref.description
                ),
                Bundle("Analysis",
                       Analysis.logic_name
                )
            )
            .join(TranscriptORM.display_xref)
            .join(TranscriptORM.analysis)
            .filter(TranscriptORM.gene_id == gene_internal_id)
            .filter(TranscriptORM.is_current == '1')
            .all()
        )

        if not rows:
            return []
        
        transcripts: list[Transcript] = []
        for tr_row in rows:
            transcript_attribs = (
                session.query(TranscriptAttribORM, AttribTypeORM)
                .join(TranscriptAttribORM.attribs)
                .filter(TranscriptAttribORM.transcript_id == tr_row.Transcript.transcript_id)
                .all()
            )
            tr = cls._transcriptrow_to_transcript(session, tr_row, transcript_attribs)
            transcripts.append(tr)

        return tuple(transcripts)
    
    
    @classmethod
    def fetch_all_by_gene(cls, session: Session, gene: Gene) -> tuple[Transcript]:
        return cls.fetch_all_by_gene_id(session, gene.internal_id)


    @classmethod
    def _transcriptrow_to_transcript(cls, session: Session, tr_row: Row, transcript_attribs: list[Row]) -> Transcript:
        slice = SliceAdaptor.fetch_by_seq_region_id(session, tr_row.Transcript.seq_region_id)
        canonical_tr_id = cls._fetch_canonical_transcript_id(session, tr_row.Transcript.gene_id)
        biotype = BiotypeAdaptor.fetch_by_name_object_type(session, tr_row.Transcript.biotype, 'transcript')

        tr = Transcript(
            tr_row.Transcript.stable_id,
            tr_row.Transcript.version,
            tr_row.Transcript.transcript_id,
            slice,
            tr_row.Transcript.seq_region_start,
            tr_row.Transcript.seq_region_end,
            Strand(tr_row.Transcript.seq_region_strand),
            tr_row.Analysis.logic_name,
            exons=None,
            biotype=biotype,
            source=tr_row.Transcript.source,
            is_canonical=True if tr_row.Transcript.transcript_id == canonical_tr_id else False,
            external_name=tr_row.Xref.dbprimary_acc,
            description=tr_row.Transcript.description,
            created_date=tr_row.Transcript.created_date,
            modified_date=tr_row.Transcript.modified_date
        )

        for ta in transcript_attribs:
            tr.add_attrib(ta.AttribType.code, ta.TranscriptAttrib.value)

        return tr
    
    
    @classmethod
    def _fetch_canonical_transcript_id(cls, session: Session, gene_internal_id: int) -> int:
        res = (
            session.query(
                GeneORM.canonical_transcript_id
            )
            .where(GeneORM.gene_id == gene_internal_id)
            .first()
        )
        if not res:
            return None
        return res.canonical_transcript_id
