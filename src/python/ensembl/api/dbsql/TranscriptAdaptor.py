from sqlalchemy import and_
from sqlalchemy.orm import Session, Bundle
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.engine.row import Row

from ensembl.core.models import Transcript as TranscriptORM, TranscriptAttrib as TranscriptAttribORM, AttribType as AttribTypeORM
from ensembl.core.models import Biotype, Xref, Analysis

from typing import Union, List

from ensembl.api.core import Transcript, Location, Strand, ValueSetMetadata
from ensembl.api.dbsql.SliceAdaptor import SliceAdaptor
from enum import Enum
import warnings

__all__ = ['TranscriptAdaptor', 'TSLVERSION']

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
        Arg [1]    : str stable_id
            the stable id of the exon to retrieve
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
        Status     : Alpha
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
                Bundle("Biotype",
                       Biotype.so_acc,
                       Biotype.so_term
                ),
                Bundle("Xref",
                       Xref.dbprimary_acc,
                       Xref.description
                ),
                Bundle("Analysis",
                       Analysis.logic_name
                )
            )
            .join(TranscriptORM.display_xref)
            .join(Biotype, and_(TranscriptORM.biotype == Biotype.name, Biotype.object_type == 'transcript'))
            .join(TranscriptORM.analysis)
            .filter(TranscriptORM.stable_id == unversioned_stable_id)
            .filter(TranscriptORM.is_current == '1')
            .first()
        )
        
        # res = (session.execute(stmt)
        # .first())

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
        return TranscriptAdaptor._transcriptrow_to_transcript(session, res, transcript_attribs)

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
    def _transcriptrow_to_transcript(cls, session: Session, tr_row: Row, transcript_attribs: List[Row]) -> Transcript:
        slice = SliceAdaptor.fetch_by_seq_region_id(session, tr_row.Transcript.seq_region_id)
        sr_len = tr_row.Transcript.seq_region_end - tr_row.Transcript.seq_region_start
        slice.location = Location(tr_row.Transcript.seq_region_start, tr_row.Transcript.seq_region_end, sr_len)
        slice.strand = Strand.REVERSE if tr_row.Transcript.seq_region_strand == -1 else Strand.FORWARD

        tr = Transcript(
            '.'.join((str(tr_row.Transcript.stable_id), str(tr_row.Transcript.version))),
            tr_row.Xref.dbprimary_acc,
            slice
        )

        tr.add_metadata('source', ValueSetMetadata('transcript.source', tr_row.Transcript.source))
        tr.add_metadata('biotype', ValueSetMetadata(f"{tr_row.Transcript.biotype.lower()}", tr_row.Biotype.so_term))

        tr.add_metadata('description', ValueSetMetadata('coredb.description', tr_row.Transcript.description))
        tr.add_metadata('analysis', ValueSetMetadata('coredb.analysis', tr_row.Analysis.logic_name))
        tr.add_metadata('created_date', ValueSetMetadata('coredb.created_date', tr_row.Transcript.created_date))
        tr.add_metadata('modified_date', ValueSetMetadata('coredb.modified_date', tr_row.Transcript.modified_date))

        for ta in transcript_attribs:
            if ta.AttribType.code.lower() == 'gencode_basic':
                tr.add_metadata('gencode_basic', ValueSetMetadata('gencode_basic.true', 'true', f'{ta.TranscriptAttrib.value}'))
                continue
            if ta.AttribType.code.lower() == 'is_canonical':
                tr.add_metadata('canonical', ValueSetMetadata('canonical_transcript.true', 'true', f'{ta.AttribType.name}'))
                continue
            if ta.AttribType.code.lower() == 'tsl':
                tslv = TSLVERSION(ta.TranscriptAttrib.value[:5])
                tr.add_metadata('tsl', ValueSetMetadata(f'tsl.{tslv.name}', tslv.value, f'{ta.TranscriptAttrib.value}'))
                continue
            if 'mane' in ta.AttribType.code.lower():
                if ta.AttribType.code == 'MANE_Select':
                    # ncbi_t['ncbi_transcript'] = {'id':f'{ta.TranscriptAttrib.value}', 'url':f"https://www.ncbi.nlm.nih.gov/nuccore/{ta.TranscriptAttrib.value}"}
                    tr.add_metadata('mane.select', ValueSetMetadata('mane.select', f'{ta.TranscriptAttrib.value}', f'{ta.AttribType.name}'))
                    # value: select, label: MANE Select, ncbi_transcript: {id:"NM_1234", url:"<URL>"}
                    continue
            if ta.AttribType.code.lower() == 'appris':
                tr.add_metadata('appris', ValueSetMetadata('appris', f'{ta.TranscriptAttrib.value}', f'{ta.AttribType.name} {ta.TranscriptAttrib.value}'))
                # value: alternative1, label: APPRIS alternative:1
            if ta.AttribType.code.lower() == 'ccds_transcript':
                tr.add_metadata('ccds', ValueSetMetadata('ccds_transcript.true', f'{ta.TranscriptAttrib.value}', f'{ta.TranscriptAttrib.value}'))
                continue
            # BUG - this below must be a list!!!
            tr.add_metadata('transcript.attrib', ValueSetMetadata(ta.AttribType.code.lower(), f'{ta.TranscriptAttrib.value}', f'{ta.AttribType.name}'))
        return tr
