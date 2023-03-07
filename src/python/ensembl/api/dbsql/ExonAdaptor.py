__all__ = ['ExonAdaptor']

from sqlalchemy import select, and_
from sqlalchemy.orm import Session, Bundle
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.engine.row import Row

from ensembl.core.models import Exon as ExonORM, ExonTranscript as ExonTranscriptORM, Transcript as TranscriptORM

from typing import Union, List

from ensembl.api.core import Exon, SplicedExon, Slice, Location, Strand, ValueSetMetadata, Transcript
from ensembl.api.dbsql.SliceAdaptor import SliceAdaptor

class ExonAdaptor():
    """Contains all the exon related functions over Exon ORM
    """

    @classmethod
    def fetch_by_stable_id(cls, session: Session, stable_id: str) -> Exon:
        """
        Arg [1]    : ensembl.database.dbconnection.DBConnection dbconnection
                     The DB Connection object
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
        return res

    @classmethod
    def fetch_by_stable_id_version(cls, session: Session, unversioned_stable_id: str, version: Union[str, int]) -> Exon:
        """
        Arg [1]    : ensembl.database.dbconnection.DBConnection dbconnection
                     The DB Connection object
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
        slice = SliceAdaptor.fetch_by_seq_region_id(session, res[0].seq_region_id)
        return ExonAdaptor._exonrow_to_exon(slice, res[0])

    @classmethod
    def fetch_all_by_Transcript(cls, session: Session, transcript: Transcript) -> List[SplicedExon]:
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
        
        tr_stable_id = transcript.unversioned_stable_id
        exon_rows = (
             session.query(Bundle("Transcript",
                            TranscriptORM.stable_id,
                            TranscriptORM.is_current
                           ),
                           Bundle("ExonTranscript", ExonTranscriptORM.rank),
                           ExonORM)
                .join(TranscriptORM.exons)
                .join(ExonTranscriptORM.exon)
                .filter(and_(TranscriptORM.stable_id == tr_stable_id, TranscriptORM.is_current == 1))
                .order_by(ExonTranscriptORM.rank)
                .all()
        )
        exons: List[SplicedExon] = []
        for er in exon_rows:
            # exons.append(ExonAdaptor._exonrow_to_splicedexon(transcript.get_slice(), er))
            exons.append(ExonAdaptor._exonrow_to_splicedexon(transcript, er))
        # transcript.set_exons(exons)
        return exons
            


    @classmethod
    def _exonrow_to_exon(cls, exon_slice: Slice, row: Row) -> Exon:
        sr_len = row.seq_region_end - row.seq_region_start
        exon_slice.location = Location(row.seq_region_start, row.seq_region_end, sr_len)
        exon_slice.strand = Strand.REVERSE if row.seq_region_strand == -1 else Strand.FORWARD
        e = Exon('.'.join((str(row.stable_id), str(row.version))), slice, row.phase, row.end_phase)
        e.add_metadata('is_current', ValueSetMetadata('exon.is_current', row.is_current))
        e.add_metadata('is_constitutive', ValueSetMetadata('exon.is_constitutive', row.is_constitutive))
        e.add_metadata('created_date', ValueSetMetadata('exon.created_date', row.created_date))
        e.add_metadata('modified_date', ValueSetMetadata('exon.modified_date', row.modified_date))
        return e
    
    @classmethod
    def _exonrow_to_splicedexon(cls, transcript: Transcript, row: Row) -> SplicedExon:
        slice = Slice(region=transcript.get_slice().region, strand=transcript.get_slice().strand)
        sr_len = row.Exon.seq_region_end - row.Exon.seq_region_start
        slice.location = Location(row.Exon.seq_region_start, row.Exon.seq_region_end, sr_len)
        strand = Strand.REVERSE if row.Exon.seq_region_strand == -1 else Strand.FORWARD
        if slice.strand != strand:
            raise Exception("STRAND!!!!") 
        e = SplicedExon('.'.join((str(row.Exon.stable_id), str(row.Exon.version))), slice, row.Exon.phase, row.Exon.end_phase, row.ExonTranscript.rank)
        e.add_metadata('source', ValueSetMetadata('exon.source', transcript.get_metadata('source').value))
        e.add_metadata('is_current', ValueSetMetadata('exon.is_current', row.Exon.is_current))
        e.add_metadata('is_constitutive', ValueSetMetadata('exon.is_constitutive', row.Exon.is_constitutive))
        e.add_metadata('created_date', ValueSetMetadata('exon.created_date', row.Exon.created_date))
        e.add_metadata('modified_date', ValueSetMetadata('exon.modified_date', row.Exon.modified_date))
        return e
