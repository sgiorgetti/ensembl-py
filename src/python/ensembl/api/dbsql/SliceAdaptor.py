from sqlalchemy import select, func, text, and_
from sqlalchemy.orm import Session, Bundle
from sqlalchemy.orm.exc import NoResultFound

from ensembl.core.models import SeqRegion as SeqRegionORM, CoordSystem as CoordSystemORM, AttribType as AttribTypeORM
from ensembl.core.models import SeqRegionAttrib as SeqRegionAttribORM

from typing import List, Optional

from ensembl.api.dbsql.Exceptions import ArgumentError
from ensembl.api.dbsql.AssemblyAdaptor import AssemblyAdaptor
from ensembl.api.core import Slice, Region, Location, Strand, Assembly

__all__ = ['SliceAdaptor']

class SliceAdaptor():
    """Contains all the slice related functions over Slice ORM
    """
    
    @classmethod
    def fetch_by_seq_region(cls,
                            session: Session,
                            seq_region_name: str,
                            coord_system: str,
                            version: Optional[str] = None,
                            start: int = 1,
                            end: Optional[int] = None,
                            strand: int = 1
                            ) -> List[Slice]:
        """
        Arg [1]    : str seq_region_name
            the seq_region name as specified in the core DB
        Arg [2]    : str coord_system
            the coordinate system name of the seq_region to retrieve
        Arg [3]    : str version
            the coordinate system version as specified in the core DB - e.g. GRCh38
        Arg [4]    : int start
            the seq_region start position
        Arg [5]    : int end
            the seq_region end position
        Arg [6]    : int strand
            the seq_region strand +1 forward, -1 reverse

        Example    : region = region_adaptor->fetch_by_region( 'chromosome', 'X' );
        Description: Retrieves a slice on the requested region.  At a minimum the
                     name the name of the seq_region to fetch, alomg with the coordinate
                     system name, must be provided.
        Returntype : List[ensembl.api.core.Slice]
        Exceptions : NoResultFound, ArgumentError
        Caller     : general
        Status     : Alpha
        """
        if not session:
            raise ArgumentError()
        
        slices = []
        assert isinstance(start, int)
        if end is not None:
            assert isinstance(end, int)
        if not seq_region_name:
            raise ArgumentError('seq_region_name argument is required')
        if not coord_system:
            raise ArgumentError('coord_system argument is required')
        
        stmt = (select(
                Bundle("SeqRegion",
                        SeqRegionORM.seq_region_id,
                        SeqRegionORM.name,
                        SeqRegionORM.length
                ),
                Bundle("CoordSystem",
                        CoordSystemORM.coord_system_id,
                        CoordSystemORM.name,
                        CoordSystemORM.version,
                        CoordSystemORM.rank
                ))
        .join(SeqRegionORM.coord_system)
        .where(and_(CoordSystemORM.name == coord_system, SeqRegionORM.name == seq_region_name))
        )
        if version:
            stmt = stmt.filter(CoordSystemORM.version == version)
        
        rows = session.execute(stmt)
        if rows.rowcount == 0:
            raise NoResultFound(f'Could not find any region with name {seq_region_name} and coordinate {coord_system}')
        
        for row in rows:
            topology = 'CIRCULAR' if cls._is_circular(session, row.SeqRegion.seq_region_id) else 'LINEAR'
            assembly = AssemblyAdaptor.fetch_by_coord_system_id(session, row.CoordSystem.coord_system_id)
            region = Region(row.SeqRegion.name,
                            row.CoordSystem.name,
                            is_top_level = cls._is_top_level(session, row.SeqRegion.seq_region_id),
                            rank = cls._karyotype_rank(session, row.SeqRegion.seq_region_id),
                            topology=topology,
                            length=row.SeqRegion.length,
                            assembly=assembly)
            location = Location(start, end, row.SeqRegion.length)
            st = Strand.REVERSE if strand == -1 else Strand.FORWARD
            slices.append(Slice(location, region, st))
            
            return slices
    
    @classmethod
    def fetch_by_seq_region_id(cls, 
                               session: Session,
                               seq_region_id: int
                               ) -> Slice:
        if not session:
            raise ArgumentError()
        
        if not isinstance(seq_region_id, int) or seq_region_id <= 0:
            raise ArgumentError()

        # stmt = text('SELECT seq_region.seq_region_id, seq_region.name AS sr_name, seq_region.length AS sr_length,'
        #                 ' coord_system.name AS cs_name, coord_system.version AS cs_version, coord_system.rank AS cs_rank'
        #                 ' FROM seq_region JOIN coord_system ON seq_region.coord_system_id = coord_system.coord_system_id'
        #                 ' WHERE lower(coord_system.name)=:cs AND lower(seq_region.name)=:sr')

        result = (session.query(
                Bundle("SeqRegion",
                        SeqRegionORM.seq_region_id,
                        SeqRegionORM.name,
                        SeqRegionORM.length
                ),
                Bundle("CoordSystem",
                        CoordSystemORM.coord_system_id,
                        CoordSystemORM.name,
                        CoordSystemORM.version,
                        CoordSystemORM.rank
                ))
        .join(SeqRegionORM.coord_system)
        .where(SeqRegionORM.seq_region_id == seq_region_id)
        .first())

        topology = 'CIRCULAR' if cls._is_circular(session, result.SeqRegion.seq_region_id) else 'LINEAR'
        assembly = AssemblyAdaptor.fetch_by_coord_system_id(session, result.CoordSystem.coord_system_id)
        region = Region(result.SeqRegion.name,
                            result.CoordSystem.name,
                            is_top_level = cls._is_top_level(session, result.SeqRegion.seq_region_id),
                            rank = cls._karyotype_rank(session, result.SeqRegion.seq_region_id),
                            topology=topology,
                            length=result.SeqRegion.length,
                            assembly=assembly)
        # location = Location(start, end, row.sr_length)
        # st = Strand.REVERSE if strand == -1 else Strand.FORWARD
        return Slice(region=region)
        
        

    @classmethod
    def _is_circular(cls, session: Session, seq_region_id: int) -> bool:
        # stmt = (select(func.count(SeqRegionAttribORM.seq_region_id))
        # .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
        # .where(AttribTypeORM.code == 'circular_seq')
        # .where(SeqRegionAttribORM.seq_region_id == seq_region_id)
        # )
        # res = session.execute(stmt).first()
        res = cls._get_seq_region_attrib(session, seq_region_id, 'circular_seq')
        if res:
            return True
        return False
    
    @classmethod
    def _is_top_level(cls, session: Session, seq_region_id: int) -> bool:
        res = cls._get_seq_region_attrib(session, seq_region_id, 'toplevel')
        if res:
            return True
        return False

    @classmethod
    def _karyotype_rank(cls, session: Session, seq_region_id: int) -> int:
        return cls._get_seq_region_attrib(session, seq_region_id, 'karyotype_rank')

    @classmethod
    def _get_seq_region_attrib(cls, session: Session, seq_region_id: int, seq_region_attrib_code: str) -> str:
        res = (session.query(SeqRegionAttribORM.value)
        .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
        .where(AttribTypeORM.code == seq_region_attrib_code)
        .where(SeqRegionAttribORM.seq_region_id == seq_region_id)
        .first())
        if res:
            return res
        return ""
    