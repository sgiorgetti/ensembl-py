__all__ = ['SliceAdaptor']

from sqlalchemy import select, func, text
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound

from ensembl.core.models import SeqRegion as SeqRegionORM, CoordSystem as CoordSystemORM, AttribType as AttribTypeORM
from ensembl.core.models import SeqRegionAttrib as SeqRegionAttribORM

from typing import List, Optional

from ensembl.database.dbconnection import DBConnection
from ensembl.api.dbsql.DBAdaptor import ArgumentError
from ensembl.api.core import Slice, Region, Location, Strand

class SliceAdaptor():
    """Contains all the slice related functions over Slice ORM
    """
    
    @classmethod
    def fetch_by_seq_region(cls,
                            dbconnection: DBConnection,
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
        if not dbconnection:
            raise ArgumentError()
        
        slices = []
        assert isinstance(start, int)
        if end is not None:
            assert isinstance(end, int)
        if not seq_region_name:
            raise ArgumentError('seq_region_name argument is required')
        if not coord_system:
            raise ArgumentError('coord_system argument is required')
        
        if not version:
            stmt = text('SELECT seq_region.seq_region_id, seq_region.name AS sr_name, seq_region.length AS sr_length,'
                        ' coord_system.name AS cs_name, coord_system.version AS cs_version, coord_system.rank AS cs_rank'
                        ' FROM seq_region JOIN coord_system ON seq_region.coord_system_id = coord_system.coord_system_id'
                        ' WHERE lower(coord_system.name)=:cs AND lower(seq_region.name)=:sr')
            params = {"cs": coord_system, "sr": seq_region_name}
        else:
            stmt = text('SELECT seq_region.seq_region_id, seq_region.name AS sr_name, seq_region.length AS sr_length,'
                        ' coord_system.name AS cs_name, coord_system.version AS cs_version, coord_system.rank AS cs_rank'
                        ' FROM seq_region JOIN coord_system ON seq_region.coord_system_id = coord_system.coord_system_id'
                        ' WHERE lower(coord_system.name) = :cs AND lower(seq_region.name) = :sr'
                        ' AND lower(coord_system.version) = :v')
            params = {"cs": coord_system, "sr": seq_region_name, "v": version}
        
        with dbconnection.session_scope() as session:
            rows = session.execute(stmt, params)
            if rows.rowcount == 0:
                raise NoResultFound(f'Could not find any region with name {seq_region_name} and coordinate {coord_system}')
            for row in rows:
                topology = 'CIRCULAR' if cls._is_circular(session, row.seq_region_id) else 'LINEAR'
                region = Region(row.sr_name, row.cs_name, topology, row.sr_length, row.cs_version)
                location = Location(start, end, row.sr_length)
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

        stmt = text('SELECT seq_region.seq_region_id, seq_region.name AS sr_name, seq_region.length AS sr_length,'
                        ' coord_system.name AS cs_name, coord_system.version AS cs_version, coord_system.rank AS cs_rank'
                        ' FROM seq_region JOIN coord_system ON seq_region.coord_system_id = coord_system.coord_system_id'
                        ' WHERE lower(coord_system.name)=:cs AND lower(seq_region.name)=:sr')
        stmt = (select(
            SeqRegionORM.seq_region_id,
            SeqRegionORM.name.label('sr_name'),
            SeqRegionORM.length.label('sr_length'),
            CoordSystemORM.name.label('cs_name'),
            CoordSystemORM.version.label('cs_version'),
            CoordSystemORM.rank.label('cs_rank')
        )
        .join(CoordSystemORM, SeqRegionORM.coord_system_id == CoordSystemORM.coord_system_id)
        .where(SeqRegionORM.seq_region_id == seq_region_id)
        )

        result = session.execute(stmt).first()
        topology = 'CIRCULAR' if cls._is_circular(session, result.seq_region_id) else 'LINEAR'
        region = Region(result.sr_name, result.cs_name, topology, result.sr_length, result.cs_version)
        # location = Location(start, end, row.sr_length)
        # st = Strand.REVERSE if strand == -1 else Strand.FORWARD
        return Slice(region=region)
        
        

    @classmethod
    def _is_circular(cls, session: Session, seq_region_id: int) -> bool:
        stmt = (select(func.count(SeqRegionAttribORM.seq_region_id))
        .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
        .where(AttribTypeORM.code == 'circular_seq')
        .where(SeqRegionAttribORM.seq_region_id == seq_region_id)
        )
        res = session.execute(stmt).first()
        if res[0] > 0:
            return True
        return False