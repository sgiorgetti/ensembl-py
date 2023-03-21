from sqlalchemy import select, and_
from sqlalchemy.orm import Session, Bundle

from ensembl.core.models import SeqRegion as SeqRegionORM, CoordSystem as CoordSystemORM, AttribType as AttribTypeORM
from ensembl.core.models import SeqRegionAttrib as SeqRegionAttribORM

from typing import Optional

from ensembl.api.dbsql.Exceptions import ArgumentError
from ensembl.api.dbsql.CoordSystemAdaptor import CoordSystemAdaptor
from ensembl.api.core import Location, Strand, Slice, Region, RegionTopology, CoordSystem
import warnings

__all__ = [ 'SliceAdaptor' ]

class SliceAdaptor():
    """Contains all the slice related functions over Slice ORM
    """
    
    @classmethod
    def fetch_by_seq_region(cls,
                            session: Session,
                            seq_region_name: str,
                            coord_system: str, # this is the Region.code in CDM
                            version: Optional[str] = None, # this is Assembly.name in CDM
                            start: int = 1,
                            end: Optional[int] = None,
                            strand: int = 1
                            ) -> list[Slice]:
        """
        Arg [1]    : Session session
            the session to connect to a core DB
        Arg [2]    : str seq_region_name
            the seq_region name as specified in the core DB
        Arg [3]    : str coord_system
            the coordinate system name of the seq_region to retrieve
        Arg [4]    : str version
            the coordinate system version as specified in the core DB - e.g. GRCh38
        Arg [5]    : int start
            the seq_region start position
        Arg [6]    : int end
            the seq_region end position
        Arg [7]    : int strand
            the seq_region strand +1 forward, -1 reverse

        Example    : region = region_adaptor.fetch_by_region( session, 'chromosome', 'X' );
        Description: Retrieves a slice on the requested region.  At a minimum the
                     name the name of the seq_region to fetch, alomg with the coordinate
                     system name, must be provided.
        Returntype : List[ensembl.api.core.Slice]
        Exceptions : throw if no or wrong argument is provided
        Caller     : general
        Status     : At Risk
                   : under development
        """
        if not session:
            raise ArgumentError()
        
        slices = []
        if start is None:
            start = 1
        if not isinstance(start, int):
            raise ArgumentError('start argument must be an int')
        if end and not isinstance(end, int):
            raise ArgumentError('start argument must be an int or None')
        if not seq_region_name:
            raise ArgumentError('seq_region_name argument is required')
        if not coord_system:
            raise ArgumentError('coord_system argument is required')
        
        stmt = (select(
                Bundle("SeqRegion",
                        SeqRegionORM.seq_region_id,
                        SeqRegionORM.name.label('sr_name'),
                        SeqRegionORM.length,
                        SeqRegionORM.coord_system_id
                ))
        .join(SeqRegionORM.coord_system)
        .where(and_(CoordSystemORM.name == coord_system, SeqRegionORM.name == seq_region_name))
        )
        if version:
            stmt = stmt.filter(CoordSystemORM.version == version)
        
        rows = session.execute(stmt).all()
        
        for row in rows:
            topology = 'CIRCULAR' if cls._is_circular(session, row.SeqRegion.seq_region_id) else 'LINEAR'
            cs = CoordSystemAdaptor.fetch_by_dbID(session, row.SeqRegion.coord_system_id)
            region = Region(row.SeqRegion.sr_name,
                            topology=topology,
                            length=row.SeqRegion.length,
                            coord_system=cs)
            end = end if end else row.SeqRegion.length
            location = Location(start, end, row.SeqRegion.length)
            st = Strand(strand)
            slices.append(Slice(region, location, st))
            
        if len(slices) == 0:
            warnings.warn(f'Could not find any region with name {seq_region_name} and coordinate {coord_system}', UserWarning)

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
        cs = CoordSystemAdaptor.fetch_by_dbID(session, result.CoordSystem.coord_system_id)
        region = Region(result.SeqRegion.name,
                        topology=topology,
                        length=result.SeqRegion.length,
                        coord_system=cs)
        # location = Location(start, end, row.sr_length)
        # st = Strand.REVERSE if strand == -1 else Strand.FORWARD
        return Slice(region=region)
        
    @classmethod
    def fetch_by_name(cls,
                      session: Session,
                      name: str
                      ) -> Slice:
        """
        Arg [1]    : Session session - sqlalchemy.orm.Session object to connect to DBs
        Arg [2]    : String name - name of the Slice to look for
        Example    : name  = 'chromosome:NCBI34:X:1000000:2000000:1'
                     slice = slice_adaptor.fetch_by_name(name);
                     slice2 = slice_adaptor.fetch_by_name(slice3.name);
        Description: Fetches a slice using a slice name (i.e. the value returned by
                     the Slice.name property).  This is useful if you wish to 
                     store a unique identifier for a slice in a file or database or
                     pass a slice over a network.
                     Slice.name allows you to serialise/marshall a slice and this
                     method allows you to deserialise/unmarshal it.
                     Returns None if no seq_region with the provided name exists in
                     the database.
        Returntype : ensembl.api.core.Slice or None
        Exceptions : throw if incorrent arg provided
        Caller     : Pipeline
        Status     : At Risk
                   : under development
        """
        if not ':' in name:
            raise ArgumentError('Malfomred name for Slice')
        
        array = name.split(':')

        if len(array) < 3 or len(array) > 6:
            raise ArgumentError(
                f'Malfomred slice name [{name}]. Format is "coord_system:version:name:start:end:strand"')
        
        # make len(array) = 6, adding None as appropriate
        array.extend((None,)*(6-len(array)))

        return cls.fetch_by_seq_region(session,
                                       seq_region_name=array[2],
                                       coord_system=array[0],
                                       version=array[1],
                                       start=array[3],
                                       end=array[4],
                                       strand=array[5]
                                       )

    
    @classmethod
    def fetch_all(cls, 
                  session: Session,
                  coord_system_name: str,
                  coord_system_version: str = None,
                  include_non_reference: bool = False,
                  include_duplicates: bool = False,
                  include_lrg: bool = False
                  ) -> list[Slice]:
        """
        Arg [1]    : Session session - sqlalchemy.orm.Session object to connect to DBs
        Arg [2]    : Str coord_system_name
                     The name of the coordinate system to retrieve slices of.
                     This may be a name of an acutal coordinate system or an alias
                     to a coordinate system.  Valid aliases are 'seqlevel' or
                     'toplevel'.
        Arg [3]    : Str coord_system_version (optional - default None)
                     The version of the coordinate system to retrieve slices of
        Arg [4]    : bool include_non_reference (optional - default False)
                     If this argument is not provided then only reference slices
                     will be returned. If set, both reference and non reference
                     slices will be returned.
        Arg [5]    : bool include_duplicates (optional - default False)
                     If set duplicate regions will be returned.
                    
                     NOTE: if you do not use this option and you have a PAR
                     (pseudo-autosomal region) at the beginning of your seq_region
                     then your slice will not start at position 1, so coordinates
                     retrieved from this slice might not be what you expected.
        Arg [6]    : bool include_lrg (optional - default False)
                     If set lrg regions will be returned aswell.
        Example    : chromos = SliceAdaptor.fetch_all('chromosome','NCBI33')
                     contigs = SliceAdaptor.fetch_all('contig')
                     # get even non-reference regions
                     slices = SliceAdaptor.fetch_all('toplevel',undef,1)
                     # include duplicate regions (such as pseudo autosomal regions)
                     slices = SliceAdaptor.fetch_all('toplevel', undef,0,1)
        Description: Retrieves slices of all seq_regions for a given coordinate
                     system. Slices fetched span the entire seq_regions and are 
                     on the forward strand.
                     If the coordinate system with the provided name and version
                     does not exist an empty list is returned.
                     If the coordinate system name provided is 'toplevel', all
                     non-redundant toplevel slices are returned (note that any
                     coord_system_version argument is ignored in that case).
        Returntype : list[ensembl.api.core.Slice]
        Exceptions : none
        Caller     : general
        Status     : At Risk
                   : under development
        """
        orig_cs = CoordSystemAdaptor.fetch_by_name(session, coord_system_name, coord_system_version)
        if not orig_cs:
            return []
        
        bad_vals = {}
        if not include_non_reference:
            bad_vals.update(dict.fromkeys(cls._fetch_all_seq_region_ids_by_at_code(session, 'non_ref'), 1))
        
        if not include_lrg:
            bad_vals.update(dict.fromkeys(cls._fetch_all_seq_region_ids_by_at_code(session, 'lrg'), 1))

        if orig_cs.is_toplevel:
            seq_regs = cls._fetch_all_toplevel_seq_regions(session)
        else:
            seq_regs = cls._fetch_all_seq_regions_by_coord_system_id(session, orig_cs.internal_id)

        for (seq_region_id, name, length, cs_id) in seq_regs:
            if seq_region_id in bad_vals.keys():
                continue
            cs = CoordSystemAdaptor.fetch_by_dbID(session, cs_id)
            if not cs:
                Exception(f'seq_region {name} references non-existent coord_system {cs_id}.')

            # insert here cache mgmt - cache values for future reference

            loc = Location(1, length, length)
            tpl = RegionTopology.CIRCULAR if cls._is_circular(session, seq_region_id) else RegionTopology.LINEAR
            reg = Region(name,
                        topology=tpl,
                        length=length,
                        coord_system=cs)
            slice = Slice(reg, loc, None)
    

    @classmethod
    def _fetch_all_toplevel_seq_regions(cls, session: Session, species_id: int = 1) -> tuple:
        stmt = (
            select(
                SeqRegionORM.seq_region_id,
                SeqRegionORM.name,
                SeqRegionORM.length,
                SeqRegionORM.coord_system_id
            )
            .join(SeqRegionAttribORM.seq_region)
            .join(SeqRegionORM.coord_system)
            .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
            .where(and_(
                AttribTypeORM.code == 'toplevel',
                CoordSystemORM.species_id == species_id
            )
            )
        )
        return tuple(session.scalars(stmt).all())
    

    @classmethod
    def _fetch_all_seq_regions_by_coord_system_id(cls, session: Session, coord_system_id: int) -> tuple:
        stmt = (
            select(
                SeqRegionORM.seq_region_id,
                SeqRegionORM.name,
                SeqRegionORM.length,
                SeqRegionORM.coord_system_id
            )
            .where(SeqRegionORM.coord_system_id == coord_system_id)
        )
        return tuple(session.scalars(stmt).all())


    @classmethod
    def _fetch_all_seq_region_ids_by_at_code(cls, session: Session, at_code: str, species_id: int = 1) -> tuple[int]:
        if not at_code:
            raise ValueError('Attribute Type code must be provided')
        stmt = (
            select(
                SeqRegionORM.seq_region_id
            )
            .join(SeqRegionAttribORM.seq_region)
            .join(SeqRegionORM.coord_system)
            .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
            .where(and_(
                AttribTypeORM.code == at_code.lower(),
                CoordSystemORM.species_id == species_id
            )
            )
        )
        return tuple(session.scalars(stmt).all())
    

    @classmethod
    def _is_circular(cls, session: Session, seq_region_id: int) -> bool:
        res = cls._fetch_seq_region_attrib(session, seq_region_id, 'circular_seq')
        if res:
            return True
        return False
    
    @classmethod
    def _is_top_level(cls, session: Session, seq_region_id: int) -> bool:
        res = cls._fetch_seq_region_attrib(session, seq_region_id, 'toplevel')
        if res:
            return True
        return False

    @classmethod
    def _karyotype_rank(cls, session: Session, seq_region_id: int) -> int:
        return cls._fetch_seq_region_attrib(session, seq_region_id, 'karyotype_rank')

    @classmethod
    def _fetch_seq_region_attrib(cls, session: Session, seq_region_id: int, seq_region_attrib_code: str) -> str:
        res = (session.query(SeqRegionAttribORM.value)
        .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
        .where(AttribTypeORM.code == seq_region_attrib_code)
        .where(SeqRegionAttribORM.seq_region_id == seq_region_id)
        .first())
        if res:
            return res
        return ""
    