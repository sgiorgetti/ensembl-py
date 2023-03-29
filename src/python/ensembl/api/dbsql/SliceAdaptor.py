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
from sqlalchemy.orm import Session, Bundle, aliased

from ensembl.core.models import SeqRegion as SeqRegionORM, CoordSystem as CoordSystemORM, AttribType as AttribTypeORM
from ensembl.core.models import SeqRegionAttrib as SeqRegionAttribORM

from typing import Optional

from ensembl.api.dbsql.CoordSystemAdaptor import CoordSystemAdaptor
from ensembl.api.core.Strand import Strand
from ensembl.api.core import Slice, RegionTopology

import warnings

__all__ = [ 'SliceAdaptor' ]

class SliceAdaptor():
    """Contains all the slice related functions over Slice ORM
    """
    # _asm_exc_cache: tuple = ()
    @classmethod
    def fetch_by_seq_region(cls,
                            session: Session,
                            seq_region_name: str,
                            coord_system_name: str,
                            version: Optional[str] = None, # this is Assembly name
                            start: int = 1,
                            end: Optional[int] = None,
                            strand: int = 1
                            ) -> list[Slice]:
        """
        Arg [1]    : Session session
            the session to connect to a core DB
        Arg [2]    : str seq_region_name
            the seq_region name as specified in the core DB
        Arg [3]    : str coord_system_name
            the coordinate system name of the seq_region to retrieve
        Arg [4]    : str version
            the coordinate system version as specified in the core DB - e.g. GRCh38
        Arg [5]    : int start
            the seq_region start position
        Arg [6]    : int end
            the seq_region end position
        Arg [7]    : int strand
            the seq_region strand +1 forward, -1 reverse

        Example    : region = region_adaptor.fetch_by_seq_region( session, 'X', 'chromosome' );
        Description: Retrieves a slice on the requested region.  At a minimum the
                     name the name of the seq_region to fetch, alomg with the coordinate
                     system name, must be provided.
        Returntype : list[ensembl.api.core.Slice]
        Exceptions : throw if no or wrong argument is provided
        Caller     : general
        Status     : At Risk
                   : under development
        """
        if not session:
            raise ValueError()
        
        slices = []
        if not start:
            start = 1
        start = int(start)
        if end:
            end = int(end)
        if not seq_region_name:
            raise ValueError('seq_region_name argument is required')
        if not coord_system_name:
            raise ValueError('coord_system_name argument is required')
        
        subq = (select(SeqRegionAttribORM.seq_region_id, SeqRegionAttribORM.value)
                .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
                .where(AttribTypeORM.code == 'circular_seq')
                .subquery()
                )
        sr_attr_subq = aliased(SeqRegionAttribORM, subq, name="seq_region_attrib")

        stmt = (select(
                Bundle("SeqRegion",
                        SeqRegionORM.seq_region_id,
                        SeqRegionORM.name.label('sr_name'),
                        SeqRegionORM.length,
                        SeqRegionORM.coord_system_id
                ),
                Bundle("SRAttrCircular",
                       sr_attr_subq.value))
        .join(SeqRegionORM.coord_system)
        .join(sr_attr_subq, SeqRegionORM.seq_region_id == sr_attr_subq.seq_region_id, isouter=True)
        .where(and_(CoordSystemORM.name == coord_system_name, SeqRegionORM.name == seq_region_name))
        )
        
        if version:
            if version == 'default':
                version = CoordSystemAdaptor.fetch_default_version(session)
            stmt = stmt.filter(CoordSystemORM.version == version)
        
        rows = session.execute(stmt).all()
        
        for row in rows:
            cs = CoordSystemAdaptor.fetch_by_dbID(session, row.SeqRegion.coord_system_id)
            strand = Strand(strand) if strand else Strand.UNDEFINED
            topology = 'CIRCULAR' if cls._is_circular(session, row.SRAttrCircular.value) else 'LINEAR'
            s = Slice(
                cs,
                row.SeqRegion.sr_name,
                row.SeqRegion.length,
                start,
                end if end else row.SeqRegion.length,
                strand,
                internal_id=row.SeqRegion.seq_region_id,
                topology=topology
            )   
            slices.append(s)
            
        if len(slices) == 0:
            warnings.warn(f'Could not find any region with name {seq_region_name} and coordinate {coord_system_name}', UserWarning)

        return slices
    
    @classmethod
    def fetch_by_seq_region_id(cls, 
                               session: Session,
                               seq_region_id: int,
                               start: int = 1,
                               end: int = None,
                               strand: Strand = Strand.FORWARD
                               ) -> Slice:
        if not session:
            raise ValueError()
        
        if not isinstance(seq_region_id, int) or seq_region_id <= 0:
            raise ValueError()

        subq = (select(SeqRegionAttribORM.seq_region_id, SeqRegionAttribORM.value)
                .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
                .where(AttribTypeORM.code == 'circular_seq')
                .subquery()
                )
        sr_attr_subq = aliased(SeqRegionAttribORM, subq, name="seq_region_attrib")
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
                ),
                Bundle("SRAttrCircular",
                sr_attr_subq.value))
        .join(SeqRegionORM.coord_system)
        .join(sr_attr_subq, SeqRegionORM.seq_region_id == sr_attr_subq.seq_region_id, isouter=True)
        .where(SeqRegionORM.seq_region_id == seq_region_id)
        )

        result = session.execute(stmt).first()

        topology = 'CIRCULAR' if cls._is_circular(session, result.SRAttrCircular.value) else 'LINEAR'
        cs = CoordSystemAdaptor.fetch_by_dbID(session, result.CoordSystem.coord_system_id)
        return Slice(
            cs,
            result.SeqRegion.name,
            result.SeqRegion.length,
            start,
            end if end else result.SeqRegion.length,
            strand,
            internal_id=result.SeqRegion.seq_region_id,
            topology=topology
        )  

        
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
            raise ValueError('Malfomred name for Slice')
        
        array = name.split(':')

        if len(array) < 3 or len(array) > 6:
            raise ValueError(
                f'Malfomred slice name [{name}]. Format is "coord_system:version:name:start:end:strand"')
        
        # make len(array) = 6, adding None as appropriate
        array.extend((None,)*(6-len(array)))

        sl = cls.fetch_by_seq_region(session,
                                       seq_region_name=array[2],
                                       coord_system_name=array[0],
                                       version=array[1],
                                       start=array[3],
                                       end=array[4],
                                       strand=Strand(int(array[5]))
                                       )
        if len(sl) == 1:
            return sl[0]
        warnings.warn(f'Could find more than one region with name {name}', UserWarning)
        return None

    
    @classmethod
    def fetch_all(cls, 
                  session: Session,
                  coord_system_name: str,
                  coord_system_version: str = None,
                  include_non_reference: bool = False,
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
        Arg [5]    : bool include_lrg (optional - default False)
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

        out_slices = []
        for (seq_region_id, name, length, circular_attr_value) in seq_regs:
            if seq_region_id in bad_vals.keys():
                continue
            # This below makes little sense, as we "fetch_all" by coord system
            # Hence, the CS can only be the original one.
            # cs = CoordSystemAdaptor.fetch_by_dbID(session, cs_id)
            # if not cs:
            #     Exception(f'seq_region {name} references non-existent coord_system {cs_id}.')
            cs = orig_cs

            # insert here cache mgmt - cache values for future reference

            tpl = RegionTopology.CIRCULAR if circular_attr_value else RegionTopology.LINEAR
            slice = Slice(
                cs,
                name,
                length,
                1,
                length,
                topology = tpl
            )
            out_slices.append(slice)

        return out_slices          


    @classmethod
    def fetch_by_region_unique(cls,
                               session: Session,
                               coord_system_name: str,
                               seq_region_name: str,
                               start: int = 1,
                               end: int = None,
                               strand: Strand = Strand.FORWARD,
                               version: str = '') -> list[Slice]:
        """
        Arg [1]    : session: Session - sqlalchemy.orm.Session object to connect to DBs
        Arg [2]    : coord_system_name: str (optional)
                     The name of the coordinate system of the slice to be created
                     This may be a name of an actual coordinate system or an alias
                     to a coordinate system.  Valid aliases are 'seqlevel' or
                     'toplevel'.
        Arg [3]    : seq_region_name: str
                     The name of the sequence region that the slice will be
                     created on.
        Arg [4]    : start: int (optional, default = 1)
                     The start of the slice on the sequence region
        Arg [5]    : end: int (optional, default = seq_region length)
                     The end of the slice on the sequence region
        Arg [6]    : strand: Strand (optional, default = 1)
                     The orientation of the slice on the sequence region
        Arg [7]    : version: str (optional, default = default version)
                     The version of the coordinate system to use (e.g. NCBI33)
        Example    : slice = SliceAdaptor.fetch_by_region_unique('chromosome', 'HSCHR6_MHC_COX');
        Description: Synonym to fetch_by_seq_region. Provided for convenience only

        Returntype : list[ensembl.api.core.Slice]
        Exceptions : See SliceAdaptor.fetch_by_seq_region
        Caller     : general
        Status     : At Risk
                   : under development
        """
        return cls.fetch_by_seq_region(session,
                            seq_region_name,
                            coord_system_name,
                            version, # this is Assembly name
                            start,
                            end,
                            strand)

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
        return tuple(session.execute(stmt).all())
    

    @classmethod
    def _fetch_all_seq_regions_by_coord_system_id(cls, session: Session, coord_system_id: int) -> tuple:
        subq = (select(SeqRegionAttribORM.seq_region_id, SeqRegionAttribORM.value)
                .join(AttribTypeORM, SeqRegionAttribORM.attrib_type_id == AttribTypeORM.attrib_type_id)
                .where(AttribTypeORM.code == 'circular_seq')
                .subquery()
                )
        sr_attr_subq = aliased(SeqRegionAttribORM, subq, name="seq_region_attrib")
        stmt = (
            select(
                SeqRegionORM.seq_region_id,
                SeqRegionORM.name,
                SeqRegionORM.length,
                sr_attr_subq.value
            )
            .join(sr_attr_subq, SeqRegionORM.seq_region_id == sr_attr_subq.seq_region_id, isouter=True)
            .where(SeqRegionORM.coord_system_id == coord_system_id)
        )
        return tuple(session.execute(stmt).all())


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
    
    # @classmethod
    # def _build_exception_cache(cls, session: Session, species_id: int) -> None:
    #     # build up a cache of the entire assembly exception table
    #     # it should be small anyway
    #     stmt = (select(AssemblyExceptionORM)
    #         .join(SeqRegionORM, SeqRegionORM.seq_region_id == AssemblyExceptionORM.seq_region_id)
    #         .join(CoordSystemORM, CoordSystemORM.coord_system_id == SeqRegionORM.coord_system_id)
    #         .where(CoordSystemORM.species_id == species_id)
    #     )
    #     res = session.scalars(stmt).all()
    #     myasm = []
    #     for r in res:
    #         myasm.append(r)
    #     cls._asm_exc_cache = tuple(myasm)
