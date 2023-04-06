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

from ensembl.api.dbsql.CoordSystemAdaptor import CoordSystemAdaptor

from ensembl.api.core.Assembly import CoordSystem
from ensembl.api.core.AssemblyMapper import AssemblyMapper
from ensembl.api.core.Slice import Slice
from ensembl.api.core.Slice import MappedSlice
from ensembl.api.core.Strand import Strand

from ensembl.core.models import Assembly as AssemblyORM, SeqRegion as SeqRegionORM

from functools import lru_cache
import warnings

__all__ = [ 'AssemblyMapperAdaptor' ]

class AssemblyMapperAdaptor():

    _asm_mapper_cache = {}

    def __init__(self, session: Session, cs1: CoordSystem, cs2: CoordSystem) -> None:
        if not session or not isinstance(session, Session):
            ValueError(f'session must be a Session object')
        if not cs1 or not isinstance(cs1, CoordSystem):
            ValueError(f'cs1 argument must be a ensembl.api.core.CoordSystem.')
        if not cs2 or not isinstance(cs2, CoordSystem):
            ValueError(f'cs2 argument must be a ensembl.api.core.CoordSystem.')

        if cs1.is_toplevel:
            raise NotImplementedError()
        if cs2.is_toplevel:
            raise NotImplementedError()
        
        self._session = session
        self._asm_cs = cs1
        self._cmp_cs = cs2
    

    # @lru_cache
    def fetch_by_CoordSystems(self, load_mappings: bool = False):        
        mapping_path = CoordSystemAdaptor.get_mapping_path(self._session, self._asm_cs, self._cmp_cs)
        if not mapping_path:
            return None
        
        key = ':'.join(map(lambda cs: str(cs.internal_id) if cs else '-', mapping_path))

        asm_mapper = self._asm_mapper_cache.get(key)
        if asm_mapper:
            return asm_mapper

        if len(mapping_path) == 1:
            raise Exception(f"""Incorrect mapping path defined in meta table.\n
            0 step mapping encountered between:\n
            {self._asm_cs.name} {self._asm_cs.version} and {self._cmp_cs.name} {self._cmp_cs.version}""")
        
        if len(mapping_path) == 2:
            #1 step regular mapping
            asm_mapper = AssemblyMapper(mapping_path[0], mapping_path[1])
            if load_mappings:
                asm_mapper.set_asm_cmp_mappings(self.fetch_assembly_mappings())
            self._asm_mapper_cache[key] = asm_mapper
            return asm_mapper
        
        if len(mapping_path) == 3:
            #two step chained mapping
            
            #in multi-step mapping it is possible get requests with the
            #coordinate system ordering reversed since both mappings directions
            #cache on both orderings just in case
            #e.g.   chr <-> contig <-> clone   and   clone <-> contig <-> chr
            raise NotImplementedError()

        raise NotImplementedError()
    

    # @lru_cache
    def fetch_assembly_mappings_by_slice(self, slice_asm: Slice) -> dict:
        if slice_asm and not isinstance(slice_asm, Slice):
            ValueError(f'slice_asm argument must be a ensembl.api.core.Slice or None.')
        if slice_asm.coord_system != self._asm_cs:
            ValueError(f"slice_asm's coordinate system and cs1 must match")
        
        stmt = (
            select(
            Bundle("Assembly",
                   AssemblyORM.cmp_start,
                   AssemblyORM.cmp_end,
                   AssemblyORM.cmp_seq_region_id,
                   AssemblyORM.ori,
                   AssemblyORM.asm_start,
                   AssemblyORM.asm_end
                   ),
            Bundle("SeqRegion",
                   SeqRegionORM.name,
                   SeqRegionORM.length,
                   )
            )
            .join(SeqRegionORM, AssemblyORM.cmp_seq_region_id == SeqRegionORM.seq_region_id)
            .where(
            and_(
                AssemblyORM.asm_seq_region_id == slice_asm.internal_id,
                AssemblyORM.asm_start <= slice_asm.end,
                AssemblyORM.asm_end >= slice_asm.start,
                SeqRegionORM.coord_system_id == self._cmp_cs.internal_id
            )
            )
            .order_by(AssemblyORM.asm_start)
        )

        asm_mappings = dict(slice_asm.internal_id, [])

        for row in self._session.execute(stmt).all():
            # if (len(asm_mappings[slice_asm.internal_id]) > 0 
            #     and row.Assembly.asm_start > old_asm_end + 1):
            #     asm_mappings[slice_asm.internal_id].append(Gap(old_asm_end + 1, row.Assembly.asm_start - 1))
            asm_mappings[slice_asm.internal_id].append(
                MappedSlice(
                            row.Assembly.cmp_seq_region_id,
                            row.Assembly.cmp_start,
                            row.Assembly.cmp_end,
                            self._cmp_cs,
                            row.SeqRegionCmp.name,
                            row.SeqRegionCmp.length,
                            Strand(row.Assembly.ori),
                            row.Assembly.asm_start,
                            row.Assembly.asm_end,
                            self._asm_cs
                            )
            )
        
        return asm_mappings
    
    @lru_cache(maxsize=1)
    def fetch_assembly_mappings(self) -> dict[str,tuple]:
        SeqRegionAsm = aliased(SeqRegionORM)
        stmt = (
            select(
            Bundle("Assembly",
                   AssemblyORM.cmp_start,
                   AssemblyORM.cmp_end,
                   AssemblyORM.cmp_seq_region_id,
                   AssemblyORM.ori,
                   AssemblyORM.asm_seq_region_id,
                   AssemblyORM.asm_start,
                   AssemblyORM.asm_end
                   ),
            Bundle("SeqRegionCmp",
                   SeqRegionORM.name,
                   SeqRegionORM.length,
                   ),
            Bundle("SeqRegionAsm",
                   SeqRegionAsm.name.label('asm_name')
                   )
            )
            .join(SeqRegionORM, AssemblyORM.cmp_seq_region_id == SeqRegionORM.seq_region_id)
            .join(SeqRegionAsm, AssemblyORM.asm_seq_region_id == SeqRegionAsm.seq_region_id)
            .where(
            and_(
                SeqRegionAsm.coord_system_id == self._asm_cs.internal_id,
                SeqRegionORM.coord_system_id == self._cmp_cs.internal_id
            )
            )
            .order_by(SeqRegionAsm.seq_region_id, AssemblyORM.asm_start)
        )

        asm_mappings = {}
        for row in self._session.execute(stmt).all():
            if not asm_mappings.get(row.Assembly.asm_seq_region_id):
                asm_mappings[row.Assembly.asm_seq_region_id] = []
            asm_mappings[row.Assembly.asm_seq_region_id].append(
                MappedSlice(
                            row.Assembly.cmp_seq_region_id,
                            row.Assembly.cmp_start,
                            row.Assembly.cmp_end,
                            self._cmp_cs,
                            row.SeqRegionCmp.name,
                            row.SeqRegionCmp.length,
                            Strand(row.Assembly.ori),
                            row.Assembly.asm_start,
                            row.Assembly.asm_end,
                            self._asm_cs
                            )
            )

        return asm_mappings
