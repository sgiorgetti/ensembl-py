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

from .Assembly import CoordSystem
from .Slice import Slice
from .Strand import Strand

__all__ = [ 'AssemblyMapper' ]

class AssemblyMapper():

    def __init__(self, cs_from: CoordSystem, cs_to: CoordSystem, asm_cmp_mappings: dict = None) -> None:
        if cs_from is None or not isinstance(cs_from, CoordSystem):
            raise ValueError(f'Need an assembled coordinate system to start from')
        if cs_to is None or not isinstance(cs_to, CoordSystem):
            raise ValueError(f'Need a component coordinate system to go to')

        self._asm_cmp_mappings = asm_cmp_mappings if asm_cmp_mappings else {}
        self._asm_cs = cs_from
        self._cmp_cs = cs_to

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._asm_cs}->{self._cmp_cs})'

    #     my @coords = $asm_mapper->map($normal_slice->seq_region_name(),
	# 			  $normal_slice->start(),
	# 			  $normal_slice->end(),
	# 			  $normal_slice->strand(),
	# 			  $slice_cs, undef, undef, 1);

    def set_asm_cmp_mappings(self, raw_mappings) -> None:
        # DANGEROUS - can lead to inconsistencies!!!
        self._asm_cmp_mappings = raw_mappings

    def map(self, slice_from: Slice, include_original_coords: bool = False):
        if not slice_from or include_original_coords is None:
            raise ValueError('Wrong number of arguments')
        if not isinstance(slice_from, Slice):
            raise ValueError('You must provide slices as argument')
        if slice_from.coord_system != self._asm_cs:
            raise ValueError(f'Coordinate systems mismatch between from slice {slice_from.coord_system} and AssemblyMapper {self._asm_cs}')
        if not self._asm_cmp_mappings:
            raise Exception(f'asm_cmp_mapping cache is missing.')
        
        # Select cmp slices (gaps are still implicit)
        mapped = list(filter(lambda m: m in slice_from, self._asm_cmp_mappings.get(slice_from.seq_region_id)))

        # adjust first slice's start and end slice's end to fit slice_from
        start_offset = slice_from.start - mapped[0].asm_start
        end_offset = mapped[-1].asm_end - slice_from.end
        # Adjust first and last MappedSlice objs to match slice_from
        mapped[0].start += start_offset
        mapped[-1].end -= end_offset
        # mapped[0].asm_start += start_offset
        # mapped[-1].asm_end -= end_offset
        
        return mapped

        



        
        
        