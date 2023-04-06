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

from ensembl.api.core.Assembly import CoordSystem
from ensembl.api.core.Strand import Strand

from Bio.Seq import Seq, reverse_complement

import warnings
from enum import Enum
from copy import copy

__all__ = [ 'Slice', 'RegionCode', 'RegionTopology', 'MappedSlice' ]

class RegionCode(Enum):
    CHROMOSOME = 1
    PLASMID = 2
    SCAFFOLD = 3
    CONTIG = 4

class RegionTopology(Enum):
    LINEAR = 1
    CIRCULAR = 2

class Slice():

    def __init__(self,
                 coord_system: CoordSystem,
                 seq_region_name: str,
                 seq_region_length: int,
                 start: int,
                 end: int,
                 strand: Strand = Strand.FORWARD,
                 seq: Seq = None,
                 internal_id: int = None,
                 topology: RegionTopology = None,
                 adaptor = None
                 ) -> None:
        if not seq_region_name:
            raise ValueError('seq_region_name argument is required')
        if not start or not end:
            raise ValueError('both start and end arguments are required')
        if not seq_region_length:
            seq_region_length = end
        if seq_region_length <= 0:
            raise ValueError('seq_region_length must be > 0')
        if coord_system:
            if not isinstance(coord_system, CoordSystem):
                raise ValueError('coord_system ergument must be a CoordSystem object')
            if coord_system.is_toplevel:
                raise ValueError('Cannot create slice on toplevel CoordSystem')
        else:
            warnings.warn(f"Slice without coordinate system", UserWarning)

        self._seq = seq
        self._coord_system = coord_system
        self._seq_region_name = seq_region_name
        self._seq_region_length = seq_region_length
        self._start = start
        self._end = end
        self._strand = strand
        self._internal_id = internal_id
        self._adaptor = adaptor
        self._topology = topology


    def __repr__(self) -> str:
        if self._coord_system:
            return f'{self.__class__.__name__}({self._coord_system.version}:{self._coord_system.name}:{self.seq_region_name}:{self.start}-{self.end}:{self._strand.value})'
        return f'{self.__class__.__name__}(NONE:{self.seq_region_name}:{self.start}-{self.end}:{self._strand.value})'

    def __contains__(self, obj) -> bool:
        if obj == self:
            return True
        if isinstance(obj, MappedSlice):
            # only one level deep is handled
            if (obj.asm_cs == self.coord_system
                and obj.asm_start <= self.end and obj.asm_end >= self.start #and obj.strand == self.strand
                ):
                return True
        elif isinstance(obj, Slice):
            if ((obj.internal_id == self.internal_id or obj.seq_region_name == self.seq_region_name)
                and obj.coord_system == self.coord_system
                and obj.start <= self.end and obj.end >= self.start and obj.strand == self.strand
                ):
                return True
        return False

    @property
    def seq_region_name(self) -> str:
        return self._seq_region_name

    @property
    def seq_region_start(self) -> int:
        return self._start
    
    @property
    def start(self) -> int:
        return self._start
    
    @start.setter
    def start(self, val: int) -> None:
        self._start = val
    
    @property
    def seq_region_end(self) -> int:
        return self._end
    
    @property
    def end(self) -> int:
        return self._end
    
    @end.setter
    def end(self, val: int) -> None:
        self._end = val
    
    @property
    def seq_region_strand(self) -> Strand:
        return self._strand
    
    @property
    def strand(self) -> Strand:
        return self._strand
    
    @property
    def seq_region_length(self) -> int:
        return self._seq_region_length
    
    @property
    def seq_region_id(self) -> int:
        return self._internal_id
    
    @property
    def internal_id(self) -> int:
        return self._internal_id
    
    @property
    def coord_system(self) -> CoordSystem:
        return self._coord_system
    
    @property
    def source(self):
        raise NotImplementedError()
    
    @property
    def coord_system_name(self) -> str:
        return self._coord_system.name if self._coord_system else ''
    
    @property
    def coord_system_version(self) -> str:
        return self._coord_system.version if self._coord_system else ''
    
    def centrepoint(self) -> int:
        return (self._start + self._end)/2.0

    def name(self) -> str:
        # 'chromosome:NCBI36:X:1:1000000:-1'
        vals = (
            self._coord_system.name if self._coord_system else '',
            self._coord_system.version if self._coord_system else '',
            self._seq_region_name,
            str(self._start),
            str(self._end),
            str(self._strand.value)
        )
        return ':'.join(vals)
    
    def length(self) -> int:
        len = self._end - self._start + 1
        if self._start > self._end and self.is_circular():
            len += self._seq_region_length
        return len
    
    def is_reference(self):
        raise NotImplementedError()
    
    def is_toplevel(self):
        raise NotImplementedError()
    
    def is_circular(self) -> bool:
        if self._topology == RegionTopology.CIRCULAR:
            return True
        return False
    
    def get_all_genes(self):
        raise NotImplementedError()
    
    def subseq(self) -> str:
        raise NotImplementedError()
    
    def get_inverted(self):
        inv_s = copy(self)
        inv_s.strand = Strand(self.strand.value * -1)
        if self._seq:
            inv_s._seq = self._seq.reverse_complement()
        return inv_s
    

    # This can be simplified, using Pythonic constructs
    def _constrain_to_region(self) -> tuple:
        #if the slice has negative coordinates or coordinates exceeding the
        #exceeding length of the sequence region we want to shrink the slice to
        #the defined region
        if self._start > self._seq_region_length or self._end < 1:
            return ()

        right_contract = 0
        left_contract  = 0
        if self._end > self._seq_region_length:
            right_contract = self._seq_region_length - self._end
        if self._start < 1:
            left_contract = self._start - 1
        
        new_slice = self
        if left_contract or right_contract:
            if self._strand == Strand.FORWARD:
                (new_slice, tpref, fpref) = self.expand(left_contract, right_contract)
            elif self._strand == Strand.REVERSE:
                (new_slice, tpref, fpref) = self.expand(right_contract, left_contract)

        return (1-left_contract, self.length()+right_contract, new_slice)
        # return [bless [1-$left_contract, $self->length()+$right_contract,
        #                 $new_slice], "Bio::EnsEMBL::ProjectionSegment" ];


    def expand(self,
               five_prime_shift: int = 0,
               three_prime_shift: int = 0,
               force_expand: bool = False
               ) -> tuple:

        if self._seq:
            warnings.warn(f"Cannot expand a slice which has a manually attached sequence.", UserWarning)
            return None
        
        if abs(five_prime_shift) + abs(three_prime_shift) == 0:
            warnings.warn(f"5' and 3' shifts are zero. Nothing to do.", UserWarning)
            return self

        sshift = five_prime_shift if self._strand == Strand.FORWARD else three_prime_shift
        eshift = three_prime_shift if self._strand == Strand.FORWARD else five_prime_shift

        new_start = self._start - sshift
        new_end   = self._end + eshift

        # Wrap around on circular slices
        if self.is_circular():
            new_start %= self._seq_region_length
            new_end %= self._seq_region_length

        if new_start > new_end and not self.is_circular():
            if force_expand:
                # Apply max possible shift, if force_expand is set
                if sshift < 0:
                # if we are contracting the slice from the start - move the
                # start just before the end
                    new_start = new_end - 1
                    sshift    = self._start - new_start

                # if the slice still has a negative length - try to move the
                # end
                if new_start > new_end and eshift < 0:
                    new_end = new_start + 1
                    eshift  = new_end - self._end

                # return the values by which the primes were actually shifted
                tpref = eshift if self._strand == Strand.FORWARD else sshift
                fpref = sshift if self._strand == Strand.FORWARD  else eshift
            
            if new_start > new_end:
                raise Exception(f'Slice start cannot be greater than slice end')

        #fastest way to copy a slice is to do a shallow hash copy
        new_slice = copy(self)
        new_slice.start = int(new_start)
        new_slice.end   = int(new_end)

        return (new_slice, tpref, fpref)


class MappedSlice(Slice):
    def __init__(self,
                 internal_id: int,
                 start: int,
                 end: int,
                 coord_system: CoordSystem,
                 seq_region_name: str,
                 seq_region_length: int,
                 strand: Strand = Strand.FORWARD,
                 asm_start: int = 0,
                 asm_end: int = 0,
                 asm_cs: CoordSystem = None,
                 seq = None
                 ) -> None:
        if not coord_system:
            raise ValueError('coord_system is required')

        self._asm_start = asm_start
        self._asm_end = asm_end
        self._asm_cs = asm_cs

        super().__init__(coord_system,
                         seq_region_name,
                         seq_region_length,
                         start, end, strand,
                         seq,
                         internal_id,
                         RegionTopology.LINEAR,
                         None)

    def __repr__(self) -> str:
        rstr = '{}({}:{}:{}:{}-{}:{}:{}:{}:{}-{})'.format(
            self.__class__.__name__,
            self._coord_system.version,
            self._coord_system.name,
            self.seq_region_name,
            self.start,
            self.end,
            self._strand.value,
            self._asm_cs.version,
            self._asm_cs.name,
            self._asm_start,
            self._asm_end
        )
        return rstr
        
    def __contains__(self, obj) -> bool:
        return super().__contains__(obj)
    
    @property
    def asm_start(self) -> int:
        return self._asm_start
    
    @asm_start.setter
    def asm_start(self, val: int) -> None:
        self._asm_start = val
    
    @property
    def asm_end(self) -> int:
        return self._asm_end
    
    @asm_end.setter
    def asm_end(self, val: int) -> None:
        self._asm_end = val
    
    @property
    def asm_cs(self) -> CoordSystem:
        return self._asm_cs
