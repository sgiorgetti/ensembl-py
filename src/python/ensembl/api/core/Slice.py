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

import warnings
from enum import Enum
from copy import copy

__all__ = [ 'Slice', 'RegionCode', 'RegionTopology' ]

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
                 seq = None,
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

    @property
    def seq_region_name(self) -> str:
        return self._seq_region_name

    @property
    def seq_region_start(self) -> int:
        return self._start
    
    @property
    def start(self) -> int:
        return self._start
    
    @property
    def seq_region_end(self) -> int:
        return self._end
    
    @property
    def end(self) -> int:
        return self._end
    
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
        return 'NNN'
    

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
            warnings.warn(f"5' and 3' shifts are zer. Nothing to do.", UserWarning)
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
    
    




        

# class OldSlice():

#     def __init__(self, region: Region = None, location: Location = None, strand: Strand = None) -> None:
#         if not region:
#             raise ValueError('Region object is required to instantiate a Slice')
#         self._location = location
#         self._region = region
#         self._strand = strand

#     def __repr__(self) -> str:
#         return f'{self.__class__.__name__}({self._region._coord_system._version}:{self._region.name}:{self._location.start}-{self._location.end}:{self._strand.name})'

#     @property
#     def location(self) -> Location:
#         return self._location

#     @location.setter
#     def location(self, value: Location) -> None:
#         self._location = value

#     @property
#     def region(self) -> Region:
#         return self._region

#     @property
#     def strand(self) -> Strand:
#         return self._strand

#     @strand.setter
#     def strand(self, value: Strand) -> None:
#         self._strand = value

#     @property
#     def name(self) -> str:
#         # 'chromosome:NCBI36:X:1:1000000:-1'
#         vals = (
#             self._region.code.name.lower(),
#             self._region.coord_system.version,
#             self._region.name,
#             self._location.start,
#             self._location.end,
#             self._strand.value
#         )
#         return ':'.join(vals)
    
#     def get_all_genes(self):
#         # Arg [1]    : (optional) string $logic_name
#         #             The name of the analysis used to generate the genes to retrieve
#         # Arg [2]    : (optional) string $dbtype
#         #             The dbtype of genes to obtain.  This assumes that the db has
#         #             been added to the DBAdaptor under this name (using the
#         #             DBConnection::add_db_adaptor method).
#         # Arg [3]    : (optional) boolean $load_transcripts
#         #             If set to true, transcripts will be loaded immediately rather
#         #             than being lazy-loaded on request.  This will result in a
#         #             significant speed up if the Transcripts and Exons are going to
#         #             be used (but a slow down if they are not).
#         # Arg [4]    : (optional) string $source
#         #             The source of the genes to retrieve.
#         # Arg [5]    : (optional) string $biotype
#         #             The biotype of the genes to retrieve.
#         # Example    : @genes = @{$slice->get_all_Genes};
#         # Description: Retrieves all genes that overlap this slice, including those on
#         #             the reverse strand.
#         # Returntype : listref of Bio::EnsEMBL::Genes
#         # Exceptions : none
#         # Caller     : none
#         # Status     : Stable
#         # my ($self, $logic_name, $dbtype, $load_transcripts, $source, $biotype) = @_;
#         # if(my $adaptor = $self->_get_CoreAdaptor('Gene', $dbtype)) {
#         #     return $adaptor->fetch_all_by_Slice( $self, $logic_name, $load_transcripts, $source, $biotype);
#         # }
#         # return [];
#         with self._get_coreadaptor_conn('Gene').session_scope() as session:
#             pass

#     #     SliceAdaptor.fetch_by_name(session, self.name)
        

#     def _get_coreadaptor_conn(self, obj_type: str, dbtype: str = 'core') -> DBConnection:
#         """
#         Arg  [1]    : Str object_type to retrieve an adaptor for
#         Arg  [2]    : Str dbtype to search for the given adaptor in. Defaults to core
#         Description : Searches for the specified adaptor in the Registry and returns it. Otherwise
#                       it will return nothing if the adaptor was not found
#         ReturnType  : Bio::EnsEMBL::DBSQL::BaseAdaptor derived instance (specific to core-like dbs)
#         Exceptions  : missing object_type
#         """
#         if not obj_type:
#             raise ValueError('Object type is a required parameter')
#         if not dbtype:
#             dbtype = 'core'
#         return self._get_Adaptor(obj_type, dbtype)

#     def _get_Adaptor(self, obj_type: str, dbtype: str, check_db: bool = False) -> DBConnection:
#         """
#         Arg  [1]    : Str object_type to retrieve an adaptor for
#         Arg  [2]    : Str dbtype to search for the given adaptor in
#         Arg  [3]    : Boolean Turn off the checking of Registry->get_db() for your 
#                       adaptor.
#         Description : Searches for the specified adaptor in the Registry and returns it. Otherwise
#                       it will return nothing if the adaptor was not found. We consult the 
#                       "special" adaptors held by Bio::EnsEMBL::Registry::get_db() method and then
#                       fall back to the normal methods of finding an adaptor.
#                       This method will warn when adaptors are missing but will never through an
#                       exception. It is up to the calling code to decide how to handle the unavailablity
#                       of an adaptor.
#         ReturnType  : Bio::EnsEMBL::DBSQL::BaseAdaptor derrived instance. Otherwise it returns nothing
#         Exceptions  : none
#         """
#         ######Example:   Gene      core       ??
#         # my ($self, $object_type, $dbtype, $do_not_check_db) = @_;

#         if not obj_type:
#             raise ValueError('Object type is a required parameter')

#         species = self._region.coord_system.species
#         assembly_name = self._region._coord_system.version

#         # #First we query for the DBAdaptor using get_db(). This is a deprecated method
#         # #call but "special" adaptors can be registered via this method. We must
#         # #consult here 1st to find the possible special adaptor
#         # if(!$do_not_check_db && $dbtype) {
#         #     my $db = $registry->get_db($local_db, $dbtype);
#         #     if($db) {
#         #     # If we got a return then use this DBAdaptor's species name, group and the given object type.
#         #     # Special adaptors can have different species names
#         #     $adaptor = $registry->get_adaptor($db->species(), $db->group(), $object_type);
#         #     }
#         #     else {
#         #     #Otherwise just use the current species, dbtype and object type
#         #     $adaptor = $registry->get_adaptor($species, $dbtype, $object_type);
#         #     }
#         # }
#         # # Otherwise our query group is the one attached to the current adaptor
#         # else {
#         #     #If not set use the group attached to the local adaptor 
#         #     $dbtype ||= $local_db->group();
#         #     $adaptor = $registry->get_adaptor($species, $dbtype, $object_type);
#         # }
#         # return $adaptor if $adaptor;

#         #A string like 'anonymous@ensembldb.ensembl.org:3306/homo_sapiens_core_109_38'
#         #is returned by the (new)Registry

#         return DBConnection('mysql://ensro@mysql-ens-sta-1.ebi.ac.uk:4519/homo_sapiens_core_110_38')
        
        
