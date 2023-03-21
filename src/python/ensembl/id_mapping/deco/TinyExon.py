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


# =head1 CONTACT

#   Please email comments or questions to the public Ensembl
#   developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

#   Questions may also be sent to the Ensembl help desk at
#   <http://www.ensembl.org/Help/Contact>.

# =cut

# =head1 NAME

# Bio::EnsEMBL::IdMapping::TinyExon - lightweight exon object

# =head1 SYNOPSIS

#   # fetch an exon from the db and create a lightweight exon object
#   # from it
#   my $exon = $exon_adaptor->fetch_by_stable_id('ENSE000345437');
#   my $lightweight_exon = Bio::EnsEMBL::IdMapping::TinyExon->new_fast( [
#       $exon->dbID,
#       $exon->stable_id,
#       $exon->version,
#       $exon->created_date,
#       $exon->modified_date,
#       $exon->start,
#       $exon->end,
#       $exon->strand,
#       $exon->slice->seq_region_name,
#       $exon->slice->coord_system_name,
#       $exon->slice->coord_system->version,
#       $exon->slice->subseq( $exon->start, $exon->end, $exon->strand ),
#       $exon->phase,
#       $need_project,
#   ] );


# =head1 METHODS

#   start
#   end
#   strand
#   seq_region_name
#   coord_system_name
#   coord_system_version
#   seq
#   phase
#   need_project
#   common_start
#   common_end
#   common_strand
#   common_sr_name
#   length

# =cut


# internal data structure (array indices):
#
#  0-4 see TinyFeature
#  5  start
#  6  end
#  7  strand
#  8  seq_region_name
#  9  coord_system_name
# 10  coord_system_version
# 11  seq
# 12  phase
# 13  need_project
# 14  common_start
# 15  common_end
# 16  common_strand
# 17  common_sr_name

from .TinyFeature import TinyFeature
from dataclasses import dataclass, field

__all__ = [ 'TinyExon', 'TinyExonData' ]

@dataclass
class TinyExonData(TinyFeature):
    start: int = field(repr=False)
    end: int = field(repr=False)
    strand: int = field(repr=False)
    seq_region_name: str = field(repr=False)
    coord_system_name: str = field(repr=False)
    coord_system_version: str = field(repr=False)
    seq: str = field(repr=False, default=None)
    phase: int = field(repr=False, default=None)
    need_project: bool = field(repr=False, default=None)
    common_start: int = field(repr=False, default=None)
    common_end: int = field(repr=False, default=None)
    common_strand: int = field(repr=False, default=None)
    common_sr_name: str = field(repr=False, default=None)
    length: int = field(repr=False, init=False)

    
class TinyExon(TinyExonData):
    """
    This is a lightweight exon object for the stable Id mapping. See the
    documentation in TinyFeature for general considerations about its
    design.
    """
    @property
    def common_start(self) -> int:
        """
        Description : Getter for the exon's start in common coord_system
                      coordinates. Will return $self->start if no projection to a
                      common coord_system is required.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not self.common_start and not self.need_project:
            return self.start
        return self.common_start
    
    @property
    def common_end(self) -> int:
        """
        Description : Getter for the exon's end in common coord_system
                      coordinates. Will return $self->end if no projection to a
                      common coord_system is required.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not self.common_end and not self.need_project:
            return self.end
        return self.common_end
    
    @property
    def common_strand(self) -> int:
        """
        Description : Getter for the exon's strand in common coord_system
                      coordinates. Will return $self->strand if no projection to a
                      common coord_system is required.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not self.common_strand and not self.need_project:
            return self.strand
        return self.common_strand
    
    @property
    def common_sr_name(self) -> str:
        """
        Description : Getter for the seq_region name of the exon's slice on the
                      common coord_system coordinates. Will return
                      $self->seq_region_name if no projection to a common coord_system
                      is required.
        Return type : String
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not self.common_sr_name and not self.need_project:
            return self.seq_region_name
        return self.common_sr_name
    
    @property
    def start(self):
        return self.start
    
    @start.setter
    def start(self, value):
        if value != self.start:
            self.start = value
            self.length = None

    @property
    def end(self):
        return self.end
    
    @end.setter
    def end(self, value):
        if value != self.end:
            self.end = value
            self.length = None

    @property
    def length(self) -> int:
        """
        Description : Returns the exon length (distance between start and end).
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if self.length is None:
            self.length = self.end - self.start + 1
        return self.length
