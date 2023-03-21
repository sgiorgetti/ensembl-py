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

__all__ = [ 'TinyExon' ]

class TinyExon(TinyFeature):
    """
    This is a lightweight exon object for the stable Id mapping. See the
    documentation in TinyFeature for general considerations about its
    design.
    """
    def __init__(self,
                 internal_id: int,
                 stable_id: str,
                 version: int,
                 created_date: int,
                 modified_date: int,
                 start: int,
                 end: int,
                 strand: int,
                 seq_region_name: str,
                 coord_system_name: str,
                 coord_system_version: str,
                 seq: str = None,
                 phase: int = None,
                 need_project: bool = None,
                 common_start: int = None,
                 common_end: int = None,
                 common_strand: int = None,
                 common_sr_name: str = None
                ) -> None:
        super().__init__(internal_id,
                 stable_id,
                 version,
                 created_date,
                 modified_date)
        self._start = start,
        self._end = end,
        self._strand = strand,
        self._seq_region_name = seq_region_name,
        self._coord_system_name = coord_system_name,
        self._coord_system_version = coord_system_version,
        self._seq = seq,
        self._phase = phase,
        self._need_project = need_project,
        self._common_start = common_start,
        self._common_end = common_end,
        self._common_strand = common_strand,
        self._common_sr_name = common_sr_name

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._internal_id}:{self._stable_id}.{self._version})'
    
    @classmethod
    def fast_from_tuple(cls, data: tuple) -> None:
        if len(data) < 14:
            raise ValueError()
        cls.__init__(data)

    @property
    def start(self) -> int:
        """
        Description : Getter for the exon's start coordinate.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._start
    
    @property
    def end(self) -> int:
        """
        Description : Getter for the exon's end coordinate.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._end
    
    @property
    def strand(self) -> int:
        """
        Description : Getter for the exon's strand.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._strand
    
    @property
    def seq_region_name(self) -> str:
        """
        Description : Getter for the seq_region name of the slice the exon is on.
        Return type : String
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._seq_region_name
    
    @property
    def coord_system_name(self) -> str:
        """
        Description : Getter for the coord_system name of the slice the exon is on.
        Return type : String
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._coord_system_name
    
    @property
    def coord_system_version(self) -> str:
        """
        Description : Getter for the coord_system version of the slice the exon is on.
        Return type : String
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._coord_system_version
    
    @property
    def seq(self) -> str:
        """
        Description : Getter for the exon's sequence.
        Return type : String
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._seq
    
    @seq.setter
    def seq(self, seq: str) -> None:
        """
        Arg[1]      : Str - the exon's sequence
        Description : Setter for the exon's sequence.
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._seq = seq

    @property
    def phase(self) -> int:
        """
        Description : Getter for the exon's phase.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._phase
    
    @phase.setter
    def phase(self, phase: int) -> None:
        """
        Arg[1]      : Int - the exon's phase
        Description : Setter for the exon's sequence.
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._phase = phase

    @property
    def need_project(self) -> bool:
        """
        Description : Getter for the attribute determining whether an exon needs
                      to be projected onto a common coord_system. You don't need
                      to do so if the native coord_system is common to the source and
                      target assemblies, or if no common coord_system is found (the
                      Cache object has methods to determine this).
        Return type : Boolean
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._need_project
    
    @need_project.setter
    def need_project(self, need_project: bool) -> None:
        """
        Arg[1]      : Int - the exon's phase
        Description : Setter for the attribute determining whether an exon needs
                      to be projected onto a common coord_system. You don't need
                      to do so if the native coord_system is common to the source and
                      target assemblies, or if no common coord_system is found (the
                      Cache object has methods to determine this).
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._need_project = need_project
    
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
        if not self._common_start and not self._need_project:
            return self._start
        return self._common_start
    
    @common_start.setter
    def common_start(self, common_start: int) -> None:
        """
        Arg[1]      : Int - the exon's start in common coord_system
        Description : Setter for the exon's start in common coord_system
                      coordinates.
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._common_start = common_start
    
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
        if not self._common_end and not self._need_project:
            return self._end
        return self._common_end
    
    @common_end.setter
    def common_end(self, common_end: int) -> None:
        """
        Arg[1]      : Int - the exon's end in common coord_system
        Description : Setter for the exon's end in common coord_system
                      coordinates.
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._common_end = common_end
    
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
        if not self._common_strand and not self._need_project:
            return self._strand
        return self._common_strand
    
    @common_strand.setter
    def common_strand(self, common_strand: int) -> None:
        """
        Arg[1]      : Int - the exon's strand in common coord_system
        Description : Setter for the exon's strand in common coord_system
                      coordinates.
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._common_strand = common_strand
    
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
        if not self._common_sr_name and not self._need_project:
            return self._seq_region_name
        return self._common_sr_name
    
    @common_sr_name.setter
    def common_sr_name(self, common_sr_name: str) -> None:
        """
        Arg[1]      : Str - seq_region name of the exon's slice on the
                      common coord_system
        Description : Setter for the seq_region name of the exon's slice on the
                      common coord_system coordinates.
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._common_sr_name = common_sr_name
    
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
        return (self._end - self._start + 1)
