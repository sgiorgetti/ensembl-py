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

# =head1 SYNOPSIS
#  # fetch a gene from the db and create a lightweight gene object from it
#   my $gene = $gene_adaptor->fetch_by_stable_id('ENSG000345437');
#   my $lightweight_gene = Bio::EnsEMBL::IdMapping::TinyGene->new_fast( [
#       $gene->dbID,                   $gene->stable_id,
#       $gene->version,                $gene->created_date,
#       $gene->modified_date,          $gene->start,
#       $gene->end,                    $gene->strand,
#       $gene->slice->seq_region_name, $gene->biotype,
#       $gene->analysis->logic_name,
#   ] );
# =head1 DESCRIPTION
# This is a lightweight gene object for the stable Id mapping. See the
# documentation in TinyFeature for general considerations about its
# design.
# =head1 METHODS
#   start
#   end
#   strand
#   seq_region_name
#   biotype
#   logic_name
#   add_Transcript
#   get_all_Transcripts
#   length
# =cut

# package Bio::EnsEMBL::IdMapping::TinyGene;

# # internal data structure (array indices):
# #
# #  0-4 see TinyFeature
# #  5  start
# #  6  end
# #  7  strand
# #  8  seq_region_name
# #  9  biotype
# # 10  logic_name
# # 11  [transcripts]

from .TinyFeature import TinyFeature
from .TinyTranscript import TinyTranscript
from typing import Union
import warnings

__all__ = [ 'TinyGene' ]

class TinyGene(TinyFeature):
    """
    This is a lightweight gene object for the stable Id mapping. See
    the documentation in TinyFeature for general considerations about its
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
                 seq_region_name: str = None,
                 biotype: str = None,
                 transcripts: list(TinyTranscript) = None) -> None:
        super().__init__(internal_id, stable_id, version, created_date, modified_date)
        self._start = start
        self._end = end
        self._strand = strand
        self._transcripts = transcripts
        self._biotype = biotype
        self._seq_region_name = seq_region_name

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._internal_id}:{self._stable_id}.{self._version})'
    
    @property
    def start(self) -> int:
        """
        Description : Getter for the gene's start coordinate.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._start
    
    @start.setter
    def start(self, start: int) -> None:
        """
        Arg[1]      : Int - the gene's start coordinate
        Description : Setter for the gene's start coordinate.
        Return type : none
        Exceptions  : thrown on missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._start = start
    
    @property
    def end(self) -> int:
        """
        Description : Getter for the gene's end coordinate.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._end
    
    @end.setter
    def end(self, end: int) -> None:
        """
        Arg[1]      : Int - the gene's end coordinate
        Description : Setter for the gene's end coordinate.
        Return type : none
        Exceptions  : thrown on missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._end = end
    
    @property
    def strand(self) -> int:
        """
        Description : Getter for the gene's strand coordinate.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._strand
    
    @strand.setter
    def strand(self, strand: int) -> None:
        """
        Arg[1]      : Int - the gene's strand coordinate
        Description : Setter for the gene's strand coordinate.
        Return type : none
        Exceptions  : thrown on missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._strand = strand

    @property
    def biotype(self) -> str:
        """
        Description : Getter for the gene's biotype.
        Return type : String
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._biotype
    
    @biotype.setter
    def biotype(self, biotype: str) -> None:
        """
        Arg[1]      : String - the gene's biotype
        Description : Setter for the gene's biotype.
        Return type : String
        Exceptions  : thrown on missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._biotype = biotype

    @property
    def seq_region_name(self) -> str:
        """
        Description : Getter for the seq_region name of the slice the gene is on.
        Return type : String
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._seq_region_name
    
    @seq_region_name.setter
    def seq_region_name(self, seq_region_name: str) -> None:
        """
        Arg[1]      : String - seq_region name
        Description : Setter for the seq_region name of the slice the gene is on.
        Return type : String
        Exceptions  : thrown on missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._seq_region_name = seq_region_name

    @property
    def transcripts(self) -> list(TinyTranscript):
        """
        Example     : for tr in tiny_gene.transcripts:
                        # do something with transcript
        Description : Returns all transcripts attached to that gene.
        Return type : List of ensembl.id_mapping.TinyTranscript objects
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not self._transcripts:
            return []
        return self._transcripts

    def add_transcripts(self, transcripts: Union[TinyTranscript, list[TinyTranscript]]) -> None:
        """
        Arg[1]      : ensembl.id_mapping.TinyTranscript tr - the transcript to add
        Example     : tiny_transcript.add_transcript(tiny_transcript)
        Description : Adds an transcript(s) to this gene.
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not transcripts:
            warnings.warn(f"Provided transcript list is empty or None", UserWarning)
            return
        if isinstance(transcripts, TinyTranscript):
            self._transcripts.append(transcripts)
        if isinstance(transcripts[0], TinyTranscript):
            self._transcripts.extend(transcripts)

    def get_all_transcripts(self) -> list(TinyTranscript):
        """
        Example     : for tr in tiny_gene.get_all_transcripts():
                        # do something with transcript
        Description : Returns all transcripts attached to that gene.
                      It's synonym to tiny_gene.transcripts
        Return type : List of ensembl.id_mapping.TinyTranscript objects
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self.transcripts

    def _build_gene(self):
        # Check the integrity of the exons
        strand = self._transcripts[0].slice.strand.value
        location_name = self._transcripts[0].slice.region.name
        for transcript in self._transcripts:
            if transcript.strand != strand:
                raise Exception("Inconsistent strands on the transcripts. Transcripts should all reside on same strand")
            if transcript.location_name != location_name:
                raise Exception("Inconsistent location names for the transcripts. Transcripts should belong to the same parent sequence")

        # This is not needed really, but might be useful when clustering and doing thing like that
        if strand == '+':
            self._transcripts.sort(key=lambda x: x.start)
            self.start = self._transcripts[0].slice.location.start
            self.end = self._transcripts[-1].slice.location.end
        else:
            self._transcripts.sort(key=lambda x: x.end, reverse=True)
            self.start = self._transcripts[-1].slice.location.start
            self.end = self._transcripts[0].slice.location.end

        if self.start >= self.end:
            raise Exception("Gene start was >= end, this should not be")


    