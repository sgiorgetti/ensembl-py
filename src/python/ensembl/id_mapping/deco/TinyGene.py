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

from dataclasses import dataclass, field

__all__ = [ 'TinyGene', 'TinyGeneData' ]

@dataclass
class TinyGeneData(TinyFeature):
    start: int = field(repr=False)
    end: int = field(repr=False)
    strand: int = field(repr=False)
    seq_region_name: str = field(repr=False, default=None)
    biotype: str = field(repr=False, default=None)
    transcripts: list[TinyTranscript] = field(repr=False, default_factory=list)


class TinyGene(TinyGeneData):
    """
    This is a lightweight gene object for the stable Id mapping. See
    the documentation in TinyFeature for general considerations about its
    design.
    """
    
    def add_transcript(self, tr: TinyTranscript) -> None:
        """
        Arg[1]      : ensembl.id_mapping.TinyTranscript tr - the transcript to add
        Example     : tiny_transcript.add_transcript(tiny_transcript)
        Description : Adds an transcript to this gene.
        Return type : none
        Exceptions  : thrown on wrong or missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not isinstance(tr, TinyTranscript):
            raise ValueError(f'Need an ensembl.id_mapping.TinyTranscript object.')
        self._transcripts.append(tr)

    def get_all_transcripts(self) -> list[TinyTranscript]:
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
    