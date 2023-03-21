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
#   # fetch a transcript from the db and create a lightweight
#   # transcript object from it
#   my $tr = $transcript_adaptor->fetch_by_stable_id('ENST000345437');
#   my $lightweight_tr =
#     Bio::EnsEMBL::IdMapping::TinyTranscript->new_fast( [
#       $tr->dbID,          $tr->stable_id,
#       $tr->version,       $tr->created_date,
#       $tr->modified_date, $tr->start,
#       $tr->end,           $tr->strand,
#       $tr->length,        md5_hex( $tr->spliced_seq ),
#     ] );
#

# =head1 METHODS
#   start
#   end
#   strand
#   length
#   seq_md5_sum
#   add_Translation
#   translation
#   add_Exon
#   get_all_Exons
# =cut

# package Bio::EnsEMBL::IdMapping::TinyTranscript;

# # internal data structure (array indices):
# #
# #  0-4 see TinyFeature
# #  5  start
# #  6  end
# #  7  strand
# #  8  length
# #  9  seq_md5_sum
# # 10  translation
# # 11  [exons]
# # 12  biotype
# # 13  slice

from .TinyFeature import TinyFeature
from .TinyExon import TinyExon
from .TinyTranslation import TinyTranslation

from typing import Union
from dataclasses import dataclass, field

__all__ = [ 'TinyTranscript', 'TinyTranscriptData', 'DigestConversionDescriptor' ]

class DigestConversionDescriptor:
    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, obj, objtype=None):
        if obj:
            return vars(obj)[self.name]
        return None

    def __set__(self, obj, value):
        if not isinstance(value, Union[int,str]):
            raise ValueError(f'Need an int or str type.')
        if isinstance(value, str):
            value = int(value, 16)
        vars(obj)[self.name] = value

@dataclass
class TinyTranscriptData(TinyFeature):
    start: int = field(repr=False)
    end: int = field(repr=False)
    strand: int = field(repr=False)
    length: int = field(repr=False)
    seq_digest: Union[int, str] = field(repr=False, default=None)
    translation: TinyTranslation = field(repr=False, default=None)
    exons: list[TinyExon] = field(repr=False, default_factory=list)
    biotype: str = field(repr=False, default=None)
    seq_region_name: str = field(repr=False, default=None)

    
class TinyTranscript(TinyTranscriptData):
    """
    This is a lightweight transcript object for the stable Id mapping. See
    the documentation in TinyFeature for general considerations about its
    design.
    """
    seq_digest = DigestConversionDescriptor()

    def add_translation(self, translation: TinyTranslation) -> None:
        """
        Arg[1]      : ensembl.id_mapping.TinyTranslation - the translation to add
        Example     : tiny_transcript.add_Translation(tiny_translation);
        Description : Adds a translation to this transcript.
        Return type : none
        Exceptions  : thrown on wrong or missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not isinstance(translation, TinyTranslation):
            raise ValueError(f'Need an ensembl.id_mapping.TinyTranslation object.')
        self.translation(translation)

    def add_exon(self, exon: TinyExon) -> None:
        """
        Arg[1]      : ensembl.id_mapping.TinyExon exon - the exon to add
        Example     : tiny_transcript.add_Exon(tiny_exon)
        Description : Adds an exon to this transcript.
        Return type : none
        Exceptions  : thrown on wrong or missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not isinstance(exon, TinyExon):
            raise ValueError(f'Need an ensembl.id_mapping.TinyExon object.')
        self._exons.append(exon)

    def get_all_exons(self) -> list[TinyExon]:
        """
        Example     : for exon in tiny_transcript.get_all_exons():
                        # do something with exon
        Description : Returns all exons attached to that transcript.
                      It's synonym to tiny_transcript.exons
        Return type : List of ensembl.id_mapping.TinyExon objects
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self.exons