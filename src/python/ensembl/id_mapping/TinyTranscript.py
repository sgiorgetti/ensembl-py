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

__all__ = [ 'TinyTranscript' ]

class TinyTranscript(TinyFeature):
    """
    This is a lightweight transcript object for the stable Id mapping. See
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
                 length: int,
                 seq_digest: Union[int, str],
                 translation: TinyTranslation = None,
                 exons: list[TinyExon] = None,
                 biotype: str = None,
                 seq_region_name: str = None) -> None:
        super().__init__(internal_id, stable_id, version, created_date, modified_date)
        self._start = start
        self._end = end
        self._strand = strand
        self._length = length
        self._seq_digest = seq_digest if isinstance(seq_digest, int) else int(seq_digest, 16)
        self._translation = translation
        self._exons = exons
        self._biotype = biotype
        self._seq_region_name = seq_region_name

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._internal_id}:{self._stable_id}.{self._version})'
    
    @property
    def start(self) -> int:
        """
        Description : Getter for the transcript's start coordinate.
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
        Arg[1]      : Int - the transcript's start coordinate
        Description : Setter for the transcript's start coordinate.
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
        Description : Getter for the transcript's end coordinate.
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
        Arg[1]      : Int - the transcript's end coordinate
        Description : Setter for the transcript's end coordinate.
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
        Description : Getter for the transcript's strand coordinate.
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
        Arg[1]      : Int - the transcript's strand coordinate
        Description : Setter for the transcript's strand coordinate.
        Return type : none
        Exceptions  : thrown on missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._strand = strand

    @property
    def length(self) -> int:
        """
        Description : Getter for the transcript's length. Note that this is
                      *not* the distance between start and end, but rather the sum of
                      the lengths of all exons.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._length
    
    @length.setter
    def length(self, length: int) -> None:
        """
        Arg[1]      : Int - the transcript's exonic length
        Description : Setter for the transcript's length. Note that this is
                      *not* the distance between start and end, but rather the sum of
                      the lengths of all exons.
        Return type : none
        Exceptions  : thrown on missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._length = length
    
    @property
    def seq_digest(self) -> int:
        """
        Description : Getter for the digest of the transcript's sequence.
                      Note that when used as a setter, you are expected to pass a
                      digest, not the raw sequence (i.e. the digest is not created for
                      you). Also, the digest is stored and returned as integer!
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._seq_digest
    
    @seq_digest.setter
    def length(self, seq_digest: Union[int,str]) -> None:
        """
        Arg[1]      : String or Int - the digest of the transcript's sequence
        Description : Getter for the digest of the transcript's sequence.
                      Note that when used as a setter, you are expected to pass a
                      digest, not the raw sequence (i.e. the digest is not created for
                      you). A hex string may be passed instead of Int to the setter.
        Return type : String or Int
        Exceptions  : thrown on wrong or missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not isinstance(seq_digest, Union[int,str]):
            raise ValueError(f'Need an int or str type.')
        self._seq_digest = seq_digest if isinstance(seq_digest, int) else int(seq_digest, 16)

    @property
    def translation(self) -> TinyTranslation:
        """
        Description : Getter for the transcript's translation.
        Return type : ensembl.id_mapping.TinyTranslation
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._translation
    
    @translation.setter
    def translation(self, translation: TinyTranslation) -> None:
        """
        Description : Setter for the transcript's translation.
        Return type : none
        Exceptions  : thrown on wrong or missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not isinstance(translation, TinyTranslation):
            raise ValueError(f'Need an ensembl.id_mapping.TinyTranslation object.')
        self._translation = translation
    
    @property
    def exons(self) -> list[TinyExon]:
        """
        Example     : for exon in tiny_transcrip.exons:
                        # do something with exon
        Description : Returns all exons attached to that transcript.
        Return type : List of ensembl.id_mapping.TinyExon objects
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        if not self._exons:
            return []
        return self._exons
    
    @property
    def biotype(self) -> str:
        """
        Description : Getter for the transcript's biotype.
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
        Arg[1]      : String - the transcript's biotype
        Description : Setter for the transcript's biotype.
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
        Description : Getter for the seq_region name of the slice the transcript is on.
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
        Description : Setter for the seq_region name of the slice the transcript is on.
        Return type : String
        Exceptions  : thrown on missing argument
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._seq_region_name = seq_region_name
    
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
    
    