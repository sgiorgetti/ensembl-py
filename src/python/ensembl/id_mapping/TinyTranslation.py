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
#   if ( my $tl = $tr->translation ) {
#     my $lightweight_tl =
#       Bio::EnsEMBL::IdMapping::TinyTranslation->new_fast( [
#         $tl->dbID,          $tl->stable_id,
#         $tl->version,       $tl->created_date,
#         $tl->modified_date, $tr->dbID,
#         $tr->translate->seq,
#       ] );
#   }
# =head1 METHODS
#   transcript_id
#   seq
# =cut

# package Bio::EnsEMBL::IdMapping::TinyTranslation;

# # internal data structure (array indices):
# #
# #  0-4 see TinyFeature
# #  5  transcript_id
# #  6  seq

from .TinyFeature import TinyFeature

__all__ = [ 'TinyTranslation' ]

class TinyTranslation(TinyFeature):
    """
    This is a lightweight translation object for the stable Id mapping. See
    the documentation in TinyFeature for general considerations about its
    design.
    """
    def __init__(self,
                 internal_id: int,
                 stable_id: str,
                 version: int,
                 created_date: int,
                 modified_date: int,
                 transcript_id: int,
                 seq: str) -> None:
        super().__init__(internal_id, stable_id, version, created_date, modified_date)
        self._transcript_id = transcript_id
        self._seq = seq

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._internal_id}:{self._stable_id}.{self._version})'
    
    @property
    def transcript_id(self) -> int:
        """
        Description : Getter for the transcript internal Id this translation is
                      attached to.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        return self._transcript_id
    
    @transcript_id.setter
    def transcript_id(self, transcript_id: int) -> None:
        """
        Arg[1]      : Int - the transcript internal Id ("dbID")
        Description : Setter for the transcript internal Id this translation is
                      attached to.
        Return type : Int
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._transcript_id = transcript_id
    
    @property
    def seq(self) -> str:
        """
        Description : Getter for the translation's sequence.
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
        Arg[1]      : Str - the translation's sequence
        Description : Setter for the translation's sequence.
        Return type : none
        Exceptions  : none
        Caller      : general
        Status      : At Risk
                    : under development
        """
        self._seq = seq
