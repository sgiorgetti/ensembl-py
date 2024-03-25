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
"""GB Exon Cluster module"""

__all__ = []

import warnings
from ensembl.api.core import Exon, Strand, Transcript

class ExonCluster():
    def __init__(self) -> None:
        self._start: int = None
        self._end: int = None
        self._strand: Strand = None
        self._ec_type = None
        self._exonhash: dict[Exon, int] = {}
        self._exonidhash: dict[str, Exon] = {}
        self._transcripthash: dict[Transcript, list[Exon]] = {}
        self._internal_index: int = 0
        self._exon_2_biotype = None # holds biotype of exon
        self._all_exons = [] # holds all exons in the cluster

    def __repr__(self) -> str:
        pass

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

    @property
    def length(self) -> int:
        return self._end - self._start + 1

    @property
    def strand(self) -> Strand:
        return self._strand

    @property
    def ec_type(self) -> str:
        return ""

    def add_exon(self, exon: Exon, tr: Transcript, ignore_strand: bool) -> None:
        if len(self._all_exons) == 0:
            # Add new exon internal method
            pass
        #   $self->_add_transcript_reference($exon,$transcript);
        #   $self->_add_exon_biotype($exon,$transcript) ; 
        self._ec_type = None

    # INTERNAL METHODS
    def _add_new_exon (self, exon: Exon, ignore_strand: bool) -> None:
        if not self._start or exon.start < self._start:
            self._start = exon.start
        if not self._end or exon.end > self._end:
            self._end = exon.end
        if not ignore_strand:
            if self._strand is None:
                self._strand = exon.strand
            elif self._strand != exon.strand:
                raise ValueError(f"Trying to add exon with strand {exon.strand} \
                                 to cluster with strand {self._strand}")
        self._exonhash[exon] = self._internal_index + 1
        if exon in self._exonidhash:
            msg = f"There seem to be exons with the same dbID and dbname \
                in the databases $exon->dbID $exon->adaptor->db->dbname"
            raise Exception(msg)
        self._exonidhash[someExonID] = exon
        self._all_exons.append(exon)

    def _add_transcript_reference(self, exon: Exon, transcript: Transcript) -> None:
        """
        Name : _add_transcript_reference(exon, transcript) 
        Arg[1] : ensembl.api.core.Exon
        Arg[2] : ensembl.api.core.Transcript
        Function : called by add_exon , builds relation transcript-list[exon]
        """
        if not self._transcripthash.get(transcript):
            self._transcripthash[transcript] = []
        self._transcripthash[transcript] = self._transcripthash[transcript].append(exon)
        self._transcripthash[transcript] = transcript