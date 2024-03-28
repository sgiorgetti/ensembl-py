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
"""Gene module"""

__all__ = [ 'Gene' ]

from typing import Union, Self
import warnings
from . import Analysis, Location, Slice, \
    Sequence, EnsemblFeature, Strand, Transcript

class Gene(EnsemblFeature):

    __type = 'gene'

    def __init__(self,
                 location: Location,
                 strand: Strand,
                 reg_slice: Slice,
                 analysis: Analysis = None,
                 transcripts: list[Transcript] = None,
                 internal_id: Union[str, int] = None,
                 sequence: Sequence = None,
                 biotype: str = None,
                 source: str = 'ensembl',
                 stable_id: str = None,
                 version: int = 1,
                 is_current: bool = True) -> None:
                #  canonical_transcript_id: int = None,
                #  canonical_transcript: Transcript = None
        self._transcripts = [] if not transcripts else transcripts

        super().__init__(location, strand, reg_slice, analysis, internal_id, sequence,
                         biotype, source, stable_id, version, is_current)

    @classmethod
    def fastinit(cls, start: int, end: int, length: int, analysis_name: str,
                 strand: int, reg_slice: Slice = None, internal_id: int = None,
                 biotype: str = None, sequence: str = None, source: str = 'ensembl',
                 stable_id: str = None, version: int = 1, is_current: bool = True) -> Self:
        if not start:
            raise ValueError("Feature start must be specified")
        start = int(start)
        if not end and not length:
            raise ValueError("Either feature end or length must be specified")
        if not end:
            end = start + int(length) - 1
        end = int(end)
        an = Analysis(analysis_name) if analysis_name else None
        loc = Location(start, end)
        st = Strand(strand)
        seq = Sequence(seq_id=None, seq=sequence) if sequence else None
        return cls(loc, st, reg_slice, an, internal_id, seq, biotype, source,
                   stable_id, version, is_current)

    def get_transcripts(self) -> list[Transcript]:
        return self._transcripts

    def set_transcripts(self, transcripts: list[Transcript]) -> None:
        self._transcripts = transcripts

    def add_transcripts(self, transcripts) -> None:
        if not transcripts:
            warnings.warn("Provided transcript list is empty or None", UserWarning)
            return
        if isinstance(transcripts, Transcript):
            self._transcripts.append(transcripts)
        if isinstance(transcripts, (list, tuple)) and isinstance(transcripts[0], Transcript):
            self._transcripts.extend(transcripts)
