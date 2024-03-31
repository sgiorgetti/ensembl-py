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
"""Transcript module"""

__all__ = [ 'Transcript' ]

from typing import Union, Self
from . import Analysis, Exon, Location, Slice, \
    Sequence, EnsemblFeature, Strand, Translation

class Transcript(EnsemblFeature):
    """Representation of a transcript.

    TO DO!!!

    The location of a transcript is always provided 5' -> 3' on the forward strand.
    Thus the start should always be lower than the end.

    Attributes:
      exons: List of exons constituting the transcript.
      start: Start position of the transcript.
      end: End position of the transcript.
      strand: Strand of the transcript, either + or -.
      location_name: Name of the region the transcript is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the transcript.
      external_id: Name of the transcript.
      introns: List of introns constituting the transcript, should be 0 if exons size is one.
      cds_genomic_start: Genomic start position of the CDS, always be lower than cds_genomic_end.
      cds_genomic_end: Genomic end position of the CDS, always be greater than cds_genomic_start.
      cds_sequence: Translateable cDNA sequence.
      translation_sequence: Protein sequence translated from cds_sequence.
    """

    __type = 'transcript'

    def __init__(self,
                 location: Location,
                 strand: Strand,
                 reg_slice: Slice,
                 analysis: Analysis = None,
                 exons: list[Exon] = None,
                 internal_id: Union[str, int] = None,
                 sequence: Sequence = None,
                 biotype: str = None,
                 source: str = 'ensembl',
                 stable_id: str = None,
                 version: int = 1,
                 is_current: bool = True,
                 is_canonical: bool = False
                ) -> None:
        self._exons = [] if not exons else exons
        self._is_canonical = is_canonical
        self._translation = None
        self._coding_region_start: int = None
        self._coding_region_end: int = None

        super().__init__(location, strand, reg_slice, analysis, internal_id, sequence,
                         biotype, source, stable_id, version, is_current)

    @classmethod
    def fastinit(cls, start: int, end: int, length: int, analysis_name: str,
                 strand: int, reg_slice: Slice = None, internal_id: int = None,
                 biotype: str = None, sequence: str = None, source: str = 'ensembl',
                 stable_id: str = None, version: int = 1, is_current: bool = True,
                 is_canonical: bool = False) -> Self:
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
        return cls(loc, st, reg_slice, an, internal_id, seq, biotype, source, stable_id,
                   version, is_current, is_canonical)

    def __repr__(self) -> str:
        if self._attributes.get('is_canonical'):
            return f'{self.__class__.__name__}({self.stable_id} - canonical)'
        return super().__repr__()

    @property
    def is_canonical(self) -> bool:
        return self._is_canonical
    
    @is_canonical.setter
    def is_canonical(self, value: bool) -> None:
        self.is_canonical = value

    @property
    def translation(self) -> Translation:
        return self._translation

    @translation.setter
    def translation(self, value: Translation) -> None:
        if not isinstance(value, Translation):
            raise AttributeError(f'Provided value must be a Translation, instead of {type(value)}')
        self._translation = value

    def set_translation(self, tl: Translation) -> None:
        self.translation = tl

    def get_exons(self, constitutive: bool = False) -> list[Exon]:
        if constitutive:
            return [ e for e in self._exons if e.is_constitutive ]
        return self._exons

    def set_exons(self, exons: list[Exon]) -> None:
        self._exons = exons

    def get_start_exon(self) -> Exon:
        return self._exons[0]

    def get_end_exon(self) -> Exon:
        return self._exons[-1]

    def is_mane_select(self) -> bool:
        return 'MANE_Select' in self._attributes

    def is_mane(self) -> bool:
        return 'mane' in self._attributes.keys().lower()

    def coding_region_start(self) -> int:
        if self._coding_region_start:
            return self._coding_region_start
        if not self._translation:
            return None
        start = 0
        strand = self._translation.start_exon.strand
        if strand == Strand.FORWARD:
            start = self._translation.start_exon.start
            start += self._translation.start - 1
        else:
            start = self._translation.end_exon.end
            start -= self._translation.end - 1
        self._coding_region_start = start
        return start

    def coding_region_end(self) -> int:
        if self._coding_region_end:
            return self._coding_region_end
        if not self._translation:
            return None

        end = 0
        strand = self._translation.start_exon.strand
        if strand == Strand.FORWARD:
            end = self._translation.end_exon.start
            end += self._translation.end - 1
        else:
            end = self._translation.start_exon.end
            end -= self._translation.start - 1
        self._coding_region_end = end
        return end

    def get_all_translateable_exons(self) -> list[Exon]:
        if not self._translation:
            return []
        trl = self._translation
        ex_list = []
        for e in self._exons:
            adjust = 0
            if e == trl.start_exon:
                adjust = trl.start - 1
                if adjust != 0:
                    ex_list.append(e.adjust_start_end(adjust, 0))
                    continue
            if e == trl.end_exon:
                adjust = trl.end - e.length
                if adjust != 0:
                    ex_list.append(e.adjust_start_end(0, adjust))
                    break
            ex_list.append(e)
        return ex_list
