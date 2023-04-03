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

from . import SplicedExon, Slice, Feature, Strand, Biotype
from .Translation import Translation
from typing import Union


__all__ = [ 'Transcript' ]

class Transcript(Feature):
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
                 stable_id: str,
                 version: int = 1,
                 internal_id: Union[str, int] = None,
                 slice: Slice = None,
                 start: int = None,
                 end: int = None,
                 strand: Strand = None,
                 analysis: str = None,
                 exons: list[SplicedExon] = None,
                 biotype: Biotype = None,
                 source: str = 'ensembl',
                 is_current: bool = True,
                 is_canonical: bool = False,
                 external_name: str = None,
                 external_db: str = None,
                 external_status: str = None,
                 display_xref: str = None,
                 description: str = None,
                 created_date: str = None,
                 modified_date: str = None,
                 spliced_seq: str = None,
                 canonical_translation_id: int = None
                ) -> None:
        

        if not stable_id:
            raise ValueError()
        if not Slice:
            raise ValueError()
        if stable_id.find('.') > 0:
            raise ValueError()
        self._stable_id = stable_id
        self._version = version
        
        self._slice = slice
        self._attributes = {}
        self._exons = [] if not exons else exons
        self._internal_id = internal_id
        self._biotype = biotype
        self._source = source
        self._is_current = is_current
        self._is_canonical = is_canonical
        self._spliced_seq = spliced_seq
        self._translation = None
        self._canonical_translation_id = canonical_translation_id

        self._external_name = external_name
        self._external_db = external_db
        self._external_status = external_status
        self._display_xref = display_xref
        self._description = description

        super().__init__(start, end, strand, slice, analysis, internal_id, created_date, modified_date)
        
        

    def __repr__(self) -> str:
        if self._stable_id:
            if self._attributes.get('is_canonical'):
                return f'{self.__class__.__name__}({self.stable_id} - canonical)'
            return f'{self.__class__.__name__}({self.stable_id})'
        return f'{self.__class__.__name__}(internal_id{self._internal_id})'
    
    @property
    def name(self) -> str:
        if self._external_name:
            return self._external_name
        return self.stable_id
    
    @property
    def stable_id_version(self) -> str:
        return f"{self._stable_id}.{self._version}"
    
    @property
    def stable_id(self) -> str:
        return self._stable_id

    @stable_id.setter
    def stable_id(self, value: str) -> None:
        if value.find('.') > 0:
            (self._stable_id, self._version) = value.split('.')
        else:
            self._stable_id = value
    
    @property
    def type(self) -> str:
        return self.__type

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, value) -> None:
        self._version = value

    @property
    def internal_id(self) -> int:
        return self._internal_id
    
    @property
    def dbID(self) -> int:
        return self._internal_id
    
    @property
    def biotype(self) -> Biotype:
        return self._biotype

    @biotype.setter
    def biotype(self, value: Biotype) -> None:
        self._biotype = value

    @property
    def source(self) -> str:
        return self._source

    @source.setter
    def source(self, value: str) -> None:
        self._source = value

    @property
    def is_current(self) -> bool:
        return self._is_current

    @is_current.setter
    def is_current(self, value: bool) -> None:
        self._is_current = value

    @property
    def spliced_seq(self) -> str:
        self._spliced_seq

    @spliced_seq.setter
    def spliced_seq(self, value: str) -> None:
        self._spliced_seq = value

    @property
    def external_name(self) -> str:
        return self._external_name

    @external_name.setter
    def external_name(self, value: str) -> None:
        self._external_name = value

    @property
    def external_db(self) -> str:
        return self._external_db

    @external_db.setter
    def external_db(self, value: str) -> None:
        self._external_db = value

    @property
    def external_status(self) -> str:
        return self._external_status

    @external_status.setter
    def external_status(self, value: str) -> None:
        self._external_status = value

    @property
    def display_xref(self) -> str:
        return self._display_xref

    @display_xref.setter
    def display_xref(self, value: str) -> None:
        self._display_xref = value

    @property
    def description(self) -> str:
        return self._description

    @description.setter
    def description(self, value: str) -> None:
        self._description = value

    @property
    def canonical_translation_id(self) -> int:
        return self._canonical_translation_id

    # @property
    # def length(self) -> int:
    #     raise NotImplementedError
    
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

    def get_exons(self) -> list[SplicedExon]:
        return self._exons
    
    def set_exons(self, exons: list[SplicedExon]) -> None:
        self._exons = exons

    def get_introns(self) -> tuple:
        return self._introns
    
    def set_introns(self, introns: list) -> None:
        self._introns = introns
    
    def set_attribs(self, attribs: dict[str, str]) -> None:
        self._attributes = attribs

    def get_attribs(self, att_code: str = None) -> Union[dict, tuple]:
        if att_code is None:
            return self._attributes
        return tuple(att_code, self._attributes.get(att_code))
    
    def add_attrib(self, code: str, value: str) -> None:
        self._attributes[code] = value

    def is_canonical(self) -> bool:
        return True if self._is_canonical else False
    
    def is_mane_select(self) -> bool:
        if self._attributes.get('MANE_Select'):
            return True
        return False
    
    def is_mane(self) -> bool:
        for att_code in self._attributes.keys().lower():
            if 'mane' in att_code:
                return True
        return False


    def get_summary(self) -> dict[str, str]:
        """
        Example       : transcript_summary = transcript.get_summary()
        Description   : Retrieves a textual summary of this Feature.
        Returns       : Dict[str, str]
        Status        : Alpha - Intended for internal use
        """
        summary = super().get_summary()
        summary['transcript_id'] = self._stable_id

        summary['description'] = self._description
        summary['biotype'] = self._biotype.name
        summary['source'] = self._source
        summary['symbol'] = self._external_name
        summary['logic_name'] = summary['analysis']
        if self._attributes.get('ccds_transcript'):
            summary['ccdsid'] = self._attributes.get('ccds_transcript').value
        if self._attributes.get('tsl'):
            summary['transcript_support_level'] = self._attributes.get('tsl').value
        if self._attributes.get('gencode_basic'):
            summary['tag'] = 'basic'

        return summary
    
# 1       havana  mRNA    65419   71585   .       +       .       ID=transcript:ENST00000641515;Parent=gene:ENSG00000186092;Name=OR4F5-201;biotype=protein_coding;tag=basic,Ensembl_canonical,MANE_Select;transcript_id=ENST00000641515;version=2
    def gff3_qualifiers(self) -> dict[str, Union[str, tuple[str]]]:

        tags = ['basic']
        if self.is_canonical():
            tags.append('Ensembl_canonical')
        if self.is_mane_select():
            tags.append('MANE_Select')
        
        qualifiers = {
            'source': self._source,
            'score': ".",
            'ID': f"{self.__type}:{self._stable_id}",
            'Parent': "None",
            'Name': self._external_name,
            'biotype': self._biotype.name,
            'transcript_id': self._stable_id,
            'version': self._version
        }
        qualifiers['tag'] = tuple(tags)

        return qualifiers
