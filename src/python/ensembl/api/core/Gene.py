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

from typing import Union

from . import Transcript, Slice, Feature, Strand, Biotype

import warnings

__all__ = [ 'Gene' ]

class Gene(Feature):

    __type = 'gene'
    
    def __init__(self,
                 stable_id: str,
                 version: int,

                 internal_id: Union[str, int] = None,
                 slice: Slice = None,
                 start: int = None,
                 end: int = None,
                 strand: Strand = None,
                 analysis: str = None,

                 transcripts: list[Transcript] = None,
                 
                 biotype: Biotype = None,
                 source: str = 'ensembl',
                 is_current: bool = True,

                 external_name: str = None,
                 external_db: str = None,
                 external_status: str = None,
                 display_xref: str = None,
                 description: str = None,
                 created_date: str = None,
                 modified_date: str = None,
                 
                 canonical_transcript_id: int = None,
                 canonical_transcript: Transcript = None,
                ) -> None:
        if not stable_id:
            raise ValueError()
        if not slice:
            raise ValueError()
        if stable_id.find('.') > 0:
            raise ValueError()
        self._stable_id = stable_id
        self._version = version

        self._slice = slice
        self._internal_id = internal_id

        self._attributes = {}
        self._transcripts = [] if not transcripts else transcripts
        
        self._biotype = biotype
        self._source = source
        self._is_current = is_current

        self._external_name = external_name
        self._external_db = external_db
        self._external_status = external_status
        self._display_xref = display_xref
        self._description = description

        super().__init__(start, end, strand, slice, analysis, internal_id, created_date, modified_date)


    def __repr__(self) -> str:
        if self._stable_id:
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


    
    
    def get_transcripts(self) -> list[Transcript]:
        return self._transcripts
    
    def set_transcripts(self, transcripts: list[Transcript]) -> None:
        self._transcripts = transcripts

    def add_transcripts(self, transcripts) -> None:
        if not transcripts:
            warnings.warn(f"Provided transcript list is empty or None", UserWarning)
            return
        if isinstance(transcripts, Transcript):
            self._transcripts.append(transcripts)
        if (isinstance(transcripts, list) or isinstance(transcripts, tuple)) and isinstance(transcripts[0], Transcript):
            self._transcripts.extend(transcripts)

    def set_attribs(self, attribs: dict[str, str]) -> None:
        self._attributes = attribs

    def get_attribs(self, att_code: str = None) -> Union[dict, tuple]:
        if att_code is None:
            return self._attributes
        return tuple(att_code, self._attributes.get(att_code))
    
    def add_attrib(self, code: str, value: str):
        self._attributes[code] = value

    
    def get_summary(self) -> dict:
        summary = super().get_summary()
        # update here
        return summary
    

# 1       ensembl_havana  gene    65419   71585   .       +       .       ID=gene:ENSG00000186092;Name=OR4F5;biotype=protein_coding;description=olfactory receptor family 4 subfamily F member 5 [Source:HGNC Symbol%3BAcc:HGNC:14825];gene_id=ENSG00000186092;logic_name=ensembl_havana_gene_homo_sapiens;version=7
    def gff3_qualifiers(self) -> dict[str, Union[str, tuple[str]]]:

        # tags = ['basic']
        # if self.is_canonical():
        #     tags.append('Ensembl_canonical')
        # if self.is_mane_select():
        #     tags.append('MANE_Select')
        
        qualifiers = {
            'source': self._source,
            'score': ".",
            'ID': f"{self.__type}:{self._stable_id}",
            'Name': self._external_name,
            'biotype': self._biotype.name,
            'gene_id': self._stable_id,
            'version': self._version,
            'logic_name': self._analysis
        }
        # qualifiers['tag'] = tuple(tags)

        return qualifiers

    