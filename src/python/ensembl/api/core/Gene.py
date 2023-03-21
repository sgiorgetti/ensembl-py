from typing import Union, Optional
from ensembl.api.core.Transcript import Transcript
from ensembl.api.core.Slice import Slice

import warnings

__all__ = [ 'Gene' ]

class Gene():

    __type = 'gene'
    
    def __init__(self,
                 stable_id: str,
                 slice: Slice,
                 transcripts: list[Transcript] = None,
                 symbol: str = None,
                 internal_id: Union[str, int] = None,
                 biotype: str = None,
                 source: str = None
                ) -> None:
        if not stable_id:
            raise ValueError()
        if not Slice:
            raise ValueError()
        if stable_id.find('.') < 0:
            raise ValueError()
        (self._unversioned_stable_id, self._version) = stable_id.split('.')
        self._symbol = symbol
        self._slice = slice
        self._metadata = {}
        self._attributes = {}
        self._transcripts = [] if not transcripts else transcripts
        self._internal_id = internal_id
        self._biotype = biotype
        self._source = source

    def __repr__(self) -> str:
        if self._unversioned_stable_id:
            return f'{self.__class__.__name__}({self.stable_id})'
        return f'{self.__class__.__name__}(internal_id{self._internal_id})'

    @property
    def name(self) -> str:
        if self.get_metadata('name'):
            return self.get_metadata('name').value
        return self.stable_id
    
    @property
    def stable_id(self) -> str:
        return f"{self._unversioned_stable_id}.{self._version}"

    @stable_id.setter
    def stable_id(self, value: str) -> None:
        (self._unversioned_stable_id, self._version) = value.split('.')
    
    @property
    def symbol(self) -> str:
        return self._symbol

    @symbol.setter
    def symbol(self, value: str) -> None:
        self._symbol = value
    
    @property
    def type(self) -> str:
        return self.__type

    @property
    def unversioned_stable_id(self) -> str:
        return self._unversioned_stable_id

    @unversioned_stable_id.setter
    def unversioned_stable_id(self, value) -> None:
        self._unversioned_stable_id = value

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
    def biotype(self) -> str:
        return self._biotype

    @biotype.setter
    def biotype(self, value: str) -> None:
        self._biotype = value

    @property
    def source(self) -> str:
        return self._source

    @source.setter
    def source(self, value: str) -> None:
        self._source = value

    def get_slice(self) -> Slice:
        return self._slice
    
    def set_slice(self, slice: Slice) -> None:
        self._slice = slice
    
    def get_transcripts(self) -> list[Transcript]:
        return self._transcripts
    
    def set_transcripts(self, transcripts: list[Transcript]) -> None:
        self._transcripts = transcripts

    def add_transcripts(self, transcripts: Union[Transcript, list[Transcript]]) -> None:
        if not transcripts:
            warnings.warn(f"Provided transcript list is empty or None", UserWarning)
            return
        if isinstance(transcripts, Transcript):
            self._transcripts.append(transcripts)
        if isinstance(transcripts[0], Transcript):
            self._transcripts.extend(transcripts)
    
    def get_metadata(self, md: Optional[str] = None) -> Union[dict, tuple]:
        if md is None:
            return self._metadata
        return tuple(md, self._metadata.get(md))
    
    def set_metadata(self, metadata: dict[str, str]) -> None:
        self._metadata = metadata

    def add_metadata(self, meta_key: str, meta_value) -> None:
        self._metadata[meta_key] = meta_value

    def set_attribs(self, attribs: dict[str, str]) -> None:
        self._attributes = attribs

    def get_attribs(self, att_code: str = None) -> Union[dict, tuple]:
        if att_code is None:
            return self._attributes
        return tuple(att_code, self._attributes.get(att_code))
    
    def add_attrib(self, code: str, value: str):
        self._attributes[code] = value

    