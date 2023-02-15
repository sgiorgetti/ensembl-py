from typing import List, Dict

class Gene(object):

    __type = 'GeneAllele'
    __internal_identifier = 1

    def __init__(self, *args):
        if len(args) > 0:
            self._symbol = args[0]

    def __repr__(self):
        return f"I am an Ensembl gene! I am {self._symbol}!!!"
    
    @property
    def alternative_names(self) -> List[str]:
        return self._alternative_names

    @property
    def alternative_symbols(self) -> List[str]:
        return self._alternative_symbols

    @property
    def internal_identifier(self) -> int:
        return self.__internal_identifier
    
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

    def get_metadata(self):
        return self._metadata
    
    def set_metadata(self, metadata_obj) -> None:
        self._metadata = metadata_obj

    def get_slice(self):
        return self._slice

    def set_slice(self, slice_obj) -> None:
        self._slice = slice_obj

    def get_transcripts(self) -> List:
        return self._transcripts

    def set_transcripts(self, transcripts) -> None:
        self._transcripts = transcripts

    def get_xrefs(self) -> List:
        return self._xrefs

    def set_xrefs(self, xrefs) -> None:
        self._xrefs = xrefs

