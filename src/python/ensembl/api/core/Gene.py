from typing import List, Dict, Any
from ensembl.api.core.Metadata import Metadata
from ensembl.api.core.Transcript import Transcript

class Gene(object):

    __type = 'Gene'
    
    def __init__(self) -> None:
        pass

    def __repr__(self) -> str:
        return f"I am an Ensembl gene! I am {self._symbol}!!!"

    @property
    def alternative_symbols(self) -> List[str]:
        return self._alternative_symbols
    
    @property
    def external_id(self) -> str:
        return self._external_id

    @external_id.setter
    def external_id(self, value: str) -> None:
        self._external_id = value

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        self._name = value
    
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

    def get_slice(self) -> Any:
        return self._slice
    
    def set_slice(self, slice: Any) -> None:
        self._slice = slice
    
    def get_transcripts(self) -> List[Transcript]:
        return self._transcripts
    
    def set_transcripts(self, transcripts: List[Transcript]) -> None:
        self._transcripts = transcripts
    
    def get_metadata(self) -> Dict[str, Metadata]:
        return self._metadata
    
    def set_metadata(self, metadata_item: Dict[str, Metadata]) -> None:
        self._metadata = metadata_item


def main():
    # g = Gene('BRCA2')
    g = Gene()
    g.symbol = 'BRCA2'
    g.stable_id = 'ENSG000999888.34'
    # Gene._type = 'sss'
    print(f"Object type: {g.type} - {g.unversioned_stable_id}")
    print(f"{g}")

if __name__ == '__main__':
    main()