from ensembl.api.core.Metadata import Metadata
from typing import Any, Optional

class XrefMetadata(Metadata):
    def __init__(self, accession_id=None, value=None, url:Optional[str]=None, source:Optional[Any]=None) -> None:
        super().__init__(accession_id, value)
        if url is not None:
            self._url = url
        if source is not None:
            self._source = source
    
    def __repr__(self) -> str:
        return f"{super().__repr__()}-{self._url}-{self._source}"

    @property
    def url(self) -> str:
        return self._url

    @url.setter
    def url(self, url: str) -> None:
        self._url = url

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, source) -> None:
        self._source = source