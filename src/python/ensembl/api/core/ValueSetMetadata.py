from ensembl.api.core.Metadata import Metadata
from typing import Optional

class ValueSetMetadata(Metadata):
    def __init__(self, accession_id:str, value:str=None,
                 label:Optional[str]=None, definition:Optional[str]=None, 
                 description:Optional[str]=None) -> None:
        super().__init__(accession_id, value)
        if label:
            self._label = label
        if definition:
            self._definition = definition
        if description:
            self._description = description

    @property
    def label(self) -> str:
        return self._label

    @label.setter
    def label(self, label:str) -> None:
        self._label = label

    @property
    def definition(self) -> str:
        return self._definition

    @definition.setter
    def definition(self, definition:str) -> None:
        self._definition = definition

    @property
    def description(self) -> str:
        return self._description

    @description.setter
    def description(self, description:str) -> None:
        self._description = description