from ensembl.api.core.Metadata import Metadata

__all__ = ['ValueSetMetadata']

class ValueSetMetadata(Metadata):
    def __init__(self,
                 accession_id: str,
                 value: str = '',
                 label: str = '',
                 definition: str = '', 
                 description: str = '' 
        ) -> None:
        if label is not None:
            self._label = label
        if definition is not None:
            self._definition = definition
        if description is not None:
            self._description = description
        super().__init__(accession_id, value)

    @property
    def accession_id(self) -> str:
        return super().accession_id
    
    @property
    def value(self) -> str:
        return super().value
    
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