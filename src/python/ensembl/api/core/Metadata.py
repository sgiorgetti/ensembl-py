
class Metadata():

    def __init__(self, accession_id, value=None) -> None:
        if accession_id is None:
            raise ValueError('An accession_id value must be set')
        if value is not None:
            self._value = value

    def __repr__(self) -> str:
        return f"{self._accession_id}-{self._value}"
        
    @property
    def value(self) -> str:
        return self._value

    @value.setter
    def value(self, value: str) -> None:
        self._value = value

    @property
    def accession_id(self) -> str:
        return self._accession_id

    @accession_id.setter
    def accession_id(self, accession_id: str) -> None:
        self._accession_id = accession_id
