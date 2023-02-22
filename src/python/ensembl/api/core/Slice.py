from ensembl.api.core.Location import Location
from ensembl.api.core.Region import Region
from ensembl.api.core.Strand import Strand

class Slice(object):

    def __init__(self, region: Region = None, location: Location = None, strand: Strand = None) -> None:
        if not region:
            raise ValueError('Region object is required to instantiate a Slice')
        self._location = location
        self._region = region
        self._strand = strand

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._region.assembly}:{self._region.name}:{self._location.start}-{self._location.end}:{self._strand.name})'

    @property
    def location(self) -> Location:
        return self._location

    @location.setter
    def location(self, value: Location) -> None:
        self._location = value

    @property
    def region(self) -> Region:
        return self._region

    @property
    def strand(self) -> Strand:
        return self._strand

    @strand.setter
    def strand(self, value: Strand) -> None:
        self._strand = value
