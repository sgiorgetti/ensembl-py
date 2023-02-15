from Location import Location
from Region import Region
from Strand import Strand

class Slice(object):

    def __init__(self, location: Location=None, region: Region=None, strand: Strand=None) -> None:
        if not location:
            raise ValueError('Location object is required to instantiate a Slice')
        if not region:
            raise ValueError('Region object is required to instantiate a Slice')
        if not strand:
            raise ValueError('Strand object is required to instantiate a Slice')
        self._location = location
        self._region = region
        self._strand = strand

    def __repr__(self) -> str:
        pass

    @property
    def location(self) -> Location:
        return self._location

    @property
    def region(self) -> Region:
        return self._region

    @property
    def strand(self) -> Strand:
        return self._strand
