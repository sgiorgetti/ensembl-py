from typing import Any, Optional

class Location(object):
    """Representation of a location.

    The Location data type stores coordinates that position a Feature on a Region.

    Attributes:
      start: start coordinate, closest to the 5'-end of the Feature/Region.
      end: end coordinate, closest to the 3'-end of the Feature/Region.
      length: The number of nucleotides between start and end coordinates, inclusive.
      location_modifier: used to define whether or not the feature is a "partial model" with undefined start and/or end coordinates. - TO ADD
    Raises:
      Exception: ???
    """
 
    def __init__(self, start: int, end: int, length: int, location_modifier: Optional[Any] = None) -> None:
        self._start = start
        self._end = end
        self._length = length
        self._location_modifier = location_modifier

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._start},{self._end},{self._length})'

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end
    
    @property
    def length(self) -> int:
        return self._length
    

