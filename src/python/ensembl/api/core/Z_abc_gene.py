import abc
from typing import List

class GeneI(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def __init__(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def __repr__(self):
        raise NotImplementedError('users must define __repr__ to use this base class')
    
    @property
    def alternative_names(self) -> List[str]:
        return self._alternative_names

    @property
    def alternative_symbols(self) -> List[str]:
        return self._alternative_symbols

    @property
    @abc.abstractmethod
    def metadata(self):
        raise NotImplementedError('')

    @property
    def slice(self):
        return self._slice

    @property
    @abc.abstractmethod
    def symbol(self) -> str:
        return self._symbol

    @property
    @symbol.setter
    @abc.abstractmethod
    def symbol(self, value: str) -> None:
        raise NotImplementedError('')

    @property
    def type(self):
        return 'Gene'

class Gene(GeneI):

    def __init__(self, *kargs):
        if len(kargs) > 0:
            self._symbol = kargs[0]

    def __repr__(self):
        return f"I am an Ensembl gene! I am {self._symbol}!!!"

    def metadata(self):
        return self._metadata

    # @property
    # def symbol(self) -> str:
    #     return self._symbol

    # @symbol.setter
    # def symbol(self, value) -> None:
    #     self._symbol = value


def main():
    # g = Gene('BRCA2')
    g = Gene()
    g.symbol = 'BRCA2'
    print(f"Object type: {g.type} - {g.symbol}")
    print(f"{g}")

if __name__ == '__main__':
    main()