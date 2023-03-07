

__all__ = ['Assembly']

class Assembly():
    def __init__(self,
                 id: str,
                 name: str,
                 accession_id: str,
                 accessioning_body: str,
                 species: str,
                 gb_last_geneset_update: str,
                 assembly_date: str) -> None:
        self._id = id
        self._name = name
        self._accession_id = accession_id
        self._accessioning_body = accessioning_body
        self._species = species
        self._gb_last_geneset_update = gb_last_geneset_update
        self._assembly_date = assembly_date

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._species}-{self._id}-{self._accession_id})'

    @property
    def id(self) -> str:
        return self._id
    
    @property
    def name(self) -> str:
        return self._name
    
    @property
    def accession_id(self) -> str:
        return self._accession_id
    
    @property
    def accessioning_body(self) -> str:
        return self._accessioning_body
    
    @property
    def species(self) -> str:
        return self._species
    
    @property
    def gb_last_geneset_update(self) -> str:
        return self._gb_last_geneset_update
    
    @property
    def assembly_date(self) -> str:
        return self._assembly_date
    
    



# {
#   "id": "GRCh38.p13",
#   "name": "GRCh38",
#   "accession_id": "GCA_000001405.14",
#   "accessioning_body": "EGA",
#   "organism": { ... },
#   "regions": [ ... ],
#   "default": true,
#   "tolid": null
# }