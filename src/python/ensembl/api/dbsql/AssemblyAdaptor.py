from sqlalchemy import and_
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound

from ensembl.core.models import CoordSystem as CoordSystemORM
from ensembl.core.models import SeqRegion as SeqRegionORM
from ensembl.core.models import Meta as MetaORM

from ensembl.api.core.Assembly import Assembly


__all__ = ['AssemblyAdaptor']

class AssemblyAdaptor():

    __assembly_meta_keys = (
        'assembly.provider_name',
        'assembly.provider_url',
        'assembly.name',
        'assembly.default',
        'assembly.date',
        'assembly.accession',
        'genebuild.last_geneset_update',
        'species.production_name'
    )

    @classmethod
    def fetch_by_seq_region_id(cls, session: Session, seq_region_id: int):
        species_id = (session.query(CoordSystemORM.species_id)
        .join(SeqRegionORM.coord_system)
        .where(SeqRegionORM.seq_region_id == seq_region_id)
        .first())
        return cls.fetch_by_species_id(session, species_id[0])

    @classmethod
    def fetch_by_coord_system_name(cls, session: Session, cs_name: str):
        species_id = (session.query(CoordSystemORM.species_id)
        .where(and_(CoordSystemORM.coord_system_id == cs_name, CoordSystemORM.rank == 1))
        .first())
        return cls.fetch_by_species_id(session, species_id[0])

    @classmethod
    def fetch_by_coord_system_id(cls, session: Session, coord_system_id: int) -> Assembly:
        species_id = (session.query(CoordSystemORM.species_id)
        .where(CoordSystemORM.coord_system_id == coord_system_id)
        .first())
        return cls.fetch_by_species_id(session, species_id[0])

    @classmethod
    def fetch_by_species_production_name(cls, session: Session, species_production_name: str) -> Assembly:
        species_id = (session.query(MetaORM.species_id)
        .where(MetaORM.meta_key == 'species.production_name')
        .first())
        if not species_id:
            raise NoResultFound(f"Could not find 'species.production_name' in meta table for species {species_production_name}")
        return cls.fetch_by_species_id(session, species_id[0])

    @classmethod
    def fetch_by_species_id(cls, session: Session, species_id: int) -> Assembly:
        res = (session.query(MetaORM)
        .where(MetaORM.species_id == species_id)
        .filter(MetaORM.meta_key.in_(cls.__assembly_meta_keys))
        .order_by(MetaORM.meta_id)
        .all()
        )
        if not res.count:
            raise NoResultFound(f"Could not find meta keys in meta table for species id {species_id}")
        
        (
            last_geneset_update,
            assembly_date,
            assembly_name,
            assembly_default,
            assembly_accession,
            assembly_provider_name,
            assembly_provider_url,
            species_production_name
        ) = tuple( r.meta_value for r in res )
        
        return Assembly(assembly_name,
                        assembly_default,
                        assembly_accession,
                        assembly_provider_name,
                        species_production_name,
                        last_geneset_update,
                        assembly_date)