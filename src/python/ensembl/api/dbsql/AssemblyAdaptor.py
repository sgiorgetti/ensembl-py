# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from sqlalchemy import and_, select, func
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound

from ensembl.core.models import CoordSystem as CoordSystemORM, Meta as MetaORM
from ensembl.core.models import AssemblyException as AssemblyExceptionORM, SeqRegion as SeqRegionORM

from ensembl.api.dbsql.CoordSystemAdaptor import CoordSystemAdaptor

from ensembl.api.core.Assembly import Assembly

__all__ = [ 'AssemblyAdaptor' ]

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
    def fetch_by_coord_system_id(cls, session: Session, coord_system_id: int) -> Assembly:
        species_id = (session.query(CoordSystemORM.species_id)
        .where(CoordSystemORM.coord_system_id == coord_system_id)
        .first())
        return cls.fetch_by_species_id(session, species_id[0])

    @classmethod
    def fetch_by_coord_system_version(cls, session: Session, coord_system_version: str, species_id: int = 1) -> Assembly:
        cs_row = (
            session.query(CoordSystemORM.version, CoordSystemORM.attrib)
            .distinct()
            .where(and_(CoordSystemORM.version == coord_system_version,
                        CoordSystemORM.species_id == species_id))
            .first()
        )
        if not cs_row:
            raise NoResultFound(f"Could not find coord system with version {coord_system_version}")
        
        if 'default' in cs_row.attrib:
            return cls.fetch_by_species_id(session, species_id, unsafe=True)
        
        meta_row = (session.query(MetaORM)
                    .where(MetaORM.species_id == species_id)
                    .where(MetaORM.meta_key == 'species.production_name')
                    .first()
        )

        return Assembly(cs_row.version,
                        cs_row.version,
                        meta_row.meta_value,
                        is_default=False)

    @classmethod
    def fetch_by_species_id(cls, session: Session, species_id: int, unsafe: bool = False) -> Assembly:
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

        if not unsafe:
            default_cs: str = CoordSystemAdaptor.fetch_default_version(session)
            if default_cs.lower() != assembly_default.lower():
                raise Exception('Inconsistent default coord system configuration meta:{assembly_default} Vs coord_system:{default_cs.name}')
        
        return Assembly(assembly_name,
                        assembly_default,
                        species_production_name,
                        assembly_accession,
                        assembly_provider_name,
                        last_geneset_update,
                        assembly_date,
                        is_default=True)
    

    @classmethod
    def check_assembly_exceptions(cls, session: Session, species_id: int = 1) -> bool:
        if species_id < 0 or not isinstance(species_id, int):
            raise ValueError(f'Species_id must be a positive integer or None (defaults to 1).')
        stmt = (select(func.count(SeqRegionORM.seq_region_id).label('cnt'))
                .join(AssemblyExceptionORM.seq_region)
                .join(SeqRegionORM.coord_system)
                .where(CoordSystemORM.species_id == species_id)
                )
        res = session.scalars(stmt).first()
        if res > 0:
            return True
        return False