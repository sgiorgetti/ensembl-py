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

from sqlalchemy import select, and_
from sqlalchemy.orm import Session

from ensembl.core.models import Meta as MetaORM

from functools import cached_property, lru_cache, cache

import warnings

__all__ = [ 'MetaAdaptor' ]

class MetaAdaptor():

    def __init__(self, session: Session) -> None:
        if not session:
            raise ValueError(f'Missing argument session')
        self._session = session

    @property
    def session(self) -> Session:
        return self._session

    def is_multispecies(self) -> bool:
        """Multispecies DB flag. Uses ``species`` property

        Raises:
            KeyError: if ``meta`` table is not in the database.
            sqlalchemy.exc.NoResultFound: if meta key ``schema_version`` is not present.
            sqlalchemy.exc.MultipleResultsFound: if meta key ``schema_version`` returns multiple rows.

        """
        ss = self.species
        if len(ss) > 1:
            return True
        return False

    @cached_property
    def species(self) -> dict[str,int]:
        """Species list of the database, located in the ``meta`` table, and retrieved as dict

        Raises:
            KeyError: if ``meta`` table is not in the database.
            sqlalchemy.exc.NoResultFound: if meta key ``species.production_name`` is not present.
            sqlalchemy.exc.MultipleResultsFound: if meta key ``species.production_name`` returns multiple rows.

        """
        stmt = (
            select(MetaORM.species_id, MetaORM.meta_value)
            .where(MetaORM.meta_key == "species.production_name")
        )
        return { k:v for v,k in self._session.execute(stmt).all() }
    
    #  This is a test!!
    @lru_cache
    def get_species_prod_name_by_id(self, species_id: int) -> str:
        stmt = (
            select(MetaORM.meta_value)
            .where(and_(MetaORM.meta_key == "species.production_name",
                        MetaORM.species_id == species_id)
                   )
        )
        return self._session.scalars(stmt).first()
    

    def get_species_id_by_prod_name(self, species_prod_name: str) -> int:
        return self.species.get(species_prod_name)
    
    @staticmethod
    def fetch_meta_value_by_key(session: Session, meta_key: str) -> tuple[str]:
        if not session:
            raise ValueError(f'Missing argument "session"')
        if not meta_key:
            raise ValueError(f'Missing argument "meta_key"')
        stmt = select(MetaORM.meta_value).where(MetaORM.meta_key == meta_key)
        res = session.scalars(stmt).all()
        return res
