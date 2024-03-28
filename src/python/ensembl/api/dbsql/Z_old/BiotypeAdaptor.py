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

from sqlalchemy import select
from sqlalchemy.orm import Session
from sqlalchemy.engine.row import Row

from ensembl.core.models import Biotype as BiotypeORM

from ensembl.api.core.Biotype import Biotype
import warnings

__all__ = [ 'BiotypeAdaptor' ]

class BiotypeAdaptor():
    """
    This adaptor provides a means to retrieve and store information related
    to Biotypes.  Primarily this involves the retrieval or storage of
    ensembl.api.core.Biotype objects from a database.

    See ensembl.api.core.Biotype for details of the Biotype class.
    """

    _gene_biotypes = {}
    _transcript_biotypes = {}
    
    @classmethod
    def fetch_by_name_object_type(cls, session: Session, name: str, object_type: str) -> Biotype:
        """
        Arg [1]    : Session session
                     The session to connect to the DB
        Arg [2]    : String name
                     The name of the biotype to retrieve
        Arg [3]    : String object_type
                     The object type of the biotype to retrieve (gene or transcript)
        Example    : biotype = BiotypeAdaptor.fetch_by_name_object_type('mRNA', 'gene');
        Description: Retrieves a biotype object from the database via its combined key (name, object_type).
                     If the Biotype requested does not exist in the database, a new Biotype object is
                     created with the provided name and object_type to be returned.
        Returntype : ensembl.api.core.Biotype
        Exceptions : none
        """
        if object_type == 'gene':
            if not cls._gene_biotypes:
                cls.fetch_all_by_object_type(session, object_type, flush_cache=True)
            biotype = Biotype(name=name, object_type=object_type) if not cls._gene_biotypes.get(name) else cls._gene_biotypes.get(name)
        elif object_type == 'transcript':
            if not cls._transcript_biotypes:
                cls.fetch_all_by_object_type(session, object_type, flush_cache=True)
            biotype = Biotype(name=name, object_type=object_type) if not cls._transcript_biotypes.get(name) else cls._transcript_biotypes.get(name)
        else:
            biotype = None
        
        return biotype
    

    @classmethod
    def fetch_all_by_object_type(cls, session: Session, object_type: str, flush_cache: bool = False) -> tuple[Biotype]:
        """
        Arg [1]    : Session session
                     The session to connect to the DB
        Arg [2]    : String object_type
                     The object_type of the biotypes to retrieve (gene or transcript).
        Arg [2]    : Bool flush_cache (default False)
                     If True, it replaces the cache with data from the DB
        Example    : biotypes = BiotypeAdaptor.fetch_all_by_object_type('gene')
        Description: Retrieves a list of biotype objects from the database.
        Returntype : list[ensembl.api.core.Biotype] objects or empty list
        Warning    : If empty list is to be returned
        Exceptions : none
        """
        if flush_cache:
            if object_type == 'gene':
                cls._gene_biotypes = {}
            if object_type == 'transcript':
                cls._transcript_biotypes = {}

        if object_type == 'gene' and cls._gene_biotypes:
            return tuple(cls._gene_biotypes.values())
        if object_type == 'transcript' and cls._transcript_biotypes:
            return tuple(cls._transcript_biotypes.values())
        
        stmt = (select(BiotypeORM)
               .where(BiotypeORM.object_type == object_type)
               )
        rows = session.scalars(stmt).all()

        biotypes = []
        if not rows:
            warnings.warn(f'No objects retrieved. Check if object_type "{object_type}" is correct.', UserWarning)
            return biotypes
        
        for row in rows:
            biotype = cls._biotyperow_to_biotype(row)
            biotypes.append(biotype)
            if object_type == 'gene':
                cls._gene_biotypes[biotype.name] = biotype
            if object_type == 'transcript':
                cls._transcript_biotypes[biotype.name] = biotype
        
        return tuple(biotypes)
    

    @classmethod
    def _biotyperow_to_biotype(cls, row: Row) -> Biotype:
        return Biotype(
            row.biotype_id,
            row.name,
            row.object_type,
            row.biotype_group,
            row.so_acc,
            row.so_term,
            row.description,
            row.db_type
        )