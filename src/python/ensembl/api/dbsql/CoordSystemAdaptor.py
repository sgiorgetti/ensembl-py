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

from sqlalchemy import select, func, and_
from sqlalchemy.orm import Session

from ensembl.core.models import CoordSystem as CoordSystemORM, Meta as MetaORM

from typing import Optional

from ensembl.api.core.Assembly import CoordSystem

import warnings

__all__ = [ 'CoordSystemAdaptor' ]


class CoordSystemAdaptor():
    """Contains all the coordinate system related functions over CoordSystem ORM
    This adaptor allows the querying of information on the coordinate
    system.
    Note that many coordinate systems do not have a concept of a version
    for the entire coordinate system (though they may have a per-sequence
    version).  The 'chromosome' coordinate system usually has a version
    (i.e. the assembly version) but the clonal coordinate system does not
    (despite having individual sequence versions).  In the case where a
    coordinate system does not have a version an empty string ('') is used
    instead.
    """
    
    @classmethod
    def fetch_all(cls, session: Session, species_id: int = None) -> list[CoordSystem]:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Example    : for cs in CoordSystemAdaptor.fetch_all(dbconnection):
                     print(f"{cs.name} {cs.version}";
        Description: Retrieves every coordinate system defined in the DB.
                     These will be returned in ascending order of species_id and rank. I.e.
                     The coordinate system with lower rank would be first in the
                     array.
        Returntype : List[ensembl.dbsql.CoordSystem]
        Exceptions : ensembl.api.dbsql.DBAdaptor.ArgumentError
        Caller     : general
        Status     : At Risk
                   : under development
        """

        if not species_id:
            rows = (
                session.query(CoordSystemORM, MetaORM.meta_value.label('species_prod_name'))
                .join(MetaORM, CoordSystemORM.species_id == MetaORM.species_id)
                .where(MetaORM.meta_key == 'species.production_name')
                .order_by(CoordSystemORM.species_id, CoordSystemORM.rank)
                .all()
            )
        else:
            rows = (
                session.query(CoordSystemORM, MetaORM.meta_value.label('species_prod_name'))
                .join(MetaORM, CoordSystemORM.species_id == MetaORM.species_id)
                .where(MetaORM.meta_key == 'species.production_name')
                .where(CoordSystemORM.species_id == species_id)
                .order_by(CoordSystemORM.rank)
                .all()
            )

        cs_list = []
        for row in rows:
            toplevel = True if row.name == 'top_level' else False
            seqlevel = False
            default = False
            if row.attrib:
                if 'sequence_level' in row.attrib:
                    seqlevel = True
                if 'default' in row.attrib:
                    default = True

            cs = CoordSystem(row.name, row.version, row.rank, toplevel, seqlevel, default, row.species_id, row.coord_system_id)
            cs_list.append(cs)

        if len(cs_list) <= 0:
            warnings.warn(f'Could not find any coordinate system', UserWarning)
        
        return cs_list
    

    @classmethod
    def fetch_by_rank(cls, session: Session, rank: int) -> list[CoordSystem]:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Arg [2]    : int rank
        Example    : cs_list = CoordSystemAdaptor.fetch_by_rank(session, 1)
        Description: Retrieves a CoordinateSystem via its rank. 0 is a special
                     rank reserved for the pseudo coordinate system 'toplevel'.
                     undef is returned if no coordinate system of the specified rank
                     exists.
        Returntype : List[ensembl.dbsql.CoordSystem]
        Exceptions : ensembl.api.dbsql.DBAdaptor.ArgumentError
        Caller     : general
        Status     : At Risk
                   : under development
        """
        if rank < 0:
            raise AttributeError('Rank argument must be a non-negative integer.')
        
        if rank == 0:
            pass # return fetch_top_level

        rows = (
            session.query(CoordSystemORM)
            .where(CoordSystemORM.rank == rank)
            .order_by(CoordSystemORM.species_id, CoordSystemORM.rank)
            .all()
        )

        cs_list = []

        for row in rows:
            toplevel = True if row.name == 'top_level' else False
            seqlevel = False
            default = False
            if row.attrib:
                if 'sequence_level' in row.attrib:
                    seqlevel = True
                if 'default' in row.attrib:
                    default = True

            cs = CoordSystem(row.name, row.version, row.rank, toplevel, seqlevel, default, row.species_id, row.coord_system_id)
            cs_list.append(cs)

        if len(cs_list) <= 0:
            warnings.warn(f'Could not find any coordinate system with rank {rank}', UserWarning)
        
        return cs_list

    @classmethod
    def fetch_by_name(cls, 
                      session: Session,
                      name: str, 
                      version: Optional[str] = None, 
                      species_id: int = 1
                    ) -> CoordSystem:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Arg [2]    : str name
                     The name of the coordinate system to retrieve.  Alternatively
                     this may be an alias for a real coordinate system.  Valid
                     aliases are 'toplevel' and 'seqlevel'.
        Arg [3]    : str version
                     The version of the coordinate system to retrieve.  If not
                     specified the default version will be used.
        Arg [4]    : int species_id (default: 1)
                     The species_id the coordinate system refers to.
                     If not specified the default value 1 will be used.
        Example    : cs_list = CoordSystemAdaptor.fetch_by_name('contig')
                     cs_list = CoordSystemAdaptor.fetch_by_name('chromosome','GRCh37')
        Description: Retrieves a coordinate system by its name
        Returntype : ensembl.dbsql.CoordSystem
        Exceptions : throw for wrong or missing arguments
        Caller     : general
        Status     : At Risk
                   : under development
        """
        warn_str = f'Could not find any coordinate system with name {name}'

        if name == 'seqlevel':
            # fetch sequence level 
            return cls.fetch_sequence_level(session)

        if name == 'toplevel':
            # fetch top level 
            return cls.fetch_top_level(session)
        
        if not version:
            stmt = (select(CoordSystemORM)
                .where(
                    and_(
                        func.lower(CoordSystemORM.name) == name.lower(),
                        CoordSystemORM.attrib.like(r'%default%'),
                        CoordSystemORM.species_id == species_id
                    )
                )
                .order_by(CoordSystemORM.species_id, CoordSystemORM.rank)
            )
        else:
            stmt = (select(CoordSystemORM)
            .join(MetaORM, CoordSystemORM.species_id == MetaORM.species_id)
            .where(
                and_(
                    func.lower(CoordSystemORM.name) == name.lower(),
                    func.lower(CoordSystemORM.version) == version.lower()
                )
            )
            .order_by(CoordSystemORM.species_id, CoordSystemORM.rank)
            )
            warn_str += f" and version {version}"

        cs_row = session.scalars(stmt).first()

        if not cs_row:
            warnings.warn(warn_str, UserWarning)
            return None
        
        toplevel = True if cs_row.name == 'top_level' else False
        seqlevel = False
        default = False
        if cs_row.attrib:
            if 'sequence_level' in cs_row.attrib:
                seqlevel = True
            if 'default' in cs_row.attrib:
                default = True

        return CoordSystem(cs_row.name, cs_row.version, cs_row.rank, toplevel, seqlevel, default, cs_row.species_id, cs_row.coord_system_id)

    
    @classmethod
    def fetch_sequence_level(cls, session: Session) -> CoordSystem:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Example    : cs = CoordSystemAdaptor.fetch_sequence_level(session);
        Description: Retrieves the coordinate system at which sequence
                     is stored at.
        Returntype : ensembl.dbsql.CoordSystem
        Exceptions : throw if no sequence_level coord system exists at all
                     throw if multiple sequence_level coord systems exists
        Caller     : general
        Status     : At Risk
                   : under development
        """
        rows = (
            session.query(CoordSystemORM)
                .where(CoordSystemORM.attrib.like(r'%sequence_level%'))
                .all()
        )

        if not rows:
            raise Exception(f'No sequence_level coord_system is defined')
        
        if len(rows) > 1:
            raise Exception(f'Multiple sequence_level coord_systems are defined. Only one is currently supported')
    
        toplevel = True if rows[0].name == 'top_level' else False
        seqlevel = False
        default = False
        if rows[0].attrib:
            if 'sequence_level' in rows[0].attrib:
                seqlevel = True
            if 'default' in rows[0].attrib:
                default = True

        return CoordSystem(rows[0].name, rows[0].version, rows[0].rank, toplevel, seqlevel, default, rows[0].species_id, rows[0].coord_system_id)
        

    @classmethod
    def fetch_by_dbID(cls, session: Session, dbID: int) -> CoordSystem:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Arg [2]    : dbID: Int
                     The DB ID for the coordinate system to retrieve
        Example    : cs = CoordSystemAdaptor.fetch_by_dbID(session, 4);
        Description: Retrieves the coordinate system given its internal core DB ID.
        Returntype : ensembl.dbsql.CoordSystem
        Exceptions : none
        Caller     : general
        Status     : At Risk
                   : under development
        """
        stmt = select(CoordSystemORM).where(CoordSystemORM.coord_system_id == dbID)
        cs_row = session.scalars(stmt).first()

        if not cs_row:
            return None
    
        toplevel = True if cs_row.name == 'top_level' else False
        seqlevel = False
        default = False
        if cs_row.attrib:
            if 'sequence_level' in cs_row.attrib:
                seqlevel = True
            if 'default' in cs_row.attrib:
                default = True

        return CoordSystem(cs_row.name, cs_row.version, cs_row.rank, toplevel, seqlevel, default, cs_row.species_id, cs_row.coord_system_id)


    @classmethod
    def fetch_default_version(cls, session: Session) -> str:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Example    : cs = CoordSystemAdaptor.fetch_default_version(session);
        Description: Retrieves the default version of the assembly
        Returntype : str
        Exceptions : throw if no default version is defined
        Caller     : general
        Status     : At Risk
                   : under development
        """
        row = (
            session.query(CoordSystemORM.version.distinct())
                .where(
                    and_(
                        CoordSystemORM.attrib.like(r'%default_version%'),
                        CoordSystemORM.version != None
                    )
                )
                .first()
        )
        if not row:
            raise Exception(f'No default_version coord_system is defined')
        return row[0]
        
    
    @classmethod
    def fetch_top_level(cls, session: Session) -> str:
        raise Exception('Not implemented')
    

    @classmethod
    def get_default_version(cls, session: Session) -> str:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Description: Alias to fetch_default_version. Kept for backward compatibility
        Returntype : str
        Exceptions : throw if no default version is defined
        """
        return cls.fetch_default_version(session)
    

    @classmethod
    def fetch_all_default(cls, session: Session, species_id: int = 1) -> list[CoordSystem]:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Arg [2]    : species_id: int (default 1)
                     The species_id as defined in the core DB
        Example    : cs = CoordSystemAdaptor.fetch_all_default(session);
        Description: Retrieves the default coordinate systems ordered by rank.
        Returntype : list[ensembl.dbsql.CoordSystem]
        Exceptions : if no default coordinate system is found
        Caller     : general
        Status     : At Risk
                   : under development
        """
        rows = (
            session.query(CoordSystemORM)
                .where(CoordSystemORM.attrib.like(r'%default%'))
                .order_by(CoordSystemORM.rank)
                .all()
        )

        if not rows:
            raise Exception(f'No default coord_system is defined')
        
        cs_list = []
        for cs_row in rows:
            toplevel = False
            seqlevel = True if cs_row.attrib and 'sequence_level' in cs_row.attrib else False
            default = True
            cs = CoordSystem(cs_row.name, cs_row.version, cs_row.rank, toplevel, seqlevel, default, cs_row.species_id, cs_row.coord_system_id)
            cs_list.append(cs)

        return cs_list
    

    @classmethod
    def fetch_all_nondefault(cls, session: Session, species_id: int = 1) -> list[CoordSystem]:
        """
        Arg [1]    : session: Session
                     The session object for connecting to the DB
        Arg [2]    : species_id: int (default 1)
                     The species_id as defined in the core DB
        Example    : cs = CoordSystemAdaptor.fetch_all_default(session);
        Description: Retrieves the non-default coordinate systems ordered by rank.
        Returntype : list[ensembl.dbsql.CoordSystem]
        Exceptions : if no coordinate system is found
        Caller     : general
        Status     : At Risk
                   : under development
        """
        rows = (
            session.query(CoordSystemORM)
                .where(CoordSystemORM.attrib.not_like(r'%default%'))
                .order_by(CoordSystemORM.rank)
                .all()
        )

        if not rows:
            raise Exception(f'No non-default coord_system is defined')
        
        cs_list = []
        for cs_row in rows:
            toplevel = False
            seqlevel = True if cs_row.attrib and 'sequence_level' in cs_row.attrib else False
            default = False
            cs = CoordSystem(cs_row.name, cs_row.version, cs_row.rank, toplevel, seqlevel, default, cs_row.species_id, cs_row.coord_system_id)
            cs_list.append(cs)

        return cs_list
       
