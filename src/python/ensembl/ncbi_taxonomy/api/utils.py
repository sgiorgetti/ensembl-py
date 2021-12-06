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

from sqlalchemy import and_
from sqlalchemy.orm import as_declarative, Session, aliased
from sqlalchemy.orm.exc import NoResultFound

from ensembl.ncbi_taxonomy.models import NCBITaxaNode, NCBITaxonomy


@as_declarative()
class Taxonomy():
    """Contains all the taxonomy related functions over NCBITaxonomy ORM.

    Attributes:
        session: db Session()
    """


    @classmethod
    def fetch_node_by_id(cls, session: Session, taxon_id: int) -> NCBITaxonomy:
        """Returns taxonomy node object by taxon_id

        Args:
            taxon_id: Unique taxonomy identifier in database

        Raises:
            NoResultFound() if taxon_id does not exist
        """
        q = (
            session.query(NCBITaxonomy)
            .filter(NCBITaxonomy.taxon_id == taxon_id)
            .first()
        )
        if not q:
            raise NoResultFound()
        return q


    @classmethod
    def fetch_taxon_by_species_name(cls, session: Session, name: str) -> NCBITaxonomy:
        """Returns taxon_id integer matching name

        Args:
            name: Scientific ncbi_taxa_name.name in database

        Raises:
            NoResultFound() if taxon_id does not exist or returns multiple results
        """
        name.replace("_", " ")
        return (
            session.query(NCBITaxonomy)
            .filter(NCBITaxonomy.name.like(name))
            .filter(NCBITaxonomy.name_class == "scientific name")
            .one()
        )


    @classmethod
    def parent(cls, session: Session, taxon_id: int) -> NCBITaxonomy:
        """Returns taxonomy node object for parent node

        Args:
            taxon_id: Unique taxonomy identifier in database

        Raises:
            NoResultFound() if taxon_id does not exist
        """
        ParentTaxonomy = aliased(NCBITaxonomy, name="parent_ncbi_taxonomy")
        q = (
            session.query(NCBITaxonomy, ParentTaxonomy)
            .outerjoin(
                (ParentTaxonomy, NCBITaxonomy.parent_id == ParentTaxonomy.taxon_id)
            )
            .filter(NCBITaxonomy.taxon_id == taxon_id)
            .filter(ParentTaxonomy.name_class == "scientific name")
            .first()
        )
        try:
            return q[1]
        except TypeError:
            raise NoResultFound()


    @classmethod
    def children(cls, session: Session, taxon_id: int) -> tuple:
        """Returns taxonomy node object for children nodes

        Args:
            taxon_id: Unique taxonomy identifier in database

        Raises:
            NoResultFound() if taxon_id does not exist or has no children
        """
        q = (
            session.query(NCBITaxonomy)
            .filter(NCBITaxonomy.parent_id == taxon_id)
            .filter(NCBITaxonomy.name_class == "scientific name")
            .all()
        )
        results = list(q)
        rows = [x.__dict__ for x in results]
        q = tuple(rows)
        if not q:
            raise NoResultFound()
        return q


    @classmethod
    def is_root(cls, session: Session, taxon_id: int) -> bool:
        """Returns true if taxon_id is a root and false if not

        Args:
            taxon_id: Unique taxonomy identifier in database
        """
        try:
            if (
                session.query(NCBITaxaNode)
                .filter(NCBITaxaNode.root_id == taxon_id, NCBITaxaNode.taxon_id == taxon_id)
                .one()
            ):
                return True
        except NoResultFound:
            return False


    @classmethod
    def num_descendants(cls, session: Session, taxon_id: int) -> int:
        """Returns number of descendants from taxon_id

        Args:
            taxon_id: Unique taxonomy identifier in database

        Raises:
            NoResultFound() if taxon_id does not exist
        """
        session.query(NCBITaxaNode).filter(NCBITaxaNode.taxon_id == taxon_id).one()
        right_index = (
            session.query(NCBITaxaNode.right_index)
            .filter(NCBITaxaNode.taxon_id == taxon_id)
            .scalar()
        )
        left_index = (
            session.query(NCBITaxaNode.left_index)
            .filter(NCBITaxaNode.taxon_id == taxon_id)
            .scalar()
        )
        return (right_index - left_index - 1) / 2


    @classmethod
    def is_leaf(cls, session: Session, taxon_id: int) -> bool:
        """Returns true if taxon_id is a leaf and false if not

        Args:
            taxon_id: Unique taxonomy identifier in database

        Raises:
            NoResultFound() if taxon_id does not exist
        """
        if cls.num_descendants(session, taxon_id) == 0:
            return True
        return False

    @classmethod
    def fetch_ancestors(cls, session: Session, taxon_id: int) -> tuple:
        """Returns a tuple of ancestor node objects from taxon_id

        Args:
            taxon_id: Unique taxonomy identifier in database

        Raises:
            NoResultFound() if taxon_id does not exist or has no ancestors
        """
        ParentTaxaNode = aliased(NCBITaxaNode)
        q = (
            session.query(ParentTaxaNode, NCBITaxaNode)
            .outerjoin(
                NCBITaxaNode,
                and_(
                    NCBITaxaNode.left_index.between(
                        ParentTaxaNode.left_index, ParentTaxaNode.right_index
                    ),
                    ParentTaxaNode.taxon_id != NCBITaxaNode.taxon_id,
                ),
            )
            .filter(NCBITaxaNode.taxon_id == taxon_id)
            .all()
        )
        if not q:
            raise NoResultFound()
        results = []
        for row in q:
            taxon = row[0].__dict__
            results.append(taxon)
        ordered_results = sorted(results, key=lambda x: x["taxon_id"])
        q = tuple(ordered_results)
        return q


    @classmethod
    def all_common_ancestors(
        cls, session: Session, taxon_id_1: int, taxon_id_2: int
    ) -> tuple:
        """Returns a tuple of common ancestor node objects shared between taxa

        Args:
            taxon_id_1: Unique taxonomy identifier in database
            taxon_id_2: Unique taxonomy identifier in database

        Raises:
            NoResultFound() if taxon_ids do not exist or have no common ancestors
        """
        ancestors_1 = cls.fetch_ancestors(session, taxon_id_1)
        ancestors_2 = cls.fetch_ancestors(session, taxon_id_2)
        if ancestors_1 is None or ancestors_2 is None:
            raise NoResultFound()
        ancestors_1 = list(ancestors_1)
        ancestors_2 = list(ancestors_2)
        ancestors_ids_1 = [taxon["taxon_id"] for taxon in ancestors_1]
        ancestors_ids_2 = [taxon["taxon_id"] for taxon in ancestors_2]
        common_ancestors = list(set(ancestors_ids_1).intersection(ancestors_ids_2))
        common_ancestors.sort(
            key=lambda taxon_id: (-cls.num_descendants(session, taxon_id), taxon_id)
        )
        all_common_ancs = [
            cls.fetch_node_by_id(session, taxon_id) for taxon_id in common_ancestors
        ]
        return tuple(all_common_ancs)


    @classmethod
    def last_common_ancestor(
        cls, session: Session, taxon_id_1: int, taxon_id_2: int
    ) -> NCBITaxonomy:
        """Returns most recent common ancestor node object shared between taxa

        Args:
            taxon_id_1: Unique taxonomy identifier in database
            taxon_id_2: Unique taxonomy identifier in database

        Raises:
            NoResultFound() if taxon_ids do not exist or have no common ancestors
        """
        common_ancestors = cls.all_common_ancestors(session, taxon_id_1, taxon_id_2)
        return common_ancestors[0]
