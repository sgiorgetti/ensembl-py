from sqlalchemy.orm import Session

from ensembl.core.models import CoordSystem, SeqRegion

def fetch_species_id_by_seq_region_id(session: Session, seq_region_id: int) -> int:
    species_id = (
        session.query(CoordSystem.species_id)
        .join(SeqRegion.coord_system)
        .where(SeqRegion.seq_region_id == seq_region_id)
        .first()
    )
    if species_id:
        return species_id[0]
    return None

def fetch_species_id_by_seq_region_id(session: Session, seq_region_id: int) -> int:
    species_id = (
        session.query(CoordSystem.species_id)
        .join(SeqRegion.coord_system)
        .where(SeqRegion.seq_region_id == seq_region_id)
        .first()
    )
    if species_id:
        return species_id[0]
    return None