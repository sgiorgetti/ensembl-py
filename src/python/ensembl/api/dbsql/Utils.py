from time import time
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

def timeme(func):
    # This function shows the execution time of 
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s')
        return result
    return wrap_func