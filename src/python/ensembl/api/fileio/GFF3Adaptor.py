__all__ = ['GFF', 'GFFAdaptor']

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

"""
TO BE IMPLEMENTED!!!
"""

class GFF3():
    def __init__(self,
                 seq_region,

    ) -> None:
        pass

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(????)'



def main():
    out_file = "your_file.gff"
    seq = Seq("GATCGATCGATCGATCGATC")
    # rec = SeqRecord(seq, "ID1")
    rec = SeqRecord(seq, "ID1", "name", "descr", ["xref1", "xref2"])
    qualifiers = {
        "source": "prediction",
        "score": 10.0,
        "other": ["Some", "annotations"],
        "ID": "gene1",
    }
    sub_qualifiers = {"source": "prediction"}
    top_feature = SeqFeature(
        FeatureLocation(0, 20), type="gene", strand=1, qualifiers=qualifiers
    )
    top_feature.sub_features = [
        SeqFeature(FeatureLocation(0, 5), type="exon", strand=1, qualifiers=sub_qualifiers),
        SeqFeature(
            FeatureLocation(15, 20), type="exon", strand=1, qualifiers=sub_qualifiers
        ),
    ]
    rec.features = [top_feature]

    with open(out_file, "w") as out_handle:
        GFF.write([rec], out_handle)

if __name__ == '__main__':
    main()