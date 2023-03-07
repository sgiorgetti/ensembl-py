__all__ = ['GFF', 'GFFAdaptor']

from BCBio import GFF
from BCBio.GFF.GFFOutput import GFF3Writer as BCBioGFF3, write, _IdHandler
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from ensembl.api.core.Strand import Strand
from ensembl.api.core.Slice import Slice

from typing import Optional

"""
TO BE IMPLEMENTED!!!
"""

class GFF3Writer(BCBioGFF3):
    def __init__(self,
                 custom_dirs: dict[str, str] = None
    ) -> None:
        BCBioGFF3.__init__(self)
        self._custom_dirs = None
        if custom_dirs:
            self._custom_dirs = custom_dirs

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(????)'
    
    def _write_header(self, out_handle):
        """Write out standard header directives.
           This overrides the Biopython one by adding
           Ensembl standard directives. 
        """
        # out_handle.write("##gff-version 3\n")
        super()._write_header(out_handle)
        if self._custom_dirs:
            # WATCH THE ORDER!!!
            for k, v in self._custom_dirs.items():
                out_handle.write(f"{k} {v}\n")
        
    def _write_rec(self, rec, out_handle):
        pass


def write(recs, out_handle, include_fasta=False, custom_dirs: dict[str, str] = None):
    """High level interface to write GFF3 files from SeqRecords and SeqFeatures.
    If include_fasta is True, the GFF3 file will include sequence information
    using the ##FASTA directive.
    """
    if custom_dirs:
        writer = GFF3Writer(custom_dirs)
    else:
        writer = GFF3Writer()
    return writer.write(recs, out_handle, include_fasta)


def main():
    out_file = "your_file.gff"
    ensembl_dirs = {
        "##sequence-region  ": "ID1 1 999999",
        "#!genome-build": "accessioning_body assembly.id",
        "#!genome-version": "assembly.name",
        "#!genome-date": "assembly_date",
        "#!genome-build-accession": "accession_id",
        "#!genebuild-last-updated": "gb_last_geneset_update",
    }
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
        write([rec], out_handle, custom_dirs=ensembl_dirs)

if __name__ == '__main__':
    main()