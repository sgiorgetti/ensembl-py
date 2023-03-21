from ensembl.api.core.Exon import SplicedExon
from ensembl.api.core.Slice import Slice
from typing import Optional, Union

__all__ = ['Transcript']

class Transcript(object):
    """Representation of a transcript.

    TO DO!!!

    The location of a transcript is always provided 5' -> 3' on the forward strand.
    Thus the start should always be lower than the end.

    Attributes:
      exons: List of exons constituting the transcript.
      start: Start position of the transcript.
      end: End position of the transcript.
      strand: Strand of the transcript, either + or -.
      location_name: Name of the region the transcript is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the transcript.
      external_id: Name of the transcript.
      introns: List of introns constituting the transcript, should be 0 if exons size is one.
      cds_genomic_start: Genomic start position of the CDS, always be lower than cds_genomic_end.
      cds_genomic_end: Genomic end position of the CDS, always be greater than cds_genomic_start.
      cds_sequence: Translateable cDNA sequence.
      translation_sequence: Protein sequence translated from cds_sequence.
    """

    __type = 'transcript'
 
    def __init__(self,
                 stable_id: str,
                 symbol: str,
                 slice: Slice,
                 internal_id: Union[str, int] = None,
                 biotype: str = None,
                 source: str = None
                ) -> None:
        if not stable_id:
            raise ValueError()
        if not symbol:
            raise ValueError()
        if not Slice:
            raise ValueError()
        if stable_id.find('.') < 0:
            raise ValueError()
        (self._unversioned_stable_id, self._version) = stable_id.split('.')
        self._symbol = symbol
        self._slice = slice
        self._metadata = {}
        self._attributes = {}
        self._internal_id = internal_id
        self._biotype = biotype
        self._source = source

    def __repr__(self) -> str:
        if self._unversioned_stable_id:
            if self._attributes.get('is_canonical'):
                return f'{self.__class__.__name__}({self.stable_id} - canonical)'
            return f'{self.__class__.__name__}({self.stable_id})'
        return f'{self.__class__.__name__}(internal_id{self._internal_id})'
    
    @property
    def stable_id(self) -> str:
        return f"{self._unversioned_stable_id}.{self._version}"

    @stable_id.setter
    def stable_id(self, value: str) -> None:
        (self._unversioned_stable_id, self._version) = value.split('.')
    
    @property
    def symbol(self) -> str:
        return self._symbol

    @symbol.setter
    def symbol(self, value: str) -> None:
        self._symbol = value
    
    @property
    def type(self) -> str:
        return self.__type

    @property
    def unversioned_stable_id(self) -> str:
        return self._unversioned_stable_id

    @unversioned_stable_id.setter
    def unversioned_stable_id(self, value) -> None:
        self._unversioned_stable_id = value

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, value) -> None:
        self._version = value

    @property
    def slice(self) -> Slice:
        return self._slice
    
    @slice.setter
    def slice(self, value: Slice) -> None:
        self._slice = value

    @property
    def internal_id(self) -> int:
        return self._internal_id
    
    @property
    def dbID(self) -> int:
        return self._internal_id
    
    @property
    def biotype(self) -> str:
        return self._biotype

    @biotype.setter
    def biotype(self, value: str) -> None:
        self._biotype = value

    @property
    def source(self) -> str:
        return self._source

    @source.setter
    def source(self, value: str) -> None:
        self._source = value

    def get_exons(self) -> list[SplicedExon]:
        return self._exons
    
    def set_exons(self, exons: list[SplicedExon]) -> None:
        self._exons = exons

    def get_introns(self) -> tuple:
        return self._introns
    
    def set_introns(self, introns: list) -> None:
        self._introns = introns
    
    def get_metadata(self, md: Optional[str] = None) -> Union[dict, tuple]:
        if md is None:
            return self._metadata
        return tuple(md, self._metadata.get(md))
    
    def add_metadata(self, code: str, value: str):
        self._attributes[code] = value
    
    def set_attribs(self, attribs: dict[str, str]) -> None:
        self._attributes = attribs

    def get_attribs(self, att_code: str = None) -> Union[dict, tuple]:
        if att_code is None:
            return self._attributes
        return tuple(att_code, self._attributes.get(att_code))
    
    def add_attrib(self, code: str, value: str):
        self._attributes[code] = value

    def is_canonical(self) -> bool:
        if self._attributes.get('is_canonical'):
            return True
        return False
    
    def is_mane_select(self) -> bool:
        if self._attributes.get('MANE_Select'):
            return True
        return False
    
    def is_mane(self) -> bool:
        for att_code in self._attributes.keys().lower():
            if 'mane' in att_code:
                return True
        return False


    def get_summary(self) -> dict[str, str]:
        """
        Example       : transcript_summary = transcript.get_summary()
        Description   : Retrieves a textual summary of this Feature.
        Returns       : Dict[str, str]
        Status        : Alpha - Intended for internal use
        """
        summary = {}
        summary['id'] = self.unversioned_stable_id
        summary['transcript_id'] = self.unversioned_stable_id
        summary['version'] = self.version if self.version else ''
        summary['start'] = self._slice.location.start
        summary['end'] = self._slice.location.end
        summary['strand'] = self._slice.strand.value
        summary['seq_region_name'] = self._slice.region.name

        summary['description'] = self._attributes.get('description').value
        summary['biotype'] = self._biotype
        summary['source'] = self._source
        summary['symbol'] = self._symbol
        summary['logic_name'] = self._attributes.get('analysis').value
        if self._attributes.get('ccds_transcript'):
            summary['ccdsid'] = self._attributes.get('ccds_transcript').value
        if self._attributes.get('tsl'):
            summary['transcript_support_level'] = self._attributes.get('tsl').value
        if self._attributes.get('gencode_basic'):
            summary['tag'] = 'basic'

        return summary
    
# 1       havana  mRNA    65419   71585   .       +       .       ID=transcript:ENST00000641515;Parent=gene:ENSG00000186092;Name=OR4F5-201;biotype=protein_coding;tag=basic,Ensembl_canonical,MANE_Select;transcript_id=ENST00000641515;version=2
    def gff3_qualifiers(self) -> dict[str, Union[str, tuple[str]]]:

        tags = ['basic']
        if self.is_canonical():
            tags.append('Ensembl_canonical')
        if self.is_mane_select():
            tags.append('MANE_Select')
        
        qualifiers = {
            'source': self._source,
            'score': ".",
            'ID': f"{self.__type}:{self._unversioned_stable_id}",
            'Parent': "gene:ENSG00000186092",
            'Name': self._symbol,
            'biotype': self._biotype,
            'transcript_id': self._unversioned_stable_id,
            'version': self._version
        }
        qualifiers['tag'] = tuple(tags)

        return qualifiers
