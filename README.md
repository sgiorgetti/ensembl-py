# ensembl-py

Python Ensembl Main CORE code source repository


# Aim

Centralise generic code use across all other project within Ensembl, more particularly DBS access layers, Reusable Hive 
Component etc.. 
 

# Ensembl Python API (experimental/basic-genes branch)?

This is a super-early, exploratory work for assessing the viablity of "porting" the Perl API into Python.
A few words of warning are due, at this point. First and foremost, naive "porting" would make no sense:
- some code features would not be needed anymore
- we'd like to leave some historical hackings and peculiarities behind
- the target (microservices) architecture differs significantly from the original one
- the target system/process architecture for the Ensembl release lifecycle (handover) is still to be fully defined

## Scope

Scope is limited to the few core features, and related ancillary information, we have to deal with. Namely,
- Gene
- Transcript
- Exon
- Translation/Protein
- Limited subset of metadata/Xrefs (TBD)
- Link to sequence
- Adaptors to MySQL core schema
- Adaptors to std file formats (GTF/GFF3/FASTA ... and - if any?)

## Goals

Achieve an understanding of
- Comprehensiveness of the new Core Data Model [CDM] for internal purposes (not web/externally-facing GraphQL)
- I/O needed complexity, especially when dealing with files as storage system
- Possible scoping and phasing "porting" features from Perl API to Python
- Try to mitigate the proliferation of "Python APIs" doing pretty much the same thing over and over again

## Approach (first steps so far)

### Unfruitful starting points so far

Started initially from `ensembl-io` with the idea to port/re-implement simple I/O features.
This turned out to be tricky, because `ensembl-io` is comfortably sitting on and embedding the core API, including the part of the glorious "projection" complexity.
Anyway, the idea is still that one, but some "foundation" is required to start building and thinking on I/O.

Thus, I came to `ensembl-py` [ensembl-py-url]. This repo contains - among a few other bits - an ORM implementation [SQLAlchemy] for the core schema.
I thought I could start there, but it turned out to be too "poor" in features: the ORM maps SQL tables into Python classes and provides convenient framework for I/O with RDBMS (MySQL for us).
Extending the Gene class, which mirrors the Gene table, is feasible, but may result in tying up a "Gene" abstract concept to a very concrete and storage-aware implementation of it. Not very nice.

On the other hand, the `ensembl-genes-api` [ensembl-genes-api-url] repo contains better structured code for defining an abstract Gene - unsurprisingly it's GBuilders' repo - but it is very focussed and related to file-based storage only.

To make things more complicated, neither of these approaches complies with the new [CDM] specs.

### Current/latest try

Decided to go for a familiar approach Object -- ObjectAdaptor(s), instead of something new, which would require a learning/acceptance curve. Let's start from - say - a gene in a GTF/GFF3/FASTA(/BED?) sort of context.
I think I need:
- Gene class - modelled w.r.t. CDM, but must be storage-agnostic, similarly to the `Gene.pm`
- Min set of ancillary classes to Gene - e.g. Metadata - required to sensibly define a Gene.
- Transcript class: a gene seems to be a set of transcripts, after all
- Sequence (class?)
- Adaptor classes for RBDMS and files
- Standard libraries (Biopython, Pysam, ...) for basic features
- ???

> Watch the open PRs on the CDM repository
> Gene and Metadata concepts may be (are) affected
> There are still inconsistencies - possibly minor ones - around RE types and constraints


### Current "sprint" goal
- Gene class compliant with `CDM`, `ensembl-genes-api`, and `ensembl` Gene concepts
- Minimal Gene adaptor for MySQL, using SQL Alchemy's one from `ensembl-py`
- Minimal Gene adaptor for GTF file, using something as close as possible to `ensembl-genes-api`'s
- What's missing? Does this approach make sense?

[//]: # (List of references used before)

  [CDM]: <https://github.com/Ensembl/ensembl-cdm-docs>
  [ensembl-genes-api-url]: <https://github.com/Ensembl/ensembl-genes-api>
  [ensembl-py-url]: <https://github.com/Ensembl/ensembl-py>
  [SQLAlchemy]: <https://www.sqlalchemy.org>