# ensembl-py

Python Ensembl Main CORE code source repository


# Aim

Centralise generic code use across all other project within Ensembl, more particularly DBS access layers, Reusable Hive 
Component etc.. 
 

# Ensembl Python API?

This is a super-early, exploratory work for assessing the viablity of "porting" the Perl API into Python.
A few words of warning are due, at this point. First and foremost, naive "porting" would make no sense:
- some code features would not be needed anymore
- we'd like to leave some historical hackings and peculiarities behind
- the target (microservices) architecture differs significantly from the original one
- the target system/process architecture for the Ensembl release lifecycle (handover) is still to be fully defined

## Scope

Scope is limited to a few core features, and related ancillary information, we have to deal with. Namely,
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

### experimental/basic-genes branch

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

### Goals and results
- Gene class compliant with `CDM`, `ensembl-genes-api`, and `ensembl` Gene concepts
- Minimal Gene adaptor for MySQL, using SQL Alchemy's one from `ensembl-py`
- Minimal Gene adaptor for GTF file, using something as close as possible to `ensembl-genes-api`'s
- What's missing? Does this approach make sense?

#### On the CDM
Because of its design, the new `CDM` is very functional, end-user oriented: it tries to capture and present information ... provided the information are there.

This is absolutely fine for the final data consumer (the web site), but it adds some un-needed complexity in the process of preparing data, during the Release lifecycle.

#### On the flat files IO
With reference to IO libraries available in Python, it looks like there are some quite adequate to use a (R)DBMS, but things are different when going to flat files (GFF, FASTA, ...).

In truth there are a few approaches, but could not find many or well-supported.
`gffutils` for instance looks nice, but it has some disadvantages like being based on `SQLite`, which needs to be created, maintained ..., and provides little or no features for dumping data out of the `SQLite` DB back into flat file format.

#### What's missing? Does this approach make sense?
Well, a lot is missing ... notably the sequences and related (basic) functionalities, as well as many slice/projection-related tools.

The approach does not make a lot of sense:
- use of `CDM` makes the data building phase cumbersome
- use of a single `data model` for building/presenting data seems not working very well
- some information are intrinsically related ... is the flat-file-only the only and right way to go? Isn't it a one-size-fits-all approach?

### Latest try on `experimental/genes-api`
In this exercise, the target scenario is related to the `IDMapping` re-implementation, which comprises
- Re-write the Perl scripts/pipeline with Python
- Avoid RDBMS core DBs, as much as possible
- Replace `exonerate` with `minimap2` directly or via `Liftoff`
- Basic feature classes compatible with `ensembl`  concepts
- Minimal feature adaptors for MySQL, using SQL Alchemy's one from `ensembl-py`
- Building the ID Mapping Cache file ==> THIS MAKES NO SENSE, as it would require the full mapping and CS features of the API. Better keep the cache building as it is, and re-implement the ID Mapping logic with new tech! Also, `Liftoff` works with standard GFF files, and there is no need to use a custom dump from core DBs.
- Basic sequence support


### Repurposing/refocussing
- Minimal feature adaptors for GFF/GTF file, using something as close as possible to `ensembl-genes-api`'s
- This should be preliminary to create GFF/GTF serializers and the mapping between storage systems (files and DB) and the data model


[//]: # (List of references used before)

  [CDM]: <https://github.com/Ensembl/ensembl-cdm-docs>
  [ensembl-genes-api-url]: <https://github.com/Ensembl/ensembl-genes-api>
  [ensembl-py-url]: <https://github.com/Ensembl/ensembl-py>
  [SQLAlchemy]: <https://www.sqlalchemy.org>