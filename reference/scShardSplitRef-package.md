# scShardSplitRef: Reference Genome Splitting for Single-Cell Multiome Data

Splits chromosomes of large-genome species into smaller shards that
respect the maximum contig size supported by 10x Genomics' 'Cell Ranger
ARC' reference builder, avoiding breaks inside gene bodies, and rewrites
the corresponding 'GTF' annotation coordinates to match. This enables
researchers to build custom 'Cell Ranger ARC' references for Multiome
'ATAC' and Gene Expression data from genomes that would otherwise exceed
the tool's chromosome size limits.

## See also

Useful links:

- <https://irinavkuznetsova.github.io/scShardSplitRef/>

- <https://github.com/IrinaVKuznetsova/scShardSplitRef>

- Report bugs at
  <https://github.com/IrinaVKuznetsova/scShardSplitRef/issues>

## Author

**Maintainer**: Irina Kuznetsova <irina.kuznetsova@curtin.edu.au>
([ORCID](https://orcid.org/0009-0001-5411-6063))

Authors:

- Irina Kuznetsova <irina.kuznetsova@curtin.edu.au>
  ([ORCID](https://orcid.org/0009-0001-5411-6063))

- Luke Pembleton <lpembleton@barenbrug.com.au>
  ([ORCID](https://orcid.org/0000-0001-8131-8210))

- Luca Curci <pasqualeluca.curci@cnr.it>
  ([ORCID](https://orcid.org/0000-0002-5780-3029))

- Adam H. Sparks <adamhsparks@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-0061-8359))

Other contributors:

- Curtin University ([ROR](https://ror.org/02n415q13)) \[funder,
  copyright holder\]

- Grains Research and Development Corporation
  ([ROR](https://ror.org/02xwr1996)) (GRDC Project CUR2210-005OPX
  (AAGI-CU)) \[funder, copyright holder\]
