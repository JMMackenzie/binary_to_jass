# Binary_to_JASS
This simple tool converts input binary sequences often used in the ds2i/pisa codebase
into an index that JASS can read. This is useful for running experiments where both
indexes are derived from the same base collection, allowing head-to-head comparisons.

This tool is not meant to be well engineered, but is a simple hack that works well
enough. Please adapt to your own use.


Acknowledgements
----------------
This codebase contains functionality found in other codebases. Some of these
codebases, or other related work, is shown below.

- [pisa](https://github.com/pisa-engine/pisa)
- [JASSv2](https://github.com/andrewtrotman/JASSv2)
- [ds2i](https://github.com/ot/ds2i)
- [ATIRE](https://github.com/snapbug/atire)
- [JASS](https://github.com/lintool/JASS)
- [FastPFor](https://github.com/lemire/FastPFor)
- [FastDifferentialCoding](https://github.com/lemire/FastDifferentialCoding)

Note that the real purpose of this library is to take the same input format
used by pisa/ds2i and convert it for usage in JassV2

## Collection input format

A _binary sequence_ is a sequence of integers prefixed by its length, where both
the sequence integers and the length are written as 32-bit little-endian
unsigned integers.

A _collection_ consists of 3 files, `<basename>.docs`, `<basename>.freqs`,
`<basename>.sizes`.

* `<basename>.docs` starts with a singleton binary sequence where its only
  integer is the number of documents in the collection. It is then followed by
  one binary sequence for each posting list, in order of term-ids. Each posting
  list contains the sequence of document-ids containing the term.

* `<basename>.freqs` is composed of a one binary sequence per posting list, where
  each sequence contains the occurrence counts of the postings, aligned with the
  previous file (note however that this file does not have an additional
  singleton list at its beginning).

* `<basename>.sizes` is composed of a single binary sequence whose length is the
  same as the number of documents in the collection, and the i-th element of the
  sequence is the size (number of terms) of the i-th document.

