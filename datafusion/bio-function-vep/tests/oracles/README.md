# REF mismatch source oracle

`ref_mismatch.pl` invokes the unmodified Ensembl VEP `codon`, `display_codon`,
`peptide` and `_get_alternate_cds` methods on minimal coordinate/cache objects.
It does not replace the sequence methods under test. This isolates the CDS and
RNA-edit rules independently of vepyr, for porting-tests issue 221.

The checked-in TSV was produced by VEP 116.0 (`57ea5c52340acc1f156267f810ad162e26597082`)
and ensembl-variation `2fb834b987ede3824e200197a838ce11e91aeb4b` in this image:

```sh
docker run --rm --platform linux/amd64 \
  -v "$PWD/datafusion/bio-function-vep/tests/oracles:/oracle:ro" \
  --entrypoint perl \
  ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11 \
  /oracle/ref_mismatch.pl > /tmp/ref_mismatch.tsv
diff -u datafusion/bio-function-vep/tests/oracles/ref_mismatch.tsv /tmp/ref_mismatch.tsv
```

The first 16 rows are strand, CDS start, REF, ALT, display codons, amino acids.
The last four rows are attribute case, reference/alternate codon, and their
display forms. Explicit RNA-edit attributes include poly-A; BAM state alone
does not activate the reference replacement branch. Rust regression tests
assert these exact results. End-to-end VCF, HGVSp and DOMAINS coverage uses the
full Docker oracle and a trimmed cache in the companion vepyr PR.
