"""Regenerate the golden files for tests/golden.rs.

Needs seqkit and Biopython on PATH; the checked-in files were made with
seqkit 2.14.0 and Biopython 1.88:

    conda create -n seqgolden -c conda-forge -c bioconda seqkit=2.14.0 biopython=1.88
    conda run -n seqgolden python make_golden.py

Outputs (next to this script):
  corpus.fastq     random reads: IUPAC codes in both cases plus `-` and `.` gaps
  gc.tsv           seqkit fx2tab -g (acgt) and -B GC (length), in percent
  qual.tsv         seqkit fx2tab -q (error-probability mean)
  revcomp.tsv      seqkit seq -r -p and Biopython reverse_complement
  coding.fasta     random sequences over ACGT, some IUPAC codes and lowercase
  translate.tsv    seqkit translate and Biopython Seq.translate per frame/table

Known oracle differences, and how they are resolved here:
  - Biopython returns B, Z or J for an ambiguous codon whose expansions are
    {D,N}, {E,Q} or {I,L}; seqkit returns X. The kernel follows seqkit, so
    Biopython's output has B, Z and J mapped to X.
  - Biopython rejects codons with a gap or an unknown byte; the coding corpus
    has neither (unit tests cover them).
  - Biopython refuses to_stop for tables 27, 28 and 31, and seqkit has no
    to_stop; those combinations are skipped.
  - seqkit 2.14.0 has no tables 15, 32 and 33; its column is empty for them.
  - seqkit 2.14.0 translate does not complement IUPAC ambiguity codes in
    negative frames (`RTT` at frame -1 gives K, from AAR, instead of N, from
    AAY), although `seqkit seq -r -p` complements them correctly. The kernel
    follows Biopython there; golden.rs compares negative frames to Biopython
    only.
"""

import random
import subprocess
import sys
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
rng = random.Random(20261115)

IUPAC = "ACGTRYSWKMBDHVN"
SEQ_ALPHABET = "ACGT" * 6 + IUPAC + IUPAC.lower() + "acgt" * 3 + "-."
CODING_ALPHABET = "ACGT" * 10 + "acgt" * 2 + "NRYKMSW"
SEQKIT_MISSING_TABLES = {15, 32, 33}


def seqkit(args, stdin):
    return subprocess.run(
        ["seqkit", *args], input=stdin, capture_output=True, text=True, check=True
    ).stdout


def write_corpus():
    records = []
    for i in range(300):
        n = rng.choice([1, 2, 3, 5, 10, 50, 151, 300])
        seq = "".join(rng.choice(SEQ_ALPHABET) for _ in range(n))
        qual = "".join(chr(33 + rng.randint(0, 41)) for _ in range(n))
        records.append((f"r{i}", seq, qual))
    fastq = "".join(f"@{r}\n{s}\n+\n{q}\n" for r, s, q in records)
    (HERE / "corpus.fastq").write_text(fastq)
    return records, fastq


def table_of(out):
    lines = [l.split("\t") for l in out.splitlines() if l]
    header = lines[0]
    return [dict(zip(header, l)) for l in lines[1:]]


def write_gc_qual(fastq):
    acgt = table_of(seqkit(["fx2tab", "-n", "-H", "-g"], fastq))
    length = table_of(seqkit(["fx2tab", "-n", "-H", "-B", "GC"], fastq))
    qual = table_of(seqkit(["fx2tab", "-n", "-H", "-q"], fastq))
    with open(HERE / "gc.tsv", "w") as fh:
        fh.write("id\tacgt_pct\tlength_pct\n")
        for a, b in zip(acgt, length):
            fh.write(f"{a['#name']}\t{a['GC']}\t{b['GC']}\n")
    with open(HERE / "qual.tsv", "w") as fh:
        fh.write("id\tavg_qual\n")
        for q in qual:
            fh.write(f"{q['#name']}\t{q['avg.qual']}\n")


def write_revcomp(records, fastq):
    rc = {}
    out = seqkit(["seq", "-r", "-p", "-t", "dna"], fastq)
    lines = out.splitlines()
    for i in range(0, len(lines), 4):
        rc[lines[i][1:]] = lines[i + 1]
    with open(HERE / "revcomp.tsv", "w") as fh:
        fh.write("id\tseqkit\tbiopython\n")
        for r, s, _ in records:
            fh.write(f"{r}\t{rc[r]}\t{Seq(s).reverse_complement()}\n")


def write_translate():
    seqs = []
    for i in range(40):
        n = rng.randint(0, 90) + rng.choice([0, 0, 120])
        seqs.append((f"c{i}", "".join(rng.choice(CODING_ALPHABET) for _ in range(n))))
    # Keep every sequence non-empty so seqkit emits a record for each.
    seqs = [(r, s or "ATG") for r, s in seqs]
    fasta = "".join(f">{r}\n{s}\n" for r, s in seqs)
    (HERE / "coding.fasta").write_text(fasta)

    tables = sorted(CodonTable.unambiguous_dna_by_id)
    combos = [(1, t, False) for t in tables]
    combos += [(f, t, s) for t in (1, 2, 11) for f in (1, 2, 3, -1, -2, -3) for s in (False, True)]
    combos = sorted(set(combos), key=lambda c: (c[1], c[0], c[2]))

    def bio(seq, frame, table, to_stop):
        s = Seq(seq)
        if frame < 0:
            s = s.reverse_complement()
        s = s[abs(frame) - 1 :]
        s = s[: len(s) // 3 * 3]
        aa = str(s.translate(table=table, to_stop=to_stop))
        return aa.translate(str.maketrans("BZJ", "XXX"))

    def seqkit_tr(frame, table):
        out = seqkit(["translate", "-w", "0", "-f", str(frame), "-T", str(table)], fasta)
        lines = out.splitlines()
        res = {}
        for i in range(0, len(lines), 2):
            # seqkit names a negative-frame record "id_frame=-1"; keep the id.
            res[lines[i][1:].split()[0].split("_frame=")[0]] = lines[i + 1]
        return res

    with open(HERE / "translate.tsv", "w") as fh:
        fh.write("id\tframe\ttable\tto_stop\tseqkit\tbiopython\n")
        for frame, table, to_stop in combos:
            sk = {} if to_stop or table in SEQKIT_MISSING_TABLES else seqkit_tr(frame, table)
            for r, s in seqs:
                b = bio(s, frame, table, to_stop)
                fh.write(f"{r}\t{frame}\t{table}\t{str(to_stop).lower()}\t{sk.get(r, '')}\t{b}\n")


def main():
    records, fastq = write_corpus()
    write_gc_qual(fastq)
    write_revcomp(records, fastq)
    write_translate()
    print("ok", file=sys.stderr)


if __name__ == "__main__":
    main()
