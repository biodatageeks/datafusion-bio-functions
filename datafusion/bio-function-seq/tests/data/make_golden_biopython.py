"""Regenerate the Biopython and Trimmomatic golden files for tests/golden_biopython.rs.

Needs Biopython; the checked-in files were made with Biopython 1.88:

    uv run --no-project --with biopython==1.88 python make_golden_biopython.py

The Trimmomatic file needs Java and the Trimmomatic jar (made with 0.39):

    TRIMMOMATIC_JAR=/path/to/trimmomatic-0.39.jar \
        uv run --no-project --with biopython==1.88 python make_golden_biopython.py

Without TRIMMOMATIC_JAR the existing trimmomatic.tsv is left untouched.

Inputs (next to this script, written by make_golden.py):
  corpus.fastq     random reads: IUPAC codes in both cases plus `-` and `.` gaps
  coding.fasta     random sequences over ACGT, some IUPAC codes and lowercase

Outputs:
  bio_cds.tsv            Seq.translate(cds=True): the protein, or `!` and the
                         failed check, for every NCBI table
  bio_translate_raw.tsv  Seq.translate per frame and table, with Biopython's
                         B, Z and J letters left in place
  bio_transcribe.tsv     Seq.transcribe and Seq.back_transcribe
  bio_motif.tsv          Bio.SeqUtils.nt_search positions (0-based)
  bio_weight.tsv         Bio.SeqUtils.molecular_weight, empty where it raises
  bio_gc_skew.tsv        Bio.SeqUtils.GC_skew
  trimmomatic.tsv        reads kept by Trimmomatic SE for several step lists

Notes:
  - Biopython reads `X` in a nucleotide sequence as "any base"; the kernels
    treat it as an unknown byte. The corpora do not use it.
  - Biopython rejects a codon with a gap, a dot or an unknown byte, and the
    kernels translate it to X or -. The CDS corpus has none.
"""

import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqUtils import GC_skew, molecular_weight, nt_search

HERE = Path(__file__).resolve().parent
rng = random.Random(20261124)

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
IUPAC = "ACGTRYSWKMBDHVN"
CDS_FAILURES = {
    "First codon": "no_start_codon",
    "Sequence length": "length_not_multiple_of_3",
    "Final codon": "no_stop_codon",
    "Extra in frame stop codon": "internal_stop_codon",
}


def records(name):
    lines = (HERE / name).read_text().splitlines()
    step = 4 if name.endswith(".fastq") else 2
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), step)]


def write(name, header, rows):
    with open(HERE / name, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")
    print(f"{name}: {len(rows)} rows", file=sys.stderr)


def cds_result(seq, table):
    try:
        return str(Seq(seq).translate(table=table, cds=True))
    except TranslationError as e:
        for prefix, name in CDS_FAILURES.items():
            if str(e).startswith(prefix):
                return "!" + name
        return None  # an invalid codon: outside what the kernels promise


def mutate(seq, pos, letters):
    return seq[:pos] + rng.choice(letters) + seq[pos + 1 :]


def write_cds():
    rows = []
    for table_id in sorted(CodonTable.unambiguous_dna_by_id):
        table = CodonTable.unambiguous_dna_by_id[table_id]
        sense = [c for c in CODONS if c not in table.stop_codons]
        candidates = []
        for _ in range(8):
            body = [rng.choice(sense) for _ in range(rng.randint(0, 20))]
            cds = rng.choice(table.start_codons) + "".join(body) + rng.choice(table.stop_codons)
            candidates.append(cds)
            n = len(cds)
            # One failure of each kind, then ambiguity and case variants.
            candidates.append(rng.choice([c for c in sense if c not in table.start_codons]) + cds[3:])
            candidates.append(cds[: n - rng.choice([1, 2])])
            candidates.append(cds[: n - 3] + rng.choice(sense))
            if body:
                k = 3 * rng.randint(1, len(body))
                candidates.append(cds[:k] + rng.choice(table.stop_codons) + cds[k:])
            candidates.append(mutate(cds, rng.randrange(n), IUPAC))
            candidates.append(mutate(cds, rng.randrange(3), "RYSWKMN"))
            candidates.append(mutate(cds, n - 1 - rng.randrange(3), "RYSWKMN"))
            candidates.append(cds.lower())
            candidates.append(cds.replace("T", "U"))
        candidates += ["", "A", "ATG", "ATGTAA", "TAA"]
        for seq in candidates:
            result = cds_result(seq, table_id)
            if result is not None:
                rows.append((table_id, seq, result))
    write("bio_cds.tsv", ["table", "seq", "biopython"], rows)


def write_translate_raw():
    rows = []
    for table in sorted(CodonTable.unambiguous_dna_by_id):
        for frame in (1, 2, 3, -1, -2, -3):
            for name, seq in records("coding.fasta"):
                s = Seq(seq)
                if frame < 0:
                    s = s.reverse_complement()
                s = s[abs(frame) - 1 :]
                s = s[: len(s) // 3 * 3]
                rows.append((name, frame, table, str(s.translate(table=table))))
    write("bio_translate_raw.tsv", ["id", "frame", "table", "biopython"], rows)


def write_transcribe():
    rows = [
        (name, Seq(seq).transcribe(), Seq(seq).back_transcribe())
        for name, seq in records("corpus.fastq")
    ]
    write("bio_transcribe.tsv", ["id", "transcribe", "back_transcribe"], rows)


def write_motif():
    rows = []
    for _ in range(300):
        alphabet = "ACGT" * 8 + rng.choice(["", "", "N", "RYN"])
        seq = "".join(rng.choice(alphabet) for _ in range(rng.choice([0, 1, 5, 30, 200])))
        pattern = "".join(
            rng.choice("ACGT" * 4 + IUPAC) for _ in range(rng.randint(1, 6))
        )
        # A short pattern taken from the sequence, so that hits are common.
        if len(seq) > 6 and rng.random() < 0.5:
            k = rng.randrange(len(seq) - 4)
            pattern = "".join(c if c in "ACGT" else "N" for c in seq[k : k + rng.randint(1, 4)])
        positions = nt_search(seq, pattern)[1:]
        rows.append((seq, pattern, ",".join(map(str, positions))))
    write("bio_motif.tsv", ["seq", "pattern", "positions"], rows)


def write_weight():
    amino = "ACDEFGHIKLMNPQRSTVWY"
    cases = []
    for _ in range(40):
        n = rng.choice([0, 1, 2, 10, 60, 300])
        cases.append(("DNA", "".join(rng.choice("ACGT") for _ in range(n))))
        cases.append(("RNA", "".join(rng.choice("ACGU") for _ in range(n))))
        cases.append(("protein", "".join(rng.choice(amino) for _ in range(n))))
    cases += [
        ("DNA", "acgtACGT"),
        ("DNA", "ACGN"),
        ("DNA", "ACGU"),
        ("DNA", "AC-GT"),
        ("RNA", "ACGT"),
        ("RNA", "acgu"),
        ("protein", "MKVX"),
        ("protein", "MKV*"),
        ("protein", "MKVB"),
        ("protein", "MOUK"),
        ("protein", "mkv"),
    ]
    rows = []
    for seq_type, seq in cases:
        for double_stranded in (False, True):
            if double_stranded and seq_type == "protein":
                continue
            for circular in (False, True):
                for monoisotopic in (False, True):
                    try:
                        weight = repr(
                            molecular_weight(seq, seq_type, double_stranded, circular, monoisotopic)
                        )
                    except ValueError:
                        weight = ""
                    flags = [str(f).lower() for f in (double_stranded, circular, monoisotopic)]
                    # The weight is not last: an empty one must survive the
                    # trailing-whitespace hook.
                    rows.append((seq_type, *flags, weight, seq))
    write(
        "bio_weight.tsv",
        ["seq_type", "double_stranded", "circular", "monoisotopic", "weight", "seq"],
        rows,
    )


def write_gc_skew():
    rows = []
    for name, seq in records("corpus.fastq"):
        for window in (1, 3, 10, 100):
            rows.append((name, window, ",".join(repr(v) for v in GC_skew(seq, window))))
    write("bio_gc_skew.tsv", ["id", "window", "skew"], rows)


# (leading, trailing, window size, window quality); "" for a step not used.
TRIM_CONFIGS = [
    (10, "", "", ""),
    ("", 10, "", ""),
    ("", "", 4, 15),
    ("", "", 4, 20.5),
    ("", "", 1, 10),
    ("", "", 10, 25),
    (3, 3, 4, 15),
    (20, 20, "", ""),
    (5, 30, 3, 12),
]


def write_trimmomatic():
    jar = os.environ.get("TRIMMOMATIC_JAR")
    if not jar:
        print("trimmomatic.tsv: skipped (TRIMMOMATIC_JAR is not set)", file=sys.stderr)
        return
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.fastq"
        for leading, trailing, size, quality in TRIM_CONFIGS:
            steps = []
            if leading != "":
                steps.append(f"LEADING:{leading}")
            if trailing != "":
                steps.append(f"TRAILING:{trailing}")
            if size != "":
                steps.append(f"SLIDINGWINDOW:{size}:{quality}")
            subprocess.run(
                ["java", "-jar", jar, "SE", "-phred33", str(HERE / "corpus.fastq"), str(out), *steps],
                check=True,
                capture_output=True,
            )
            lines = out.read_text().splitlines()
            for i in range(0, len(lines), 4):
                rows.append((leading, trailing, size, quality, lines[i][1:], lines[i + 1], lines[i + 3]))
    write(
        "trimmomatic.tsv",
        ["leading", "trailing", "window_size", "window_quality", "id", "seq", "qual"],
        rows,
    )


def main():
    write_cds()
    write_translate_raw()
    write_transcribe()
    write_motif()
    write_weight()
    write_gc_skew()
    write_trimmomatic()


if __name__ == "__main__":
    main()
