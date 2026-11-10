"""Tests for the local gate report rules.

Run: python3 -m unittest datafusion/bio-function-vep/scripts/local_gate/test_report.py -v
"""

import contextlib
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from report import COLUMNS, check, load_arm, main, totals  # noqa: E402


def rows(dataset, times, **overrides):
    """One row per contig with passing defaults; `overrides` maps column -> {contig: value}."""
    out = []
    for contig, time_s in times.items():
        row = {
            "dataset": dataset,
            "contig": contig,
            "workers": "8",
            "mode": "time",
            "compression": "plain",
            "exit": "0",
            "wall_s": f"{float(time_s) + 1.0:.2f}" if time_s != "NA" else "NA",
            "rss_gib": "4.000",
            "time_s": str(time_s),
            "records": "100",
            "body_md5": f"md5-{contig}",
            "verdict": "-",
            "env": "",
            "bg_cpu": "50",
        }
        for column, by_contig in overrides.items():
            if contig in by_contig:
                row[column] = str(by_contig[contig])
        out.append(row)
    return out


class ReportTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._n = 0

    def arm_dir(self, arm_rows):
        self._n += 1
        path = os.path.join(self._tmp.name, f"arm{self._n}")
        os.makedirs(path)
        with open(os.path.join(path, "contigs.tsv"), "w") as handle:
            handle.write("\t".join(COLUMNS) + "\n")
            for row in arm_rows:
                handle.write("\t".join(row[c] for c in COLUMNS) + "\n")
        return path

    def arm(self, arm_rows):
        return load_arm(self.arm_dir(arm_rows))

    def results(self, base, cand, parity=(), identity=(), base_1w=None, cand_1w=None, targets=None):
        checked = check(base, cand, list(parity), list(identity), base_1w, cand_1w, targets or {})
        return {name: ok for name, ok, _ in checked}

    def pair(self, base_times, cand_times, **cand_overrides):
        base = [self.arm(rows("hg002", base_times)) for _ in range(2)]
        cand = [self.arm(rows("hg002", cand_times, **cand_overrides)) for _ in range(2)]
        return base, cand

    def test_load_arm_and_totals(self):
        loaded = self.arm(rows("hg002", {"chr1": 10.0, "chr2": 8.5}) + rows("clinvar", {"1": 3.0}))
        self.assertEqual(len(loaded), 3)
        self.assertEqual(totals(loaded), {"hg002": 18.5, "clinvar": 3.0})

    def test_meets_target_and_passes(self):
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 8.0, "chr2": 7.0})
        results = self.results(base, cand, targets={"hg002": 10.0})
        self.assertTrue(results["target hg002"])
        self.assertTrue(all(results.values()), results)

    def test_target_missed(self):
        # 18.0 -> 17.5 is 2.8%, target 10%.
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 9.7, "chr2": 7.8})
        self.assertFalse(self.results(base, cand, targets={"hg002": 10.0})["target hg002"])

    def test_negative_target_allows_noise(self):
        # A target of -3 lets the candidate be up to 3% slower.
        base, cand = self.pair({"chr1": 100.0}, {"chr1": 102.0})
        self.assertTrue(self.results(base, cand, targets={"hg002": -3.0})["target hg002"])
        base, cand = self.pair({"chr1": 100.0}, {"chr1": 104.0})
        self.assertFalse(self.results(base, cand, targets={"hg002": -3.0})["target hg002"])

    def test_single_contig_regression(self):
        # chr2 8.0 -> 9.0: +12.5% and +1.0 s.
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 8.0, "chr2": 9.0})
        self.assertFalse(self.results(base, cand)["no contig regression hg002"])

    def test_small_absolute_regression_is_allowed(self):
        # chr2 1.0 -> 1.3: +30% but only 0.3 s.
        base, cand = self.pair({"chr1": 10.0, "chr2": 1.0}, {"chr1": 10.0, "chr2": 1.3})
        self.assertTrue(self.results(base, cand)["no contig regression hg002"])

    def test_rss_over_ten_percent_fails(self):
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 10.0, "chr2": 8.0}, rss_gib={"chr2": "4.500"})
        self.assertFalse(self.results(base, cand)["peak rss hg002"])
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 10.0, "chr2": 8.0}, rss_gib={"chr2": "4.300"})
        self.assertTrue(self.results(base, cand)["peak rss hg002"])

    def test_identity_mismatch_fails(self):
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 10.0, "chr2": 8.0})
        self.assertTrue(self.results(base, cand)["identity"])
        other = self.arm(rows("hg002", {"chr1": 10.0, "chr2": 8.0}, body_md5={"chr2": "different"}))
        self.assertFalse(self.results(base, cand, identity=[other])["identity"])

    def test_candidate_differs_from_baseline_fails_identity(self):
        base, cand = self.pair({"chr1": 10.0}, {"chr1": 10.0}, body_md5={"chr1": "changed"})
        self.assertFalse(self.results(base, cand)["identity"])

    def test_parity_fail_verdict_fails(self):
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 10.0, "chr2": 8.0})
        good = self.arm(rows("hg002", {"chr1": 10.0, "chr2": 8.0}, verdict={"chr1": "PASS", "chr2": "PASS"}))
        bad = self.arm(rows("hg002", {"chr1": 10.0, "chr2": 8.0}, verdict={"chr1": "PASS", "chr2": "FAIL"}))
        self.assertTrue(self.results(base, cand, parity=[good])["parity"])
        self.assertFalse(self.results(base, cand, parity=[bad])["parity"])
        # A time-mode arm ("-") given as parity evidence is not a pass.
        unchecked = self.arm(rows("hg002", {"chr1": 10.0, "chr2": 8.0}))
        self.assertFalse(self.results(base, cand, parity=[unchecked])["parity"])

    def test_failed_or_missing_row_fails(self):
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 10.0, "chr2": 8.0}, exit={"chr2": "1"})
        self.assertFalse(self.results(base, cand)["complete"])
        # A contig absent from one arm.
        base = [self.arm(rows("hg002", {"chr1": 10.0, "chr2": 8.0})) for _ in range(2)]
        cand = [self.arm(rows("hg002", {"chr1": 10.0, "chr2": 8.0})), self.arm(rows("hg002", {"chr1": 10.0}))]
        self.assertFalse(self.results(base, cand)["complete"])

    def test_missing_time_fails(self):
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 10.0, "chr2": "NA"})
        self.assertFalse(self.results(base, cand)["complete"])

    def test_one_worker_regression_over_three_percent_fails(self):
        base, cand = self.pair({"chr1": 10.0}, {"chr1": 10.0})
        base_1w = [self.arm(rows("hg002", {"chr1": 100.0}, workers={"chr1": "1"}))]
        slow = [self.arm(rows("hg002", {"chr1": 104.0}, workers={"chr1": "1"}))]
        fine = [self.arm(rows("hg002", {"chr1": 102.0}, workers={"chr1": "1"}))]
        self.assertFalse(self.results(base, cand, base_1w=base_1w, cand_1w=slow)["1 worker hg002"])
        self.assertTrue(self.results(base, cand, base_1w=base_1w, cand_1w=fine)["1 worker hg002"])

    def test_datasets_are_judged_separately(self):
        base = [self.arm(rows("hg002", {"chr1": 10.0})), self.arm(rows("clinvar", {"1": 50.0}))]
        cand = [self.arm(rows("hg002", {"chr1": 9.0})), self.arm(rows("clinvar", {"1": 49.0}))]
        results = self.results(base, cand, targets={"hg002": 5.0, "clinvar": 15.0})
        self.assertTrue(results["target hg002"])
        self.assertFalse(results["target clinvar"])

    def test_last_row_wins_after_a_rerun(self):
        failed = rows("hg002", {"chr1": "NA"}, exit={"chr1": "137"}, body_md5={"chr1": "NA"})
        loaded = self.arm(failed + rows("hg002", {"chr1": 10.0}) + rows("hg002", {"chr2": 8.0}))
        self.assertEqual([(r["contig"], r["exit"], r["time_s"]) for r in loaded], [("chr1", "0", "10.0"), ("chr2", "0", "8.0")])
        base = [self.arm(rows("hg002", {"chr1": 10.0, "chr2": 8.0}))]
        self.assertTrue(self.results(base, [loaded])["complete"])

    def test_missing_target_fails(self):
        base, cand = self.pair({"chr1": 10.0}, {"chr1": 5.0})
        self.assertFalse(self.results(base, cand, targets={})["target hg002"])

    def test_empty_arm_fails_complete(self):
        base, cand = self.pair({"chr1": 10.0}, {"chr1": 10.0})
        self.assertFalse(self.results(base, cand, parity=[self.arm([])])["complete"])

    def test_group_missing_a_dataset_fails(self):
        base = [self.arm(rows("hg002", {"chr1": 10.0})), self.arm(rows("clinvar", {"1": 50.0}))]
        cand = [self.arm(rows("hg002", {"chr1": 10.0})), self.arm(rows("clinvar", {"1": 50.0}))]
        parity = [self.arm(rows("hg002", {"chr1": 10.0}, verdict={"chr1": "PASS"}))]
        self.assertFalse(self.results(base, cand, parity=parity)["complete"])

    def test_one_worker_arm_with_fewer_contigs_fails(self):
        base, cand = self.pair({"chr1": 10.0, "chr2": 8.0}, {"chr1": 10.0, "chr2": 8.0})
        base_1w = [self.arm(rows("hg002", {"chr1": 100.0, "chr2": 80.0}))]
        cand_1w = [self.arm(rows("hg002", {"chr1": 100.0}))]
        self.assertFalse(self.results(base, cand, base_1w=base_1w, cand_1w=cand_1w)["complete"])

    def test_busy_machine_fails_quiet(self):
        base, cand = self.pair({"chr1": 10.0}, {"chr1": 10.0})
        self.assertTrue(self.results(base, cand)["quiet"])
        base, cand = self.pair({"chr1": 10.0}, {"chr1": 10.0}, bg_cpu={"chr1": "350"})
        self.assertFalse(self.results(base, cand)["quiet"])
        base, cand = self.pair({"chr1": 10.0}, {"chr1": 10.0}, bg_cpu={"chr1": "NA"})
        self.assertFalse(self.results(base, cand)["quiet"])

    def test_cli_rejects_an_unknown_target_name(self):
        base = self.arm_dir(rows("hg002", {"chr1": 10.0}))
        cand = self.arm_dir(rows("hg002", {"chr1": 8.0}))
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(["--baseline", base, "--candidate", cand, "--target-pct", "HG02=10"]), 2)

    def test_cli_exit_codes_and_target_key_case(self):
        base = self.arm_dir(rows("hg002", {"chr1": 10.0}))
        cand = self.arm_dir(rows("hg002", {"chr1": 8.0}))
        slow = self.arm_dir(rows("hg002", {"chr1": 9.9}))
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(main(["--baseline", base, "--candidate", cand, "--target-pct", "HG002=10"]), 0)
            self.assertIn("VERDICT: PASS", out.getvalue())
            self.assertEqual(main(["--baseline", base, "--candidate", slow, "--target-pct", "HG002=10"]), 1)
            self.assertIn("VERDICT: FAIL", out.getvalue())
        with contextlib.redirect_stderr(io.StringIO()):
            missing = os.path.join(self._tmp.name, "nowhere")
            self.assertEqual(main(["--baseline", missing, "--candidate", cand, "--target-pct", "HG002=10"]), 2)
            self.assertEqual(main(["--baseline", base, "--candidate", cand, "--target-pct", "HG002"]), 2)


if __name__ == "__main__":
    unittest.main()
