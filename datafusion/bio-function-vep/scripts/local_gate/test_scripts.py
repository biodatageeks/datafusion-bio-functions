"""Tests of the shell scripts' guards, run against fake venv, vepyr and data directories.

Run: python3 datafusion/bio-function-vep/scripts/local_gate/test_scripts.py -v
No build and no annotation is started by these tests.
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


class ScriptTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = self._tmp.name
        self.root = os.path.join(self.tmp, "gate-root")
        os.makedirs(os.path.join(self.root, "vepyr-build", "e2e-testing", "scripts"))
        open(os.path.join(self.root, "vepyr-build", "e2e-testing", "scripts", "run_comparison.py"), "w").close()
        self.data = os.path.join(self.tmp, "data")
        os.makedirs(os.path.join(self.data, "cache", "116_GRCh38_merged"))
        self.venv = os.path.join(self.tmp, "venv")
        os.makedirs(os.path.join(self.venv, "bin"))
        os.symlink(sys.executable, os.path.join(self.venv, "bin", "python"))
        self.identity("aaa")
        self.arm = os.path.join(self.tmp, "arm")

    def identity(self, digest):
        with open(os.path.join(self.venv, "GATE_IDENTITY.json"), "w") as handle:
            json.dump({"functions_sha": "f" * 40, "functions_dirty": False, "extension_sha256": digest}, handle)

    def run_script(self, name, *args):
        env = dict(os.environ, GATE_ROOT=self.root, DATA_VEPYR_DIR=self.data, PYTHONDONTWRITEBYTECODE="1")
        return subprocess.run([os.path.join(HERE, name), *args], env=env, capture_output=True, text=True)

    def run_arm(self, *args):
        # A chromosome list of one space runs no chromosome: only the set-up is exercised.
        return self.run_script("run_arm.sh", self.venv, self.arm, *args, "--chroms", " ")

    def test_header_has_the_background_cpu_column(self):
        done = self.run_arm("hg002", "8", "time")
        self.assertEqual(done.returncode, 0, done.stderr)
        header = open(os.path.join(self.arm, "contigs.tsv")).readline().rstrip("\n").split("\t")
        self.assertEqual(header[-1], "bg_cpu")
        self.assertEqual(len(header), 14)

    def test_same_build_and_arguments_resume(self):
        self.assertEqual(self.run_arm("hg002", "8", "time").returncode, 0)
        self.assertEqual(self.run_arm("hg002", "8", "time").returncode, 0)

    def test_refuses_an_arm_of_another_build(self):
        self.assertEqual(self.run_arm("hg002", "8", "time").returncode, 0)
        self.identity("bbb")
        done = self.run_arm("hg002", "8", "time")
        self.assertEqual(done.returncode, 2)
        self.assertIn("another build or other arguments", done.stderr)

    def test_refuses_other_arguments(self):
        self.assertEqual(self.run_arm("hg002", "8", "time").returncode, 0)
        for args in (("hg002", "4", "time"), ("hg002", "8", "parity"), ("hg002", "8", "time", "--bgzf"),
                     ("hg002", "8", "time", "VEP_X=1")):
            self.assertEqual(self.run_arm(*args).returncode, 2, args)

    def test_gate_rejects_bad_targets_before_building(self):
        gate = os.path.join(self.tmp, "gate")
        for targets in ("HG02=5,ClinVar=15", "HG002=5", "HG002=5,ClinVar", "HG002=5,ClinVar=15,Other=1"):
            done = self.run_script("gate.sh", "/nonexistent-a", "/nonexistent-b", gate, "--target-pct", targets)
            self.assertEqual(done.returncode, 2, targets)
            self.assertIn("--target-pct", done.stderr)
            self.assertFalse(os.path.exists(gate), targets)


if __name__ == "__main__":
    unittest.main()
