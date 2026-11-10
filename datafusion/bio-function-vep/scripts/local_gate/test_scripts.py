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
        with open(os.path.join(self.arm, "contigs.tsv")) as handle:
            header = handle.readline().rstrip("\n").split("\t")
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

    def fake_tools(self, low_readings):
        """Put `df`, `ps` and `sleep` stand-ins first on PATH.

        `df` reports 5 GiB free for the first `low_readings` calls and 50 GiB afterwards,
        `ps` reports an idle machine and `sleep` returns at once.
        """
        bin_dir = os.path.join(self.tmp, "bin")
        os.makedirs(bin_dir)
        count = os.path.join(self.tmp, "df-calls")
        scripts = {
            "df": (
                "#!/bin/sh\n"
                f'n=$(cat "{count}" 2> /dev/null || echo 0)\n'
                "n=$((n + 1))\n"
                f'echo "$n" > "{count}"\n'
                f'if [ "$n" -le {low_readings} ]; then free=5; else free=50; fi\n'
                'echo "Filesystem 1G-blocks Used Available Capacity iused ifree %iused Mounted on"\n'
                'echo "/dev/fake 100 50 $free 50% 1 1 1% /"\n'
            ),
            "ps": "#!/bin/sh\necho 1.0\n",
            "sleep": "#!/bin/sh\nexit 0\n",
        }
        for name, text in scripts.items():
            path = os.path.join(bin_dir, name)
            with open(path, "w") as handle:
                handle.write(text)
            os.chmod(path, 0o755)
        return bin_dir, count

    def run_one_chromosome(self, low_readings):
        bin_dir, count = self.fake_tools(low_readings)
        env = dict(
            os.environ,
            GATE_ROOT=self.root,
            DATA_VEPYR_DIR=self.data,
            PYTHONDONTWRITEBYTECODE="1",
            PATH=bin_dir + os.pathsep + os.environ["PATH"],
        )
        done = subprocess.run(
            [os.path.join(HERE, "run_arm.sh"), self.venv, self.arm, "hg002", "8", "time", "--chroms", "chr22"],
            env=env,
            capture_output=True,
            text=True,
        )
        with open(os.path.join(self.arm, "contigs.tsv")) as handle:
            rows = handle.read().splitlines()[1:]
        with open(count) as handle:
            return done, rows, int(handle.read())

    def test_low_disk_reading_is_rechecked_before_stopping(self):
        # The space of the output deleted a moment ago comes back with a delay.
        done, rows, _ = self.run_one_chromosome(low_readings=3)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual([row.split("\t")[1] for row in rows], ["chr22"])

    def test_lasting_low_disk_stops_with_exit_3(self):
        done, rows, readings = self.run_one_chromosome(low_readings=10**6)
        self.assertEqual(done.returncode, 3, done.stderr)
        self.assertIn("STOP: 5 GiB free", done.stderr)
        self.assertEqual(rows, [])
        self.assertEqual(readings, 61)  # the first reading and two minutes of re-checks

    def test_gate_rejects_bad_targets_before_building(self):
        gate = os.path.join(self.tmp, "gate")
        for targets in ("HG02=5,ClinVar=15", "HG002=5", "HG002=5,ClinVar", "HG002=5,ClinVar=15,Other=1"):
            done = self.run_script("gate.sh", "/nonexistent-a", "/nonexistent-b", gate, "--target-pct", targets)
            self.assertEqual(done.returncode, 2, targets)
            self.assertIn("--target-pct", done.stderr)
            self.assertFalse(os.path.exists(gate), targets)


if __name__ == "__main__":
    unittest.main()
