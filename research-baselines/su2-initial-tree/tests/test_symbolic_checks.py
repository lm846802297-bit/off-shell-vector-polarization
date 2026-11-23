import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class SymbolicChecks(unittest.TestCase):
    def run_script(self, script, result_file):
        subprocess.run([sys.executable, str(ROOT / script)], check=True, cwd=ROOT)
        data = json.loads((ROOT / result_file).read_text(encoding="utf-8"))
        return data

    def test_projector_identities(self):
        data = self.run_script(
            "scripts/verify_projector_identities.py",
            "results/projector_identities.json",
        )
        self.assertTrue(data["all_pass"])
        self.assertTrue(all(data["checks"].values()))

    def test_numeric_projector_scan(self):
        data = self.run_script(
            "scripts/numeric_projector_scan.py",
            "results/projector_numeric_scan.json",
        )
        self.assertTrue(data["all_pass"])
        self.assertEqual(data["samples"], 10000)

    def test_offshell_ambiguity(self):
        data = self.run_script(
            "scripts/verify_offshell_ambiguity.py",
            "results/off_shell_ambiguity.json",
        )
        self.assertTrue(data["all_pass"])

    def test_vertex_st_identity(self):
        data = self.run_script(
            "scripts/verify_vertex_st_identity.py",
            "results/vertex_st_identity.json",
        )
        self.assertTrue(data["pass"])
        self.assertEqual(data["residual_after_assumptions"], "0")


if __name__ == "__main__":
    unittest.main()
