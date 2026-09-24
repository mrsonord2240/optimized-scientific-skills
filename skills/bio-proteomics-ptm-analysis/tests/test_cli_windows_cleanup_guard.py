"""Source-contract test for the opt-in Windows cli teardown workaround."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
GUARD = ROOT / "scripts" / "cli_windows_cleanup_guard.Rprofile"


class CliWindowsCleanupGuardContractTest(unittest.TestCase):
    def test_is_a_narrow_onload_guard_without_computation_changes(self):
        text = GUARD.read_text(encoding="utf-8")

        self.assertIn('setHook(packageEvent("cli", "onLoad")', text)
        self.assertIn('.Platform$OS.type != "windows"', text)
        self.assertIn("cli_state$unloaded <- TRUE", text)
        self.assertIn('Sys.getenv("CLI_WINDOWS_CLEANUP_GUARD_RECEIPT"', text)
        self.assertIn('file(receipt, open = "wx")', text)
        self.assertIn('"cli_unloaded=TRUE"', text)
        self.assertNotIn(".Last", text)
        self.assertNotIn("library(", text)
        self.assertNotIn("MSstats", text)
        self.assertNotIn("groupComparisonPTM", text)
        self.assertNotIn("write.csv", text)
        self.assertNotIn("system(", text)


if __name__ == "__main__":
    unittest.main()
