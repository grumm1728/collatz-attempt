"""Replay the saved observations using single standard-map transitions."""

from functools import lru_cache
import json
from pathlib import Path
import unittest

from collatz import standard_step
import escape


@lru_cache(maxsize=216)
def standard_replay(k, seed_hex, standard_budget):
    n = int(seed_hex, 16)
    peak = n
    shortcut_steps = 0
    steps = 0
    while steps < standard_budget and n != 1:
        # For these odd multipliers each shortcut contains exactly one halving.
        if n % 2 == 0:
            shortcut_steps += 1
        n = standard_step(n, k)
        peak = max(peak, n)
        steps += 1
    return n, peak, steps, shortcut_steps


class EscapeChecks(unittest.TestCase):
    def test_controls_and_prefixes(self):
        escape.check()

    def test_saved_endpoints_and_standard_peaks(self):
        path = Path(__file__).resolve().parent / "results" / "escape.json"
        rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
        self.assertEqual(len(rows), 216)
        for row in rows:
            with self.subTest(k=row["k"], j=row["j"], u=row["u"],
                              kind=row["kind"], round=row["round"]):
                end, peak, steps, shortcuts = standard_replay(
                    row["k"], row["seed_hex"], row["standard_steps"])
                self.assertEqual(hex(end), row["terminal_hex"])
                self.assertEqual(hex(peak), row["peak_hex"])
                self.assertEqual((steps, shortcuts), (row["standard_steps"], row["shortcut_steps"]))
                self.assertEqual(end == 1, row["status"] == "reached_one")
                if row["status"] == "cap":
                    self.assertEqual(shortcuts, escape.CAP)


if __name__ == "__main__":
    unittest.main()
