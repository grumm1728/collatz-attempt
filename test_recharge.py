import unittest
import json
from pathlib import Path

import recharge
from collatz import shortcut_step


class RechargeTests(unittest.TestCase):
    def test_formula_against_standard_map(self):
        for m in range(1, 100):
            with self.subTest(m=m):
                value = (1 << m) - 1
                for _ in range(m):
                    value = 3 * value + 1
                    while value % 2 == 0:
                        value //= 2
                self.assertEqual(recharge.odd_after_climb(m), value)
                self.assertEqual(
                    recharge.immediately_recharges(m),
                    (value + 1) % (1 << m) == 0,
                )
        self.assertTrue(recharge.immediately_recharges(1))
        self.assertFalse(recharge.immediately_recharges(2))

    def test_small_certificate_and_independent_exhaustive_scan(self):
        bits = 8
        rows = [recharge.lift_row(bits, a) for a in range(1, bits)]
        self.assertTrue(recharge.certificate(bits, rows)["verified_no_immediate_recharge"])
        for m in range(3, 1 << bits):
            if m % 2 == 0:
                a = (m & -m).bit_length() - 1
                offset = a + 2
            else:
                offset = 1
            modulus = 1 << (m + offset)
            self.assertNotEqual(pow(3, m, modulus), (1 - (1 << offset)) % modulus)

    def test_verifier_rejects_bad_or_missing_rows(self):
        bits = 8
        rows = [recharge.lift_row(bits, a) for a in range(1, bits)]
        with self.assertRaises(ValueError):
            recharge.verify_rows(bits, rows[:-1])
        first = rows[0]
        altered = {**first, "least_b": first["least_b"] ^ 2}
        with self.assertRaises(ValueError):
            recharge.verify_rows(bits, [altered, *rows[1:]])

    def test_saved_certificate_and_separate_full_scan(self):
        saved = json.loads((Path(__file__).parent / "results" / "recharge.json").read_text())
        self.assertEqual(saved["bits"], 20)
        self.assertEqual(recharge.certificate(20, saved["residue_rows"]), saved)
        # Track only low bits. A zero residue would leave the valuation unknown
        # and must fail the check instead of being counted as an exclusion.
        mask = (1 << 43) - 1
        power = 3
        for m in range(2, 1 << 20):
            power = (3 * power) & mask
            divisor = (power - 1) & -(power - 1)
            self.assertTrue(divisor)
            numerator = (power + divisor - 1) & mask
            self.assertTrue(numerator)
            r = (numerator & -numerator).bit_length() - divisor.bit_length()
            self.assertLess(r, m)

    def test_two_prescribed_odd_runs(self):
        for m in range(2, 7):
            for q in (1, 2, 3):
                length = 2 * 3**(m - 1) * q
                cofactor, remainder = divmod((1 << length) - 1, 3**m)
                self.assertEqual((remainder, cofactor % 2), (0, 1))
                n = (cofactor << m) - 1
                for _ in range(m):
                    self.assertEqual(n % 2, 1)
                    n = shortcut_step(n)
                self.assertEqual(n, (1 << length) - 2)
                n = shortcut_step(n)
                self.assertEqual(n, (1 << (length - 1)) - 1)
                for _ in range(length - 1):
                    self.assertEqual(n % 2, 1)
                    n = shortcut_step(n)
                self.assertEqual(n % 2, 0)


if __name__ == "__main__":
    unittest.main()
