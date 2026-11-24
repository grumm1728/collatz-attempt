from fractions import Fraction
from itertools import permutations
from math import gcd
from random import Random
import unittest

import cycles


class CycleSearchTests(unittest.TestCase):
    def test_modular_affine_against_standard_map(self):
        prime = 31
        for motif_index, motif in enumerate(cycles.MOTIFS):
            schedule = [(motif_index, 1)]
            with self.subTest(motif=motif):
                checked = 0
                a, b, c = cycles.modular_affine(schedule, prime)
                residues = cycles.residue_map(schedule, prime)
                for start in range(1, 1000, 2):
                    end = start
                    observed = []
                    for _ in motif:
                        end, divisions = cycles.standard_odd_step(end)
                        observed.append(divisions)
                    if tuple(observed) == motif:
                        self.assertEqual((a * start + c) % prime, b * end % prime)
                        self.assertEqual(residues[start % prime], end % prime)
                        checked += 1
                self.assertGreater(checked, 0)

    def test_quotient_jumping_against_exhaustive_small_window(self):
        lower, upper = (7, 5), (8, 5)
        cap = 100
        expected = {
            (numerator, denominator)
            for denominator in range(1, cap + 1)
            for numerator in range(1, 2 * denominator)
            if 7 * denominator < 5 * numerator < 8 * denominator
            and gcd(numerator, denominator) == 1
        }
        observed = list(cycles.bounded_slopes(lower, upper, cap))
        self.assertEqual(len(observed), len(expected))
        self.assertEqual(set(observed), expected)

    def test_shuffle_retains_packet_multiset_and_lengths(self):
        counts = (1234, 5678)
        alternating, shuffled = cycles.matched_schedules(counts, Random(42))
        self.assertCountEqual(alternating, shuffled)
        self.assertEqual(
            tuple(sum(copies for index, copies in alternating if index == motif) for motif in range(2)),
            counts,
        )
        self.assertEqual(cycles.motif_counts(*cycles.lengths(counts)), counts)

    def test_packet_minimum_bound_against_exact_fraction_closing_values(self):
        tiny_packets = ((0, 1), (0, 2), (1, 1), (1, 3))
        caps = cycles.inverse_packet_weight_caps()
        self.assertEqual(caps, (Fraction(5), Fraction(1)))
        certificate = cycles.packet_minimum_certificate([(10, 17)], 1_000)
        self.assertTrue(certificate["all_covered_packet_words_excluded"])

        for packets in permutations(tiny_packets):
            with self.subTest(packets=packets):
                expanded = [cycles.MOTIFS[index] * copies for index, copies in packets]
                odd_steps = sum(map(len, expanded))
                divisions = sum(map(sum, expanded))
                self.assertEqual((odd_steps, divisions), (10, 17))

                scores = [0]
                for packet in expanded[:-1]:
                    scores.append(scores[-1] + odd_steps * sum(packet) - divisions * len(packet))
                start = scores.index(max(scores))
                rotation = expanded[start:] + expanded[:start]
                multiplier = Fraction(2**divisions, 3**odd_steps)
                prefix_weight, offset = Fraction(1), Fraction(0)
                for packet in rotation:
                    packet_start = prefix_weight
                    packet_offset = Fraction(0)
                    for exponent in packet:
                        packet_offset += prefix_weight / 3
                        prefix_weight *= Fraction(2**exponent, 3)
                    self.assertLessEqual(prefix_weight, multiplier)
                    motif_index = 0 if packet[0] == 1 else 1
                    self.assertLess(packet_offset, caps[motif_index] * max(packet_start, prefix_weight))
                    offset += packet_offset

                closing_value = offset / (multiplier - 1)
                self.assertLess(closing_value, 12 / (1 - 1 / multiplier))
                self.assertLess(closing_value, certificate["strict_upper_bound_on_minimum"])


if __name__ == "__main__":
    unittest.main()
