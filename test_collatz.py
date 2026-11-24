import unittest

from collatz import odd_step, shortcut_step, standard_step


class ReferenceChecks(unittest.TestCase):
    def test_trajectories_through_10000(self):
        for start in range(1, 10_001):
            n = start
            while n % 2 == 0:
                self.assertEqual(shortcut_step(n), standard_step(n))
                n //= 2
            while n != 1:
                expected = standard_step(n)
                next_short = standard_step(expected)
                self.assertEqual(shortcut_step(n), next_short)
                halvings = 0
                while expected % 2 == 0:
                    expected = standard_step(expected)
                    halvings += 1
                self.assertEqual(odd_step(n), (expected, halvings))
                n = expected

    def test_standard_cycle(self):
        self.assertEqual([standard_step(n) for n in (1, 4, 2)], [4, 2, 1])
        self.assertEqual(odd_step(1), (1, 2))

    def test_five_cycle(self):
        cycle = [13, 66, 33, 166, 83, 416, 208, 104, 52, 26]
        self.assertEqual([standard_step(n, 5) for n in cycle], cycle[1:] + cycle[:1])
        self.assertEqual([odd_step(n, 5) for n in (13, 33, 83)],
                         [(33, 1), (83, 1), (13, 5)])

    def test_even_multiplier_escape(self):
        # If n = 2**j - 1 then 2*n+1 = 2**(j+1) - 1, again odd.
        for j in range(2, 150):
            n = 2**j - 1
            self.assertEqual(standard_step(n, 2), 2**(j + 1) - 1)
            self.assertEqual(shortcut_step(n, 2), standard_step(n, 2))
            self.assertEqual(odd_step(n, 2), (standard_step(n, 2), 0))

    def test_invalid_domain(self):
        for step in (standard_step, shortcut_step, odd_step):
            with self.assertRaises(ValueError):
                step(0)
            with self.assertRaises(ValueError):
                step(1, 0)
        with self.assertRaises(ValueError):
            odd_step(4)


if __name__ == "__main__":
    unittest.main()
