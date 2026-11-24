# Compressed cycle search

An entry `a` denotes the odd step followed by exactly `a` halvings. Write `U=(1,2)` and `V=(2)`. A packet is a positive number of repetitions of just one of these motifs. A closed recipe with `K` odd visits and `S` halvings uses `K+S` standard steps. These are counts with multiplicity in the proposed word.

## Reproduce the two rounds

```powershell
python -m unittest -v test_collatz test_cycles
python cycles.py --round all --output results/cycles.json
```

The code fixes the PRNG seeds and records the full parameters and rejection counts in [the machine-readable report](../results/cycles.json). The model sent each prediction to the lead before running its round. The repository records the predictions alongside the results.

1. Start with 59,813,072,756 copies of U followed by 24,488,718,470 copies of V. The composed residue map is `x → x+7` modulo 31, which precludes closure. Rearrange the packets at the same totals. The split search exhausts all 90 modular exponent classes of the two-split recipe. The sampled search takes 256 matched packet multisets and tests each in alternating order and shuffled order.
2. Use the verified minimum `B=2^71` and the pilot interval `137,500,000,000 ≤ K ≤ 10^13`. Multiplication of `2^a x_next = 3x+1` around a hypothetical cycle gives `δ = S ln 2 − K ln 3 = Σ ln(1+1/(3x))`, so `0 < δ < K/(3B)`. The generator enumerates all 10,181 integer length pairs in that necessary window using bounded logarithms and rational arithmetic. Of the 1,077 pairs admitting both a discovery prime and the held-out prime 7 in the closing denominator, test the first 128 ordered by K then S. Generate eight packet multisets per pair. The modular sieve compares an affine composition with an independent permutation on every residue.

| Trials | Alternating order | Shuffled order |
|---|---:|---:|
| Round 1 rejected by 31, out of 256 | 251 | 247 |
| Round 2 rejected by discovery primes, out of 1,024 | 979 | 982 |
| Further rejections by 7 | 34 of 45 | 34 of 42 |
| Passed the modular tests | 11 | 8 |

Both pre-run predictions held. Round one had at least 90% rejection for both orders with a gap under ten percentage points. Round two's held-out prime rejected at least half the survivors of each ordering. Eligibility for that prime was selected in advance by the length pair. These rates describe the specified samples.

## Why the 19 survivors also fail

The following argument was developed **after** the modular runs. It covers every cyclic ordering of at most 16 U packets and at most 16 V packets in the stated length window, with arbitrary positive packet sizes. It uses the published verification of all starts below `B`. [Bařina](https://link.springer.com/article/10.1007/s11227-025-07337-0)

For any prefix with `j` odd steps and `s` halvings, its inverse affine multiplier is `W=2^s/3^j`. Rotate the recipe to a packet boundary maximizing `s−jS/K` over the original boundaries. The centered sum over the whole cycle is zero. After this rotation every boundary satisfies `s−jS/K ≤ 0`, even across the original seam. Hence `W ≤ exp(δj/K) ≤ exp(δ)`.

The inverse of one U is `x_start=(8/9)x_end−5/9`. For V it is `x_start=(4/3)x_end−1/3`. Let `Q` be the total subtracted offset in the composition. Summing the geometric series within a U packet bounds its contribution by `5 W_start`. A V packet contributes exactly `W_end−W_start`. Each boundary weight is at most `exp(δ)`, so `Q ≤ (5·16+16)exp(δ)`.

Closure requires `x_start=exp(δ)x_start−Q`. Therefore

```text
x_start ≤ 96/(1−exp(−δ)) ≤ 96(1+1/δ).
```

The last inequality follows from `exp(δ) ≥ 1+δ`. Over the 10,181 possible length pairs, the smallest certified lower bound for δ exceeds `6.59 × 10^-14`. Using the unrounded rational bound gives an upper bound of **1,456,116,810,424,474** on the odd value at the chosen boundary. It is below `2^71`. Every positive integer there reaches 1, which contradicts membership in a cycle excluding 1. This also explains why finding better modular survivors within the same packet budget cannot rescue this construction.

The proof uses the packet count and the bounded range of K. It leaves open more complex words, periods above the pilot ceiling, and trajectories that never repeat.
