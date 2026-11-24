# First counterexample pilot

**Follow-up.** The three starts censored by this pilot have since reached 1. The immediate-recharge question posed below has also been resolved for every exponent. This note preserves what was known at the original cutoff. [Completed paths](continuation.md) · [Recharge proof](recharge.md)

The two rounds found no positive Collatz counterexample. The strongest result rules out a whole family of cycle recipes in the pilot's length range. The simplest unfinished trajectory is `2^4096−1`. It fell from a 4,096-bit seed to a 3,089-bit integer by the step limit, so it cannot be the *least* counterexample. Its eventual destination was unresolved within the pilot.

## Which cycles can we rule out?

Use the halving motifs `U=(1,2)` and `V=(2)`. A short description can specify billions of copies of either motif. In 2,048 second-round trials, 19 descriptions passed every prime we tested. A subsequent proof rules out all cyclic arrangements of up to 16 runs of each motif, for periods between 137,500,000,000 and 10,000,000,000,000 odd visits. At a suitable boundary, the starting value would have to be below 1,456,116,810,424,474, inside the published convergence range. The proof covers the choices of packet lengths and order which the random samples missed. Its use of an externally verified bound is explicit. [Proof, assumptions and exact certificate](cycle-search.md)

## What happens after forced growth?

The escape arm fixed 18 starts of the form `u·2^(2^j)−1` for 3x+1 with `u∈{1,3,5}` and `7≤j≤12`. For each odd multiplier `k=3,5,7`, set `m=2^j`, `d=k−2`, and take the `(u+1)/2`-th positive odd `v` satisfying `2^m·v ≡ 1 (mod d)`. Start at `(2^m·v−1)/d`. An odd shortcut step multiplies `d·x+1` by `k/2`. This proves that the first `m` shortcut steps follow the odd branch and end at `(k^m·v−1)/d`.

The prediction was that most structured starts at 3 would reach 1, with larger post-climb gains and more cap outcomes at 5 and 7. Round one matched random odd controls by seed bit length. Round two followed the *same* structured paths and drew new controls matched to the first odd value after the forced climb. The table compares an equal horizon of 128 shortcut steps from that odd value. All 18 observations in each group reached that horizon.

| Multiplier | Structured starts reaching 1 / 18 | Median change in log₂(value), structured | Median change in log₂(value), round-two controls |
|---|---:|---:|---:|
| 3 | 15 | −25.77 | −24.18 |
| 5 | 0 | +21.76 | +20.60 |
| 7 | 0 | +51.67 | +44.65 |

The three remaining 3x+1 starts all had `j=12`. Each had dropped below its own seed at the 20,000-shortcut-step cap. All 36 structured traces at 5 and 7 reached the cap. Raw stopping counts use a common budget from the seed, so a longer forced prefix leaves fewer steps to observe its continuation. The 216 saved rows represent 162 distinct multiplier/seed pairs. An independent single-standard-step implementation reproduced their endpoints, peaks and step counts. [Data and stopping rules](../results/escape.json) · [Illustrated comparison](../figures/README.md)

For each of the 18 canonical starts, we checked whether it visited any of the three next-generation targets `u′·2^(2^(j+1))−1`. None of the 54 start/target comparisons hit within the observation window. For `u′=1` there is an exact obstruction. The target is divisible by 3, but after the first odd Collatz step the value is nonzero modulo 3, and halving preserves that property.

## Two wildcard tests

- **Can an initial climb immediately recharge itself?** From `2^m−1`, after `m` odd-to-odd moves the odd endpoint is the odd part of `3^m−1`. Write it as `y`. The agent predicted `v₂(y+1)<m` throughout `3≤m≤20,000`, where `v₂` counts factors of 2. Exact evaluation confirmed the prediction. The largest measured valuation was 17, first at `m=1198`. For every odd `m`, the valuation equals 1 since `3^m≡3 (mod 8)`. A valuation of 17 measures the low bits that can force a further climb.
- **Can residue classes form a closed set avoiding 1?** The reverse search reached every residue in the conservative quotient graph for 4,103 moduli, covering 9,967,359 residue classes. A quotient path can use different integer representatives at successive steps. A literature check gives the stronger conclusion for this idea: Monks proved every nonconstant arithmetic progression contains a start whose orbit merges with the orbit of 1. Thus no nonempty union of whole residue classes can consist entirely of counterexamples. [Exact wildcard output](../results/wildcard.json) · [Monks, Theorem 1.1](https://monks.scranton.edu/content/pubs/SufficiencyRev4.pdf)

## The next mathematical question

**For every even `m≥4`, is the recharge valuation strictly less than `m`?** Define `t=v₂(3^m−1)`. A failure is equivalent to the exact congruence `3^m ≡ 1−2^t (mod 2^(m+t))`. Split by `v₂(m)`, which fixes `t=2+v₂(m)` for even `m`, then seek a bound or lift a candidate residue to an actual integer satisfying the full congruence. A proof would eliminate immediate full recharge throughout this family. A failing integer would give a short, reproducible seed for testing a second recharge. Establishing infinitely many successive recharges along one orbit would still require a separate argument.

## Reproduction and model roles

Astra directed the rounds and developed the packet-boundary exclusion with the cycle agent. GPT-6.2 Sol implemented the cycle search. An independent GPT-6.2 Astra session handled escape identities and controls. GPT-6.1 Sol proposed the wildcard tests, and a fresh Sol session audited the arithmetic and the exclusion argument. Predictions were communicated in the working chat before the relevant runs. The packet exclusion was a follow-up to the observed modular survivors. The agent comparisons here are descriptive.

From the repository root, run the following with Python 3.13 and Matplotlib 3.10.8. The 15 CPU-minute allowance applied per arm per round. All default runs completed within it.

```powershell
python cycles.py --round all --output results/cycles.json
python escape.py
python wildcard.py --round 1 --max-k 20000
python wildcard.py --round 2 --max-modulus 4096
python figures.py
python -m unittest -v
```

Merged implementation milestones are [foundation, 85cd9db](https://github.com/grumm1728/collatz-attempt/pull/1), [cycle search, 177de28](https://github.com/grumm1728/collatz-attempt/pull/2), and [escape experiments and figures, 9ac845b](https://github.com/grumm1728/collatz-attempt/pull/3). The findings pull request links its final merge commit in the GitHub history.
