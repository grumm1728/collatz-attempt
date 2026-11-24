# Collatz counterexample search

We seek a positive integer whose standard Collatz trajectory avoids 1 forever, or a proof that such an integer exists. An explicit positive cycle excluding 1 has a finite certificate. An escape proof could instead give an increasing, unbounded sequence of checkpoints along one fixed orbit.

All 18 structured escape candidates from the pilot now have verified paths to 1. A follow-up argument rules out an immediately equal or longer second climb for every Mersenne seed `2^m−1` with `m≥2`. Varying the odd cofactor produces explicit seeds with any prescribed first climb length and arbitrarily long second climbs. [Recharge proof and construction](research/recharge.md) · [Completed trajectories](research/continuation.md) · [Figure collection](figures/README.md)

The first pilot also eliminated a bounded family of compressed cycle recipes. Its original report preserves the initial stopping rules and findings. [First sprint](research/pilot.md)

## Running the research code

Use Python 3.13 from the repository root. The arithmetic and reference checks use the Python standard library.

```powershell
python -m unittest -v test_collatz
```

For the multiplier experiments the standard step halves an even integer and maps an odd integer to `k*n+1`. With odd `k`, the shortcut combines an odd step and the first halving. An odd-to-odd step combines an odd step and *all* following halvings. Each result identifies its counting convention.

The reference suite checks all starting values from 1 through 10,000 against standard arithmetic. Positive controls include the standard cycle through 1, the cycle through 13 at multiplier 5, and the identity `2*(2**j-1)+1 = 2**(j+1)-1` for the odd branch at multiplier 2. The last identity proves escape from 3 for that variant by induction.

## First pilot

Two local rounds explore compressed cycle recipes, forced-growth starts and independent wildcard ideas. The numerical escape budget is 20,000 shortcut steps per trajectory. Each arm has a 15 CPU-minute ceiling per round. Fixed-seed controls and independent checks accompany the results. Budget exhaustion records an unresolved trajectory.

The published verification below `2**71` and the resulting cycle-length bounds guide the scale of the cycle arm. [Bařina](https://link.springer.com/article/10.1007/s11227-025-07337-0) · [Hercher and corrigendum](https://cs.uwaterloo.ca/journals/JIS/VOL26/Hercher/hercher5.html)
