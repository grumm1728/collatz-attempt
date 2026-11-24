"""Reproducible exact-integer wildcard checks for the local Collatz pilot.

U(n)=(3*n+1)/2**v2(3*n+1), for positive odd n. One U step is
one standard odd step followed by v2(3*n+1) standard even steps.
No finite trajectory length or size is interpreted as nontermination.
"""

import argparse
from collections import Counter, deque
import json
from math import gcd
from pathlib import Path
import sys
import time


def v2(n):
    assert n > 0
    return (n & -n).bit_length() - 1


def standard_step(n):
    return 3 * n + 1 if n & 1 else n // 2


def standard_odd_step(n):
    """Independent standard-map implementation, including standard step count."""
    assert n > 0 and n & 1
    n = standard_step(n)
    steps = 1
    while not n & 1:
        n = standard_step(n)
        steps += 1
    return n, steps


def verify():
    n = 1
    cycle = [n]
    for _ in range(3):
        n = standard_step(n)
        cycle.append(n)
    assert cycle == [1, 4, 2, 1]
    for k in range(1, 129):
        n = (1 << k) - 1
        for j in range(k):
            n, steps = standard_odd_step(n)
            if j < k - 1:
                assert steps == 2
                assert n == 3 ** (j + 1) * 2 ** (k - j - 1) - 1
        p = 3 ** k - 1
        assert n == p >> v2(p)
    return {"standard_cycle": cycle, "terminal_formula_checked_k": [1, 128]}


def round1(max_k):
    power = 1
    histogram = Counter()
    records = []
    hits = []
    odd_failures = []
    best = 0
    for k in range(1, max_k + 1):
        power *= 3
        terminal = (power - 1) >> v2(power - 1)
        recharge = v2(terminal + 1)
        histogram[recharge] += 1
        if k >= 3 and recharge >= k:
            hits.append({"k": k, "recharge": recharge})
        if k & 1 and recharge != 1:
            odd_failures.append(k)
        if recharge > best:
            best = recharge
            records.append({"k": k, "recharge": recharge})
    return {
        "hypothesis": "Some k>=3 recharges at least k: v2(U^k(2^k-1)+1)>=k",
        "prediction_before_test": "No hits through k=20000; odd k always recharge 1",
        "range": {"k_min": 3, "k_max": max_k},
        "stop_rule": "Exhaust exact k bound; do not claim general impossibility",
        "hits": hits,
        "odd_k_formula_failures": odd_failures,
        "recharge_histogram_including_k1_k2": dict(sorted(histogram.items())),
        "maximum_recharge": best,
        "record_recharges": records,
        "proof": [
            "For j<k, U^j(2^k-1)=3^j*2^(k-j)-1. The first k-1 U steps divide by exactly 2.",
            "At j=k, the odd terminal is (3^k-1)/2^v2(3^k-1).",
            "For odd k, 3^k=3 mod 8, so the terminal is 1 mod 4, and v2(terminal+1)=1.",
            "Finite failure for even k does not establish failure for all k, and a recharge hit alone would not establish escape.",
        ],
    }


def modular_predecessors(y, modulus):
    """All residue predecessors in the conservative standard-map quotient.

    Each residue denotes every positive integer in that class. For odd
    moduli both parities occur in a class. For even moduli parity is fixed.
    Dividing the even integers always gives predecessor 2*y modulo m.
    The odd case solves 3*r+1=y modulo m, retaining odd r if m is even.
    """
    yield (2 * y) % modulus
    divisor = gcd(3, modulus)
    if (y - 1) % divisor:
        return
    reduced = modulus // divisor
    base = ((y - 1) // divisor * pow(3 // divisor, -1, reduced)) % reduced
    for lift in range(divisor):
        predecessor = base + lift * reduced
        if modulus & 1 or predecessor & 1:
            yield predecessor


def modular_reachability(modulus, explicit=False):
    if explicit:
        incoming = [[] for _ in range(modulus)]
        for r in range(modulus):
            for lift in range(2):
                representative = r + lift * modulus
                incoming[standard_step(representative) % modulus].append(r)
        predecessors = incoming.__getitem__
    else:
        predecessors = lambda y: modular_predecessors(y, modulus)
    distance = [-1] * modulus
    distance[1 % modulus] = 0
    queue = deque([1 % modulus])
    while queue:
        y = queue.popleft()
        for r in predecessors(y):
            if distance[r] < 0:
                distance[r] = distance[y] + 1
                queue.append(r)
    return distance


def round2(max_modulus):
    for modulus in range(2, 129):
        assert modular_reachability(modulus) == modular_reachability(modulus, explicit=True)
    selected = sorted({2 ** a * 3 ** b for a in (0, 2, 4, 8, 12)
                       for b in range(6) if 2 <= 2 ** a * 3 ** b <= 1000000})
    moduli = sorted(set(range(2, max_modulus + 1)) | set(selected))
    failures = []
    max_distance = -1
    max_distance_moduli = []
    selected_results = []
    for modulus in moduli:
        distance = modular_reachability(modulus)
        unreachable = [r for r, value in enumerate(distance) if value < 0]
        if unreachable:
            failures.append({"modulus": modulus, "unreachable_count": len(unreachable),
                             "first_unreachable": unreachable[:100]})
        longest = max(distance)
        if longest > max_distance:
            max_distance = longest
            max_distance_moduli = [modulus]
        elif longest == max_distance:
            max_distance_moduli.append(modulus)
        if modulus in selected:
            selected_results.append({"modulus": modulus, "reachable": modulus - len(unreachable),
                                     "maximum_quotient_distance": longest})
    return {
        "hypothesis": "There is a nonempty periodic forward-closed set excluding residue 1, for some tested modulus",
        "prediction_before_test": "Every residue reaches residue 1 in the conservative quotient graph",
        "all_integer_moduli_range": [2, max_modulus],
        "number_of_moduli": len(moduli),
        "total_residue_classes": sum(moduli),
        "largest_modulus": max(moduli),
        "selected_results": selected_results,
        "unreachable_failures": failures,
        "maximum_quotient_distance": max_distance,
        "moduli_achieving_maximum": max_distance_moduli,
        "independent_graph_verification_moduli": [2, 128],
        "stop_rule": "Exhaust exact modulus list; investigate any unreachable residue",
        "interpretation": [
            "Any actual standard trajectory induces a path in this quotient graph.",
            "If any residue were unreachable from residue 1 in reverse, the unreachable residue classes would be forward closed and exclude 1, giving a genuine nonconvergence certificate.",
            "All residues reachable in this conservative graph excludes this periodic invariant mechanism only for tested moduli.",
            "A quotient path need not lift to one actual trajectory. Reachability here does not prove that integers reach 1, nor does a distance bound constrain actual stopping times.",
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--round", type=int, choices=[1, 2], default=1)
    parser.add_argument("--max-k", type=int, default=20000)
    parser.add_argument("--max-modulus", type=int, default=4096)
    parser.add_argument("--output", default="results/wildcard.json")
    args = parser.parse_args()
    if args.max_k < 3:
        parser.error("--max-k must be at least 3")
    if args.max_modulus < 2:
        parser.error("--max-modulus must be at least 2")
    path = Path(args.output)
    result = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    result.update({"arithmetic": "Python unbounded exact integers", "verification": verify()})
    started = time.process_time()
    result[f"round{args.round}"] = round1(args.max_k) if args.round == 1 else round2(args.max_modulus)
    print(f"round{args.round} cpu_seconds={time.process_time() - started:.6f}", file=sys.stderr)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
