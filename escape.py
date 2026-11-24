"""Search for exact returns after forced odd runs under odd multipliers 3, 5, 7.

Run python escape.py from the repository root. All trajectory arithmetic is
integer arithmetic. Floating point is used only for reporting base-2 logs.
"""

import json
import math
from pathlib import Path
import random
from time import process_time


CAP = 20_000
MULTIPLIERS = (3, 5, 7)
U_VALUES = (1, 3, 5)
J_VALUES = range(7, 13)
CONTROL_SEEDS = (173_205, 223_607)
SAMPLE_EVERY = 125
OUTPUT = Path(__file__).resolve().parent / "results" / "escape.json"


def family(k, u, j):
    """For k=3 this is u*2^m-1. General k takes the (u+1)/2-th odd lift."""
    m = 1 << j
    divisor = k - 2
    v = pow(1 << m, -1, divisor)
    if v % 2 == 0:
        v += divisor
    v += (u - 1) * divisor
    start = ((v << m) - 1) // divisor
    even_endpoint = (v * k**m - 1) // divisor
    halvings = (even_endpoint & -even_endpoint).bit_length() - 1
    return start, v, even_endpoint, halvings


def odd_control(bits, rng):
    return (1 << (bits - 1)) | rng.getrandbits(bits - 1) | 1


def replay_cycle(repeated, k, length):
    """Lift a repeated shortcut state to the standard map with exact arithmetic."""
    x = repeated
    steps = 0
    for _ in range(length):
        if x % 2:
            x = k * x + 1
            steps += 1
        x //= 2
        steps += 1
        assert x != 1
    assert x == repeated
    return steps


def trace(start, k, m, forced_endpoint, tail_halvings, targets, sample):
    suffix_start = m + tail_halvings
    x = start
    seen = {x: 0}
    peak = x
    standard_steps = 0
    series = []
    hits = []
    deltas = {}
    post_peak = x
    post_base = x
    cycle = None

    for step in range(CAP + 1):
        if step == m:
            assert x == forced_endpoint
        if sample and step >= m:
            plot_time = step - m
            if plot_time <= 128 or plot_time % SAMPLE_EVERY == 0:
                series.append([plot_time, round(math.log2(x) - math.log2(forced_endpoint), 6)])
        if step == suffix_start:
            post_base = post_peak = x
            assert x % 2
        if step >= suffix_start:
            elapsed = step - suffix_start
            if elapsed in (128, 1024):
                deltas[str(elapsed)] = round(math.log2(x) - math.log2(post_base), 8)

        if step and x in targets:
            hits.append({"target_u": targets[x], "shortcut_step": step})
        if x == 1:
            status = "reached_one"
            break
        if step and x in seen:
            period = step - seen[x]
            cycle = {
                "entry_shortcut_step": seen[x], "repeated_hex": hex(x),
                "shortcut_period": period, "standard_period": replay_cycle(x, k, period),
            }
            status = "repeated"
            break
        if step == CAP:
            status = "cap"
            break

        seen[x] = step
        if x % 2:
            intermediate = k * x + 1
            standard_steps += 2
            assert not m <= step < suffix_start
        else:
            intermediate = x
            standard_steps += 1
            assert step >= m
        peak = max(peak, intermediate)
        if step >= suffix_start:
            post_peak = max(post_peak, intermediate)
        x = intermediate // 2

    assert step >= suffix_start
    elapsed = step - suffix_start
    final_gain = round(math.log2(x) - math.log2(post_base), 8)
    if sample and series[-1][0] != step - m:
        series.append([step - m, round(math.log2(x) - math.log2(forced_endpoint), 6)])
    row = {
        "seed_hex": hex(start), "bits": start.bit_length(),
        "forced_steps": m, "forced_tail_halvings": tail_halvings,
        "post_start_bits": post_base.bit_length(), "post_steps": elapsed,
        "shortcut_steps": step, "standard_steps": standard_steps, "status": status,
        "peak_hex": hex(peak), "peak_log2": round(math.log2(peak), 8),
        "forced_endpoint_log2": round(math.log2(forced_endpoint), 8),
        "post_start_log2": round(math.log2(post_base), 8),
        "post_peak_hex": hex(post_peak),
        "post_peak_gain_bits": round(math.log2(post_peak) - math.log2(post_base), 8),
        "post_final_gain_bits": final_gain, "post_delta_bits": deltas,
        "terminal_hex": hex(x), "cycle": cycle, "return_hits": hits,
    }
    return row, series


def check():
    for k in MULTIPLIERS:
        for u in U_VALUES:
            for j in range(1, 5):
                x, _, end, tail = family(k, u, j)
                for _ in range(1 << j):
                    assert x % 2
                    x = k * x + 1
                    assert x % 2 == 0
                    x //= 2
                assert x == end and (x >> tail) % 2 == 1
                if k == 3 and u == 1:
                    assert tail == j + 2
    for seed, expected_steps, expected_peak in [(1, 0, 1), (3, 7, 16), (27, 111, 9232)]:
        row, _ = trace(seed, 3, 0, seed, 0, {}, True)
        assert row["standard_steps"] == expected_steps and row["status"] == "reached_one"
        assert int(row["peak_hex"], 16) == expected_peak
    row, _ = trace(13, 5, 0, 13, 0, {}, False)
    assert row["status"] == "repeated"
    assert row["cycle"]["shortcut_period"] == 7 and row["cycle"]["standard_period"] == 10


def run_round(number):
    clock_start = process_time()
    rng = random.Random(CONTROL_SEEDS[number - 1])
    rows = []
    figures = []
    for k in MULTIPLIERS:
        for u in U_VALUES:
            for j in J_VALUES:
                x, v, end, tail = family(k, u, j)
                m = 1 << j
                odd_endpoint = end >> tail
                control_bits = (x if number == 1 else odd_endpoint).bit_length()
                control = odd_control(control_bits, rng)
                targets = {}
                if number == 2 and k == 3:
                    targets = {family(k, other, j + 1)[0]: other for other in U_VALUES}
                for kind, seed in (("structured", x), ("control", control)):
                    structured = kind == "structured"
                    sample = number == 2 and structured and u == 1 and j == 8
                    row, series = trace(
                        seed, k, m if structured else 0, end if structured else seed,
                        tail if structured else 0, targets if structured else {}, sample,
                    )
                    identity = {"round": number, "k": k, "u": u, "v": v, "j": j, "kind": kind}
                    rows.append({**identity, **row})
                    if sample:
                        figures.append({"k": k, "u": u, "j": j, "points": series})
                if process_time() - clock_start > 900:
                    raise RuntimeError(f"Round {number} exceeded its 15-minute CPU budget")
    return rows, figures


def main():
    check()
    rows = []
    figures = []
    for number in (1, 2):
        batch, series = run_round(number)
        rows.extend(batch)
        figures.extend(series)
    output = {
        "metadata": {
            "shortcut_cap": CAP, "round_cpu_limit_seconds": 900, "u": U_VALUES,
            "j": list(J_VALUES), "k": MULTIPLIERS,
            "round_rng_seeds": CONTROL_SEEDS,
            "rounds": {
                "1": "Control bit length equals the structured seed's bit length.",
                "2": "Control bit length equals the first odd state after the forced climb. Test returns to any F(3,u',j+1). Structured paths repeat round 1.",
            },
            "control_sampling": "Python random.Random(seed), getrandbits(bits-1) with top and low bits set, one draw per k,u,j in listed order.",
            "stopping": "Stop on first 1, first repeated shortcut value, or exactly 20000 shortcut steps (in that priority order).",
            "map": "Even x maps to x/2, odd x maps to (k*x+1)/2. An odd shortcut step counts as two standard steps. Peaks include intermediate standard values.",
            "family": "m=2^j, d=k-2, v is the (u+1)/2-th positive odd solution of 2^m*v=1 mod d; F=(2^m*v-1)/d.",
            "suffix": "Start at the first odd value after the guaranteed m odd shortcut steps and following even steps, or at the odd control seed. Suffix step deltas use actual shortcut values, whereas post_peak_hex includes intervening standard values.",
            "figures": "Structured u=1,j=8, round 2. The forced even endpoint is time zero. Points are [elapsed shortcut steps, log2(value/forced endpoint)] at each step through 128, every 125 steps thereafter, and at the last state.",
            "return_test": "For each structured round-2 k=3 seed, test every positive shortcut time through the stop for any F(3,u',j+1), u' in {1,3,5}. All targets are odd. An empty hit list for other rows means the test was inapplicable.",
        }, "rows": rows, "figures": figures,
    }
    OUTPUT.parent.mkdir(exist_ok=True)
    OUTPUT.write_text(json.dumps(output, indent=1) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} rows and {len(figures)} traces to {OUTPUT}")


if __name__ == "__main__":
    main()
