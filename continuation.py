"""Resume the three censored 3n+1 family paths. Run python -B continuation.py."""

import json
from itertools import combinations
from math import log2
from pathlib import Path

from collatz import shortcut_step, standard_step


RESULTS = Path(__file__).resolve().parent / "results"
J = 12
U_VALUES = (1, 3, 5)
HISTORICAL_CAP = 20_000
ADDITIONAL_CAP = 100_000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def standard_checkpoints(seed):
    """Run the standard map one operation at a time and report every halving."""
    x = peak = seed
    standard_steps = 0
    yield standard_steps, x, peak
    while True:
        halving = x % 2 == 0
        x = standard_step(x)
        standard_steps += 1
        peak = max(peak, x)
        if halving:
            yield standard_steps, x, peak


def follow(history):
    old = history[1]
    u = old["u"]
    m = 1 << J
    seed = (u << m) - 1
    targets = {(other << (m * 2)) - 1: other for other in U_VALUES}
    replay = standard_checkpoints(seed)
    first_seen = {}
    hits = []
    points = []
    first_descent = None
    cycle = None
    status = "cap"

    for step in range(HISTORICAL_CAP + ADDITIONAL_CAP + 1):
        replay_standard, replay_x, replay_peak = next(replay)
        if step <= HISTORICAL_CAP:
            x, standard_steps, peak = replay_x, replay_standard, replay_peak
            if step == m:
                require(x == u * 3**m - 1, "forced endpoint differs from its formula")
        else:
            odd = x % 2
            x = shortcut_step(x)
            standard_steps += 1 + odd
            peak = max(peak, 2 * x if odd else x)
            require((x, standard_steps, peak) == (replay_x, replay_standard, replay_peak),
                    f"independent standard-map mismatch for u={u} at shortcut step {step}")

        if step and x in targets:
            hits.append({"target_u": targets[x], "shortcut_step": step, "standard_step": standard_steps})
        if step == HISTORICAL_CAP:
            for row in history:
                require((row["status"], row["seed_hex"], row["terminal_hex"], row["peak_hex"],
                         row["shortcut_steps"], row["standard_steps"])
                        == ("cap", hex(seed), hex(x), hex(peak), HISTORICAL_CAP, standard_steps),
                        f"historical round {row['round']} differs from independent replay")
            historical_hits = [{"target_u": hit["target_u"], "shortcut_step": hit["shortcut_step"]}
                               for hit in hits]
            require(historical_hits == old["return_hits"], "historical return hits differ")
            # Resume from the saved integer and counters, independently of the replay generator.
            x, standard_steps, peak = (int(old["terminal_hex"], 16), old["standard_steps"],
                                       int(old["peak_hex"], 16))

        if first_descent is None and x < seed:
            first_descent = step
        if x == 1:
            status = "reached_one"
        elif x in first_seen:
            entry_shortcut, entry_standard = first_seen[x]
            cycle = {"entry_shortcut_step": entry_shortcut, "shortcut_period": step - entry_shortcut,
                     "standard_period": standard_steps - entry_standard, "repeated_hex": hex(x)}
            status = "repeated"
        stop = status != "cap" or step == HISTORICAL_CAP + ADDITIONAL_CAP
        if step % 128 == 0 or step in (m, first_descent, HISTORICAL_CAP) or stop:
            points.append([step, round(log2(x), 6)])
        if stop:
            require(step > HISTORICAL_CAP, "historical stopping point differs from replay")
            break
        first_seen[x] = (step, standard_steps)

    row = {
        "u": u, "j": J, "status": status, "cycle": cycle,
        "old_shortcut_steps": HISTORICAL_CAP, "added_shortcut_steps": step - HISTORICAL_CAP,
        "total_shortcut_steps": step, "old_standard_steps": old["standard_steps"],
        "added_standard_steps": standard_steps - old["standard_steps"],
        "total_standard_steps": standard_steps, "first_below_seed_shortcut_step": first_descent,
        "checkpoint_bits": int(old["terminal_hex"], 16).bit_length(),
        "terminal_hex": hex(x), "peak_hex": hex(peak), "peak_bits": peak.bit_length(),
        "new_return_hits": [hit for hit in hits if hit["shortcut_step"] > HISTORICAL_CAP],
        "all_return_hits": hits, "standard_replay_verified": True,
    }
    first_seen.setdefault(x, (step, standard_steps))
    return row, {"u": u, "j": J, "points": points}, first_seen


def main():
    source = json.loads((RESULTS / "escape.json").read_text(encoding="utf-8"))
    rows, figures, visits = [], [], {}
    for u in U_VALUES:
        history = sorted((row for row in source["rows"] if row["k"] == 3 and row["j"] == J
                          and row["u"] == u and row["kind"] == "structured"), key=lambda r: r["round"])
        require([row["round"] for row in history] == [1, 2], "expected both historical rounds")
        row, figure, visits[u] = follow(history)
        rows.append(row)
        figures.append(figure)
    merges = []
    for left, right in combinations(U_VALUES, 2):
        for state, (right_shortcut, right_standard) in visits[right].items():
            if state in visits[left]:
                left_shortcut, left_standard = visits[left][state]
                merges.append({"left_u": left, "right_u": right, "shared_hex": hex(state),
                               "shared_bits": state.bit_length(), "left_shortcut_steps": left_shortcut,
                               "right_shortcut_steps": right_shortcut, "left_standard_steps": left_standard,
                               "right_standard_steps": right_standard})
                break
    output = {
        "metadata": {
            "source": "results/escape.json, structured k=3, j=12, rounds 1 and 2.",
            "prediction": "All three reach 1 in the extra budget and none visits a target.",
            "start": "For each u in {1,3,5}, resume u*2^4096-1 at its saved 20000-step endpoint.",
            "additional_shortcut_cap": ADDITIONAL_CAP, "total_shortcut_cap": HISTORICAL_CAP + ADDITIONAL_CAP,
            "stopping": "First 1, repeated shortcut state anywhere since the original seed, or total cap.",
            "return_test": "Check every positive shortcut time for u'*2^8192-1, u' in {1,3,5}. The odd targets cannot occur at a skipped even state.",
            "verification": "One continuous standard-step replay from each original seed checks both historical rows. Each continued shortcut value, cumulative standard count and peak agrees with the independent replay after its corresponding halving.",
            "figure_points": "[shortcut steps since original seed, log2(shortcut value)], every 128 steps plus first descent, the 4096-step forced endpoint, the old cutoff and terminal.",
            "merges": "First common shortcut state encountered along the right path, with exact indices and value. Determinism gives a shared suffix.",
        },
        "rows": rows, "figures": figures, "merges": merges,
    }
    target = RESULTS / "continuation.json"
    target.write_text(json.dumps(output, indent=1) + "\n", encoding="utf-8")
    for row in rows:
        print(row["u"], row["status"], row["added_shortcut_steps"], row["total_standard_steps"],
              len(row["all_return_hits"]))


if __name__ == "__main__":
    main()
