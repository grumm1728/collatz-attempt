"""Search compressed odd-to-odd Collatz words with exact modular arithmetic."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import cache, lru_cache
from itertools import chain
import json
from pathlib import Path
from random import Random


MOTIFS = ((1, 2), (2,))
ORIGINAL_COUNTS = (59_813_072_756, 24_488_718_470)
PARTS = 16
ROUND_ONE_SEED = 2026112401
ROUND_TWO_SEED = 2026112402
ROUND_ONE_SAMPLES = 256
ROUND_TWO_SAMPLES_PER_LENGTH = 8
ROUND_TWO_LENGTHS = 128
FRONTIER = 2**71
MIN_ODD_STEPS = 137_500_000_000
MAX_ODD_STEPS = 10**13
CALIBRATION_PRIME = 31
HOLDOUT_PRIME = 7
DISCOVERY_PRIMES = (
    5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
)
LOG_SCALE = 10**70
LOG_TERMS = 128


def standard_odd_step(value):
    """Take one 3x+1 step, then count the divisions until the next odd value."""
    if value <= 0 or value % 2 == 0:
        raise ValueError("expected a positive odd integer")
    value = 3 * value + 1
    divisions = 0
    while value % 2 == 0:
        value //= 2
        divisions += 1
    return value, divisions


def lengths(counts):
    """Return odd steps and divisions for a pair of motif counts."""
    u, v = counts
    return 2 * u + v, 3 * u + 2 * v


def motif_counts(odd_steps, divisions):
    return 2 * odd_steps - divisions, 2 * divisions - 3 * odd_steps


def concatenate(first, second, modulus):
    """Compose (multiplier, denominator, offset) without any modular division."""
    a, b, c = first
    d, e, f = second
    return a * d % modulus, b * e % modulus, (c * d + b * f) % modulus


@lru_cache(maxsize=8192)
def modular_power(motif, copies, modulus):
    identity = (1, 1, 0)
    base = identity
    for divisions in MOTIFS[motif]:
        base = concatenate(base, (3, 2**divisions, 1), modulus)
    result = identity
    while copies:
        if copies % 2:
            result = concatenate(result, base, modulus)
        base = concatenate(base, base, modulus)
        copies //= 2
    return result


def modular_affine(schedule, modulus):
    result = (1, 1, 0)
    for motif, copies in schedule:
        result = concatenate(result, modular_power(motif, copies, modulus), modulus)
    return result


@cache
def residue_orbits(motif, prime):
    """Build the orbit of each residue using the literal valuation transitions."""
    permutation = []
    for residue in range(prime):
        value = residue
        for divisions in MOTIFS[motif]:
            value = (3 * value + 1) % prime
            for _ in range(divisions):
                if value % 2:
                    value += prime
                value //= 2
        permutation.append(value)

    if len(set(permutation)) != prime:
        raise ValueError("the motif map must permute the residues")
    orbits = []
    for residue in range(prime):
        orbit = [residue]
        value = permutation[residue]
        while value != residue:
            orbit.append(value)
            value = permutation[value]
        orbits.append(tuple(orbit))
    return tuple(orbits)


def residue_map(schedule, prime):
    """Independently compose the complete residue permutation."""
    mapping = list(range(prime))
    for motif, copies in schedule:
        orbits = residue_orbits(motif, prime)
        for start, value in enumerate(mapping):
            orbit = orbits[value]
            mapping[start] = orbit[copies % len(orbit)]
    return mapping


def checked_translation(schedule, prime):
    """Check a closing modulus by direct affine and full residue calculations."""
    a, b, c = modular_affine(schedule, prime)
    if a != b:
        raise ValueError(f"{prime} does not divide the closing denominator")
    translation = c * pow(b, -1, prime) % prime
    expected = [(residue + translation) % prime for residue in range(prime)]
    if residue_map(schedule, prime) != expected:
        raise ArithmeticError(f"modular composition disagrees modulo {prime}")
    return translation


@cache
def log_bounds(reciprocal):
    """Bound log((1+z)/(1-z)) in units of 1/LOG_SCALE, for z=1/reciprocal."""
    lower = 0
    power = reciprocal
    square = reciprocal**2
    for index in range(LOG_TERMS):
        lower += 2 * LOG_SCALE // ((2 * index + 1) * power)
        power *= square
    tail_numerator = 2 * LOG_SCALE * square
    tail_denominator = (2 * LOG_TERMS + 1) * power * (square - 1)
    tail_ceiling = (tail_numerator + tail_denominator - 1) // tail_denominator
    return lower, lower + LOG_TERMS + tail_ceiling


def bounded_slopes(lower, upper, maximum_denominator):
    """Enumerate reduced rationals strictly in a window via Stern–Brocot."""
    lower_numerator, lower_denominator = lower
    upper_numerator, upper_denominator = upper
    stack = [(0, 1, 1, 0)]
    while stack:
        a, b, c, d = stack.pop()
        while (
            a * upper_denominator < b * upper_numerator
            and c * lower_denominator > d * lower_numerator
            and b + d <= maximum_denominator
        ):
            numerator, denominator = a + c, b + d
            if numerator * lower_denominator <= denominator * lower_numerator:
                # Jump over the entire chain of mediants below the window.
                jump = (lower_numerator * b - lower_denominator * a) // (
                    lower_denominator * c - lower_numerator * d
                )
                a, b = a + jump * c, b + jump * d
            elif numerator * upper_denominator >= denominator * upper_numerator:
                jump = (upper_denominator * c - upper_numerator * d) // (
                    upper_numerator * b - upper_denominator * a
                )
                c, d = c + jump * a, d + jump * b
            else:
                yield numerator, denominator
                stack.append((a, b, numerator, denominator))
                a, b = numerator, denominator


def parameter_pairs(minimum, maximum, frontier):
    """Exhaust the necessary cycle length window, including nonreduced slopes."""
    log2_lower, log2_upper = log_bounds(3)
    log3_lower, log3_upper = log_bounds(2)
    three_frontier = 3 * frontier
    lower = (log3_lower, log2_upper)
    upper = (three_frontier * log3_upper + LOG_SCALE, three_frontier * log2_lower)
    pairs = []
    for numerator, denominator in bounded_slopes(lower, upper, maximum):
        low_gap = numerator * log2_lower - denominator * log3_upper
        high_gap = numerator * log2_upper - denominator * log3_lower
        if not (0 < low_gap and three_frontier * high_gap < denominator * LOG_SCALE):
            raise ArithmeticError("a slope meets only the outer interval")
        first = (minimum + denominator - 1) // denominator
        last = maximum // denominator
        pairs.extend((denominator * k, numerator * k) for k in range(first, last + 1))
    return sorted(pairs)


def inverse_packet_weight_caps():
    """Bound a motif power's inverse offset per start or end boundary weight."""
    caps = []
    for motif in MOTIFS:
        multiplier = Fraction(1)
        offset = Fraction(0)
        for divisions in motif:
            offset += multiplier / 3
            multiplier *= Fraction(2**divisions, 3)
        caps.append(offset / abs(1 - multiplier))
    return tuple(caps)


def packet_minimum_certificate(pairs, frontier):
    """Bound one rotation of every word with at most PARTS packets per motif.

    Rotate to a packet boundary maximizing the centered division prefix
    s - j*S/K. Every rotated boundary then has centered prefix at most zero,
    including after a cyclic wrap. Its weight 2**s/3**j is at most exp(gap),
    where gap = S*log(2) - K*log(3) is positive.

    The inverse affine offset of a packet is a geometric sum. For a motif
    with inverse multiplier below one, bound the sum using its start
    boundary. For a multiplier above one, use the end boundary. Dividing the
    whole offset by exp(gap)-1 gives the formal closing value at the rotation.
    Consequently the minimum is below total_cap/(1-exp(-gap)), which is
    below total_cap*(1+1/gap) by exp(gap) > 1+gap.
    """
    log2_lower, _ = log_bounds(3)
    _, log3_upper = log_bounds(2)
    gap_lower, minimizing_pair = min(
        (divisions * log2_lower - odd_steps * log3_upper, (odd_steps, divisions))
        for odd_steps, divisions in pairs
    )
    if gap_lower <= 0:
        raise ArithmeticError("the minimum bound needs a positive closing gap")

    caps = inverse_packet_weight_caps()
    total_cap = PARTS * sum(caps)
    strict_bound = total_cap * (1 + Fraction(LOG_SCALE, gap_lower))
    ceiling = (strict_bound.numerator + strict_bound.denominator - 1) // strict_bound.denominator
    return {
        "parameter_pairs_covered": len(pairs),
        "maximum_packets_per_motif": PARTS,
        "motif_weight_caps": [str(cap) for cap in caps],
        "total_weight_cap": str(total_cap),
        "closing_gap_lower_bound": {
            "numerator": gap_lower,
            "denominator": LOG_SCALE,
            "length_pair": minimizing_pair,
            "log_series_terms": LOG_TERMS,
        },
        "strict_upper_bound_on_minimum": ceiling,
        "lower_bound_on_nontrivial_cycle_minimum": frontier,
        "all_covered_packet_words_excluded": strict_bound <= frontier,
    }


def positive_parts(total, parts, rng):
    """Sample a positive ordered composition uniformly through its cut positions."""
    cuts = [0, *sorted(rng.sample(range(1, total), parts - 1)), total]
    return [end - start for start, end in zip(cuts, cuts[1:])]


def matched_schedules(counts, rng):
    """Use the same packet multiset for alternating order and a sampled shuffle."""
    u_parts = positive_parts(counts[0], PARTS, rng)
    v_parts = positive_parts(counts[1], PARTS, rng)
    alternating = list(
        chain.from_iterable(((0, u), (1, v)) for u, v in zip(u_parts, v_parts))
    )
    shuffled = alternating.copy()
    rng.shuffle(shuffled)
    return alternating, shuffled


def translation_summary(translations, prime):
    samples = sum(translations.values())
    survivors = translations[0]
    return {
        "tested": samples,
        "rejected": samples - survivors,
        "surviving_this_modulus": survivors,
        "fixed_residues_for_each_survivor": prime,
        "counts_by_translation": [translations[residue] for residue in range(prime)],
    }


def round_one():
    prime = CALIBRATION_PRIME
    counts = ORIGINAL_COUNTS
    original = list(enumerate(counts))
    if standard_odd_step(1) != (1, 2):
        raise ArithmeticError("the standard-map control did not follow 1, 4, 2, 1")
    if residue_map([(1, 1)], prime)[1] != 1:
        raise ArithmeticError("the modular verifier did not reproduce the trivial cycle")
    original_translation = checked_translation(original, prime)
    if original_translation != 7:
        raise ArithmeticError("original calibration did not reproduce translation +7")

    periods = [max(map(len, residue_orbits(motif, prime))) for motif in range(2)]
    # Each motif permutation has a single common nontrivial cycle length.
    for motif, period in enumerate(periods):
        if any(period % len(orbit) for orbit in residue_orbits(motif, prime)):
            raise ArithmeticError("split representatives do not cover the motif order")
    if any(count <= period for count, period in zip(counts, periods)):
        raise ArithmeticError("a split representative exceeds its available copies")

    split_translations = Counter()
    for u in range(1, periods[0] + 1):
        for v in range(1, periods[1] + 1):
            split = [(0, u), (1, v), (0, counts[0] - u), (1, counts[1] - v)]
            split_translations[checked_translation(split, prime)] += 1

    rng = Random(ROUND_ONE_SEED)
    alternating, shuffled = Counter(), Counter()
    matched_outcomes = Counter()
    for _ in range(ROUND_ONE_SAMPLES):
        left, right = matched_schedules(counts, rng)
        left_translation = checked_translation(left, prime)
        right_translation = checked_translation(right, prime)
        alternating[left_translation] += 1
        shuffled[right_translation] += 1
        matched_outcomes[(bool(left_translation), bool(right_translation))] += 1

    odd_steps, divisions = lengths(counts)
    if (odd_steps, divisions) not in parameter_pairs(odd_steps, odd_steps, FRONTIER):
        raise ArithmeticError("the calibration falls outside the necessary length window")
    a, b, c = modular_affine(original, prime)
    split_summary = translation_summary(split_translations, prime)
    alternating_summary = translation_summary(alternating, prime)
    shuffled_summary = translation_summary(shuffled, prime)
    alternating_rejections = alternating_summary["rejected"]
    shuffled_rejections = shuffled_summary["rejected"]
    return {
        "prediction": (
            "Original translation +7 and no fixed residue. The exhaustive split classes "
            "include rejections and survivors. At least 90% rejection in each sampled "
            "ordering, with a difference of at most ten percentage points."
        ),
        "prediction_met": (
            split_summary["rejected"] > 0
            and split_summary["surviving_this_modulus"] > 0
            and 10 * min(alternating_rejections, shuffled_rejections) >= 9 * ROUND_ONE_SAMPLES
            and 10 * abs(alternating_rejections - shuffled_rejections) <= ROUND_ONE_SAMPLES
        ),
        "modulus": prime,
        "motifs": MOTIFS,
        "motif_counts": counts,
        "odd_steps": odd_steps,
        "divisions": divisions,
        "standard_steps_if_realized": odd_steps + divisions,
        "original": {
            "affine_tuple_modulo_prime": (a, b, c),
            "translation": original_translation,
            "fixed_residues": 0,
        },
        "packet_minimum_certificate": packet_minimum_certificate(
            [(odd_steps, divisions)], FRONTIER
        ),
        "two_split_exhaustive_modulo_prime": {
            "motif_orders": periods,
            **split_summary,
        },
        "sampled_matched_templates": {
            "seed": ROUND_ONE_SEED,
            "packets_per_motif": PARTS,
            "alternating": alternating_summary,
            "shuffled": shuffled_summary,
            "both_rejected": matched_outcomes[True, True],
            "only_alternating_rejected": matched_outcomes[True, False],
            "only_shuffled_rejected": matched_outcomes[False, True],
            "both_survived": matched_outcomes[False, False],
        },
        "independent_residue_permutations_checked": (
            1 + sum(split_translations.values()) + 2 * ROUND_ONE_SAMPLES
        ),
    }


def denominator_divisible(odd_steps, divisions, prime):
    return pow(2, divisions, prime) == pow(3, odd_steps, prime)


def round_two():
    pairs = parameter_pairs(MIN_ODD_STEPS, MAX_ODD_STEPS, FRONTIER)
    holdout_eligible = [
        pair for pair in pairs if denominator_divisible(*pair, HOLDOUT_PRIME)
    ]
    selected = []
    eligible_count = 0
    for pair in holdout_eligible:
        moduli = [prime for prime in DISCOVERY_PRIMES if denominator_divisible(*pair, prime)]
        if moduli:
            eligible_count += 1
            if len(selected) < ROUND_TWO_LENGTHS:
                selected.append((pair, moduli))
    if len(selected) != ROUND_TWO_LENGTHS:
        raise ArithmeticError("too few length pairs for the declared sample")

    minimum_certificate = packet_minimum_certificate(pairs, FRONTIER)
    packet_bound_excludes_word = minimum_certificate["all_covered_packet_words_excluded"]
    rng = Random(ROUND_TWO_SEED)
    outcome_keys = (
        "tested",
        "rejected_by_discovery",
        "discovery_survivors",
        "rejected_by_holdout",
        "survived_all_tested_moduli",
        "independent_permutation_checks",
    )
    outcomes = {name: dict.fromkeys(outcome_keys, 0) for name in ("alternating", "shuffled")}
    witnesses = {}
    for (odd_steps, divisions), moduli in selected:
        counts = motif_counts(odd_steps, divisions)
        for _ in range(ROUND_TWO_SAMPLES_PER_LENGTH):
            for name, schedule in zip(outcomes, matched_schedules(counts, rng)):
                stats = outcomes[name]
                stats["tested"] += 1
                for prime in moduli:
                    stats["independent_permutation_checks"] += 1
                    if checked_translation(schedule, prime):
                        stats["rejected_by_discovery"] += 1
                        break
                else:
                    stats["discovery_survivors"] += 1
                    stats["independent_permutation_checks"] += 1
                    translation = checked_translation(schedule, HOLDOUT_PRIME)
                    if translation:
                        stats["rejected_by_holdout"] += 1
                        if name not in witnesses:
                            witnesses[name] = {
                                "odd_steps": odd_steps,
                                "divisions": divisions,
                                "discovery_moduli": moduli,
                                "holdout_translation": translation,
                                "schedule": schedule,
                            }
                    else:
                        stats["survived_all_tested_moduli"] += 1

    return {
        "prediction": (
            "Each ordering has at least one discovery survivor, and the held-out "
            "prime rejects at least half of the discovery survivors in each ordering."
        ),
        "prediction_met": all(
            stats["discovery_survivors"] > 0
            and 2 * stats["rejected_by_holdout"] >= stats["discovery_survivors"]
            for stats in outcomes.values()
        ),
        "lower_bound_on_cycle_minimum": FRONTIER,
        "inclusive_odd_step_bounds": [MIN_ODD_STEPS, MAX_ODD_STEPS],
        "exhaustive_parameter_window_pairs": len(pairs),
        "holdout_eligible_pairs": len(holdout_eligible),
        "holdout_and_discovery_eligible_pairs": eligible_count,
        "selected_pairs": len(selected),
        "pair_selection": "smallest eligible odd step counts, then division counts",
        "first_selected_pair": selected[0][0],
        "last_selected_pair": selected[-1][0],
        "seed": ROUND_TWO_SEED,
        "packets_per_motif": PARTS,
        "matched_samples_per_length": ROUND_TWO_SAMPLES_PER_LENGTH,
        "motifs": MOTIFS,
        "discovery_primes": DISCOVERY_PRIMES,
        "holdout_prime": HOLDOUT_PRIME,
        "outcomes": outcomes,
        "first_holdout_rejection": witnesses,
        "packet_minimum_certificate": minimum_certificate,
        "modular_survivors_excluded_by_packet_minimum": {
            name: stats["survived_all_tested_moduli"] if packet_bound_excludes_word else 0
            for name, stats in outcomes.items()
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round", choices=("1", "2", "all"), default="all")
    parser.add_argument("--output", type=Path, help="update a JSON report at this path")
    args = parser.parse_args()

    runs = {"round_1": round_one, "round_2": round_two}
    names = runs if args.round == "all" else [f"round_{args.round}"]
    report = {name: runs[name]() for name in names}
    print(json.dumps(report, indent=2))
    if args.output:
        previous = json.loads(args.output.read_text(encoding="utf-8")) if args.output.exists() else {}
        previous.update(report)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(previous, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
