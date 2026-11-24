"""Certify finite exclusion of immediate Mersenne recharge under Collatz.

Starting at 2**m-1, the first m-1 odd-to-odd steps each divide once after
3*x+1. The next odd value is (3**m-1)/2**t, where t=v2(3**m-1). An
immediate recharge means that value is -1 modulo 2**m, so the next
one-division climb could have the same length or greater.

For odd m>=3, t=1 and 3**m is 3 modulo 8, whereas recharge requires 7.
For even m=2**a*b with b odd, t=a+2. Indeed 3**2=9 has t=3,
each squaring raises t by one, and an odd power preserves t.
The same argument gives v2((3**(2**a))**d-1)=a+2+v2(d) for d>0.
Thus at depth n there's at most one b modulo 2**n that can satisfy
3**m = 1-2**t modulo 2**(t+n). A verified least residue above the
allowed odd factors excludes all m in its row with m>=n.
"""

import argparse
import json
from pathlib import Path


DEFAULT_BITS = 20


def odd_after_climb(m):
    """Evaluate the odd state after the initial Mersenne climb."""
    odd_numerator = 3**m - 1
    return odd_numerator // (odd_numerator & -odd_numerator)


def immediately_recharges(m):
    return (odd_after_climb(m) + 1) % (1 << m) == 0


def lift_row(bits, valuation):
    """Find a residue above the entire range's attainable odd factors."""
    maximum_b = ((1 << bits) - 1) >> valuation
    exponent_offset = valuation + 2
    target = 1 - (1 << exponent_offset)
    maximum_depth = bits + 32
    working_modulus = 1 << (exponent_offset + maximum_depth)
    factor = pow(3, 1 << valuation, working_modulus)
    power = 1
    residue = 0
    for depth in range(1, maximum_depth + 1):
        modulus = 1 << (exponent_offset + depth)
        if (power - target) % modulus:
            residue += 1 << (depth - 1)
            power = power * factor % working_modulus
        if residue > maximum_b:
            return {"valuation": valuation, "depth": depth, "least_b": residue}
        factor = factor * factor % working_modulus
    raise ArithmeticError(f"no excluding residue for valuation {valuation} by depth {maximum_depth}")


def verify_rows(bits, rows):
    """Verify all row witnesses directly, then check the small exponents."""
    if bits < 2:
        raise ValueError("bits must be at least 2")
    if [row["valuation"] for row in rows] != list(range(1, bits)):
        raise ValueError("one row is required for each possible positive valuation")

    limit = 1 << bits
    for row in rows:
        valuation, depth, residue = row["valuation"], row["depth"], row["least_b"]
        maximum_b = (limit - 1) >> valuation
        if depth < 1 or not maximum_b < residue < (1 << depth) or residue % 2 != 1:
            raise ValueError(f"invalid excluding residue at valuation {valuation}")
        offset = valuation + 2
        modulus = 1 << (offset + depth)
        if pow(3, residue << valuation, modulus) != (1 - (1 << offset)) % modulus:
            raise ValueError(f"incorrect residue at valuation {valuation}")

    max_depth = max(row["depth"] for row in rows)
    scan_end = min(limit, max(3, max_depth))
    scan = range(3, scan_end)
    if any(immediately_recharges(m) for m in scan):
        raise ArithmeticError("a recharge occurred in the exact small-exponent scan")
    return max_depth, scan_end


def certificate(bits, rows):
    max_depth, scan_end = verify_rows(bits, rows)
    return {
        "prediction": (
            "Every valuation row has an excluding residue by depth bits+32, "
            "and the exact small-exponent scan finds no recharge."
        ),
        "inclusive_lower_m": 3,
        "exclusive_upper_m": 1 << bits,
        "bits": bits,
        "maximum_lift_depth": max_depth,
        "direct_scan": {
            "inclusive_start_m": 3,
            "exclusive_end_m": scan_end,
            "tested": scan_end - 3,
            "recharges": 0,
        },
        "residue_rows": rows,
        "verified_no_immediate_recharge": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=DEFAULT_BITS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path, help="verify a saved certificate directly")
    args = parser.parse_args()
    if args.verify:
        saved = json.loads(args.verify.read_text(encoding="utf-8"))
        report = certificate(saved["bits"], saved["residue_rows"])
        if saved != report:
            parser.error("the saved certificate's summary disagrees with its verified rows")
    else:
        if args.bits < 2:
            parser.error("--bits must be at least 2")
        report = certificate(args.bits, [lift_row(args.bits, a) for a in range(1, args.bits)])
    rendered = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
