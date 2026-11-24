"""Generate the three pilot figures and their exact cycle-atlas inputs."""

import json
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from collatz import odd_step

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
MULTIPLIERS = (1, 3, 5, 7)
INK, TEAL, ORANGE = "#18363f", "#087e82", "#be4a27"
PAPER = "#faf8f1"


def save(fig, name):
    for fmt in ("svg", "png"):
        metadata = {"Date": None} if fmt == "svg" else {"Software": "Collatz pilot"}
        path = OUT / f"{name}.{fmt}"
        fig.savefig(path, dpi=160, metadata=metadata,
                    facecolor=PAPER, bbox_inches="tight", pad_inches=0.25)
        if fmt == "svg":
            lines = path.read_text(encoding="utf-8").splitlines()
            path.write_text("\n".join(s.rstrip() for s in lines) + "\n", encoding="utf-8", newline="\n")
    plt.close(fig)


def observed_recipes():
    recipes = set()
    for k in MULTIPLIERS:
        for seed in range(1, 256, 2):
            seen, trail = {}, []
            n = seed
            for _ in range(129):
                if n in seen:
                    word = tuple(trail[seen[n]:])
                    recipes.add(min(word[i:] + word[:i] for i in range(len(word))))
                    break
                if len(trail) == 128:
                    break
                seen[n] = len(trail)
                n, exponent = odd_step(n, k)
                trail.append(exponent)
    return sorted(recipes, key=lambda word: (len(word), word))


def closing_seed(word, k):
    offset, power_of_two = 0, 1
    for exponent in word:
        offset = k * offset + power_of_two
        power_of_two <<= exponent
    seed = Fraction(offset, power_of_two - k**len(word))
    if seed > 0 and seed.denominator == 1:
        n = int(seed)
        for a in word:
            n, observed = odd_step(n, k)
            if observed != a:
                raise ArithmeticError("integral solution failed halving-count replay")
        if n != seed:
            raise ArithmeticError("cycle failed to close")
    return seed


def cycle_atlas():
    words = observed_recipes()
    assert (1, 1, 5) in words
    assert closing_seed((1, 1, 5), 3) == Fraction(19, 101)
    assert closing_seed((1, 1, 5), 5) == 13
    fig, ax = plt.subplots(figsize=(9, max(3.8, len(words) * 0.55 + 1.9)))
    fig.subplots_adjust(left=0.31, bottom=0.21, top=0.78)
    rows = []
    for y, word in enumerate(words):
        row = {}
        for x, k in enumerate(MULTIPLIERS):
            seed = closing_seed(word, k)
            row[str(k)] = str(seed)
            if seed <= 0:
                label, color = "−", "#826c68"
            elif seed.denominator == 1:
                label, color = str(seed), TEAL
            else:
                label, color = "○", "#486981"
            if word == (1, 1, 5) and k == 3:
                label = "19/101"
            ax.text(x, y, label, ha="center", va="center", color=color, weight="bold")
        rows.append({"halvings": list(word), "closing_seeds": row})
    ax.set(xlim=(-0.5, 3.5), ylim=(len(words)-0.5, -0.5),
           xticks=range(4), xticklabels=[f"k = {k}" for k in MULTIPLIERS],
           yticks=range(len(words)), yticklabels=[", ".join(map(str, w)) for w in words])
    for y in range(len(words)+1):
        ax.axhline(y-0.5, color="#d8d8ce", lw=0.8)
    ax.tick_params(length=0, pad=9)
    fig.text(0.04, 0.93, "Where cycles survive", fontsize=22, weight="bold")
    fig.text(0.04, 0.85, "A row fixes the halving counts. Each column changes the multiplier.")
    fig.text(0.04, 0.035, "Teal integer = verified cycle start     ○ or fraction = positive rational     − = negative root")
    save(fig, "cycle-atlas")
    (ROOT / "results" / "atlas.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


def tapestry():
    odd = range(1, 1024, 2)
    bands = [[min(odd_step(n, k)[1], 10) for n in odd] for k in MULTIPLIERS]
    expected = [256, 128, 64, 32, 16, 8, 4, 2, 1, 1]
    assert all([row.count(a) for a in range(1, 11)] == expected for row in bands)
    colors = ["#182a4c", "#20506a", "#176e78", "#048a83", "#42a18a",
              "#89bc8e", "#c0ce8a", "#edd482", "#efae61", "#d86a47"]
    cmap = ListedColormap(colors)
    fig, axes = plt.subplots(2, 1, figsize=(10.6, 5.7), height_ratios=[1.7, 1])
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.17, top=0.81, hspace=0.7)
    for ax, count in zip(axes, [512, 32]):
        ax.imshow([row[:count] for row in bands], interpolation="nearest", aspect="auto",
                  vmin=0.5, vmax=10.5, cmap=cmap)
        ax.set(yticks=range(4), yticklabels=[f"k = {k}" for k in MULTIPLIERS])
        ticks = [0, count//4, count//2, 3*count//4, count-1]
        ax.set(xticks=ticks, xticklabels=[str(2*i+1) for i in ticks])
        ax.tick_params(length=0)
        ax.set_title("All 512 odd residues" if count == 512 else "Detail of the first 32 odd seeds",
                     loc="left", fontsize=10)
    fig.text(0.04, 0.93, "The halving tapestry", fontsize=22, weight="bold")
    fig.text(0.04, 0.86, "Each row has the same color counts. The multiplier changes where the colors fall.")
    colors_ax = fig.add_axes([0.21, 0.04, 0.6, 0.035])
    colors_ax.imshow([list(range(1, 11))], cmap=cmap, aspect="auto", interpolation="nearest")
    colors_ax.set(yticks=[], xticks=range(10), xticklabels=[*map(str, range(1, 10)), "10+"])
    colors_ax.tick_params(length=0, labelsize=9)
    fig.text(0.04, 0.05, "Halvings", fontsize=10)
    save(fig, "halving-tapestry")


def after_climb():
    data = json.loads((ROOT / "results" / "escape.json").read_text(encoding="utf-8"))
    fig, ax = plt.subplots(figsize=(9.5, 5))
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.16, top=0.78)
    for row, color in zip(sorted(data["figures"], key=lambda r: r["k"]), [TEAL, ORANGE, "#65508b"]):
        points = [(x, y) for x, y in row["points"] if x <= 128]
        t, bits = zip(*points)
        outcome = next(r for r in data["rows"] if r["round"] == 2 and
                       r["kind"] == "structured" and
                       all(r[key] == row[key] for key in ("k", "u", "j")))
        status = {"reached_one": "reached 1", "repeated": "cycle",
                  "cap": "step cap"}[outcome["status"]]
        label = f"k = {row['k']}   {status} at {outcome['shortcut_steps']:,} total steps"
        ax.plot(t, bits, color=color, lw=1.9, label=label)
    ax.axhline(0, color=INK, ls=":", lw=0.9)
    ax.grid(axis="y", alpha=0.19)
    ax.set(xlim=(0, 128), xlabel="Shortcut steps after the forced prefix",
           ylabel="Change in log₂(value) from the endpoint")
    ax.legend(frameon=False)
    fig.text(0.04, 0.93, "After the climb", fontsize=22, weight="bold")
    fig.text(0.04, 0.85, "Three matched constructions, each aligned after 256 forced odd shortcut steps.")
    save(fig, "after-the-climb")


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": INK,
                         "axes.labelcolor": INK, "axes.facecolor": PAPER,
                         "xtick.color": INK, "ytick.color": INK,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.spines.left": False, "axes.spines.bottom": False,
                         "svg.hashsalt": "collatz-pilot", "svg.fonttype": "none"})
    OUT.mkdir(exist_ok=True)
    (ROOT / "results").mkdir(exist_ok=True)
    cycle_atlas()
    tapestry()
    after_climb()


if __name__ == "__main__":
    main()
