"""Exact positive-integer maps. A shortcut odd step includes one halving for odd k."""


def standard_step(n: int, k: int = 3) -> int:
    """Apply one standard kx+1 step; both arguments are positive integers."""
    if n < 1 or k < 1:
        raise ValueError("n and k must be positive")
    return k * n + 1 if n & 1 else n // 2


def shortcut_step(n: int, k: int = 3) -> int:
    """Apply one shortcut step, also supporting the odd-to-odd branch of even k."""
    if n < 1 or k < 1:
        raise ValueError("n and k must be positive")
    if n % 2 == 0:
        return n // 2
    expanded = k * n + 1
    return expanded // 2 if k % 2 else expanded


def odd_step(n: int, k: int = 3) -> tuple[int, int]:
    """Return the next odd value and number of standard halvings after k*n+1."""
    if n < 1 or n % 2 == 0 or k < 1:
        raise ValueError("expected a positive odd n and positive k")
    expanded = k * n + 1
    halvings = (expanded & -expanded).bit_length() - 1
    return expanded >> halvings, halvings
