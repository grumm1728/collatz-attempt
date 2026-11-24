# Collatz counterexample search

## Purpose

Try to disprove the Collatz conjecture by finding a positive integer whose trajectory never reaches 1, or by constructing a rigorous argument that such an integer exists. Use the standard map: divide even integers by 2 and replace odd integers with `3n + 1`.

## Research and implementation

- Treat a cycle that excludes 1 or a proof that a specific trajectory never reaches 1 as a possible disproof. A long trajectory, a timeout, or a large intermediate value is only a lead.
- Use exact integer arithmetic for trajectory calculations. If an implementation uses bounded integers, detect overflow before it can affect a result.
- When using an accelerated map or another formulation, document its relationship to the standard map and distinguish accelerated steps from standard steps.
- Keep proven statements, assumptions, and experimental observations distinct. Record the search bounds, stopping rules, and commands needed to reproduce computational results.

## Verification

- Check trajectory code on known cases, including the `1 → 4 → 2 → 1` cycle, and compare optimized calculations against a straightforward implementation.
- Independently verify any proposed counterexample against the standard map with exact arithmetic. For a theoretical argument, state the assumptions and justify each step needed to rule out reaching 1.
