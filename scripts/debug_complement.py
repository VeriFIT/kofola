#!/usr/bin/env python3.10
"""
Debug utility for kofola complement verification.

For a given input automaton, runs kofola with raw output (preserving macrostate
labels as state names), checks equivalence with Spot's complement, and if they
differ finds a violating accepted run in the kofola automaton and visualizes it
with full macrostate labels at each step.

Usage:
    python3 debug_complement.py --input <hoa_file> [options]

Example (inductive det + OR-FIN opt):
    python3 debug_complement.py \\
        --input tests/test_data/tela_det/automaton327.hoa \\
        --params "tela=yes;tela_det_alg=inductive;sd_ind_or_opt=yes" \\
        --kofola build/src/kofola

Requirements:
    - Spot Python bindings (install Spot with Python bindings or set PYTHONPATH)
    - kofola binary built with raw=yes + show-macrostate-labels=yes support
"""

import argparse
import subprocess
import sys

try:
    import spot  # noqa: F401
except ImportError:
    sys.exit(
        "ERROR: Could not import 'spot' Python bindings.\n"
        "Install Spot with Python bindings or add its site-packages to PYTHONPATH."
    )

import buddy  # BDD operations (shipped alongside spot)


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def run_kofola(kofola_path: str, input_file: str, params: str) -> str:
    """Run kofola and return stdout (HOA string)."""
    cmd = [kofola_path, "--params", params, input_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"kofola exited with code {result.returncode}.\n"
            f"stderr:\n{result.stderr}"
        )
    return result.stdout


def get_state_name(aut, state_num: int) -> str:
    """Return the macrostate label of a state, or its number as fallback."""
    names = aut.get_state_names()
    if names and state_num < len(names):
        name = names[state_num]
        if name:
            return name
    return str(state_num)


def bdd_to_str(bdd_val, d) -> str:
    """Convert a BDD to a human-readable formula string."""
    return str(spot.bdd_to_formula(bdd_val, d))


def condition_to_str(cond, d) -> str:
    """Convert a trace condition to text, preserving pre-rendered labels."""
    if isinstance(cond, str):
        return cond
    return bdd_to_str(cond, d)


def _parse_product_state(product_aut, state) -> tuple:
    """Return the integer component-state numbers encoded in a Spot product state."""
    formatted = product_aut.format_state(state)
    parts = [part.strip() for part in formatted.split(",")]
    if len(parts) != 2:
        raise RuntimeError(f"Unexpected product state format: {formatted!r}")

    try:
        return (int(parts[0]), int(parts[1]))
    except ValueError as exc:
        raise RuntimeError(f"Could not parse product state numbers from {formatted!r}") from exc


def project_product_run(product_aut, run) -> tuple:
    """Project a product accepting run to the first component automaton trace."""
    d = product_aut.get_dict()
    prefix_trace = [
        (_parse_product_state(product_aut, step.s)[0], condition_to_str(step.label, d), step.acc)
        for step in run.prefix
    ]
    cycle_trace = [
        (_parse_product_state(product_aut, step.s)[0], condition_to_str(step.label, d), step.acc)
        for step in run.cycle
    ]

    if not cycle_trace:
        raise RuntimeError("Accepting product run has an empty cycle.")

    cycle_entry = cycle_trace[0][0]
    return prefix_trace, cycle_trace, cycle_entry


# ──────────────────────────────────────────────────────────────────────────────
# Word replay
# ──────────────────────────────────────────────────────────────────────────────

def _find_transition(aut, state: int, bdd_cond):
    """Return (dst, acc) for the first transition from *state* matching *bdd_cond*."""
    for t in aut.out(state):
        if buddy.bdd_implies(bdd_cond, t.cond):
            return (t.dst, t.acc)
    return None


def replay_word(aut, word) -> tuple:
    """
    Replay *word* (a spot.twa_word) on *aut* starting from its initial state.

    Returns
    -------
    (prefix_trace, cycle_trace, cycle_entry_state)
        prefix_trace  : list of (src_state, bdd_cond, acc_marks)
        cycle_trace   : list of (src_state, bdd_cond, acc_marks)
        cycle_entry_state : state number where the cycle starts
    Raises RuntimeError if the word cannot be replayed (automaton rejects it).

    Notes
    -----
    A word's cycle (prefix)^ω does not mean the cycle letters form a self-loop in
    the automaton.  Several passes through the cycle letters may be needed before
    any automaton state is revisited.  This function unrolls the cycle letters
    until it finds a real cycle in the automaton (by detecting the first repeated
    automaton state), extending the prefix accordingly.
    """
    current = aut.get_init_state_number()
    prefix_trace = []

    for bdd_cond in word.prefix:
        result = _find_transition(aut, current, bdd_cond)
        if result is None:
            raise RuntimeError(
                f"Cannot replay word prefix: no transition from state {current} "
                f"for condition {bdd_to_str(bdd_cond, aut.get_dict())}"
            )
        dst, acc = result
        prefix_trace.append((current, bdd_cond, acc))
        current = dst

    cycle_letters = list(word.cycle)
    if not cycle_letters:
        raise RuntimeError("Word has an empty cycle – cannot find a real cycle in the automaton.")

    # Unroll cycle letters until an automaton state is revisited.
    # By the pigeonhole principle this is guaranteed within num_states steps.
    num_letters = len(cycle_letters)
    max_steps = aut.num_states() + 1   # pigeonhole bound

    visited_states: dict = {}   # automaton state → index in unrolled_trace
    unrolled_trace: list = []   # (src, bdd_cond, acc)

    for step_count in range(max_steps + 1):
        if current in visited_states:
            # Found the first repeated automaton state → real cycle starts here.
            cycle_start = visited_states[current]
            cycle_entry = current
            cycle_trace = unrolled_trace[cycle_start:]
            extended_prefix = prefix_trace + unrolled_trace[:cycle_start]
            return extended_prefix, cycle_trace, cycle_entry

        visited_states[current] = len(unrolled_trace)
        bdd_cond = cycle_letters[step_count % num_letters]
        result = _find_transition(aut, current, bdd_cond)
        if result is None:
            raise RuntimeError(
                f"Cannot replay word cycle: no transition from state {current} "
                f"for condition {bdd_to_str(bdd_cond, aut.get_dict())}"
            )
        dst, acc = result
        unrolled_trace.append((current, bdd_cond, acc))
        current = dst

    raise RuntimeError(
        f"Could not find a real cycle within {max_steps} steps – "
        "the word may not be accepted by this automaton."
    )


# ──────────────────────────────────────────────────────────────────────────────
# Violation detection
# ──────────────────────────────────────────────────────────────────────────────

def find_violation(kofola_comp, spot_comp, original_aut):
    """
    Determine the direction of the violation and find a witnessing word.

    Returns
    -------
    (direction, word, projected_trace)
        direction : "over_approx"  – kofola over-approximates the correct language
                                     (accepts a word in L(original_aut))
                    "under_approx" – kofola under-approximates the correct language
                                     (rejects a word not in L(original_aut))
        word      : spot.twa_word for the violating word
        projected_trace : projected kofola-side run trace for over-approx, else None
    """
    # Case 1: kofola_comp ∩ L(aut) ≠ ∅ → complement over-approximates.
    # Ask Spot for a witness from the actual intersection language, not from an
    # automaton derived from one accepting run of the product.
    prod = spot.product(kofola_comp, original_aut)
    if not prod.is_empty():
        run = prod.accepting_run()
        if run is not None:
            word = run.as_twa().accepting_word()
            return ("over_approx", word, project_product_run(prod, run))

        word = prod.accepting_word()
        if word is not None:
            return ("over_approx", word, None)

    # Case 2: kofola_comp misses a word from the correct complement →
    # under-approximates.  The witness must come from L(spot_comp) \ L(kofola).
    word = spot_comp.exclusive_word(kofola_comp)
    return ("under_approx", word, None)


# ──────────────────────────────────────────────────────────────────────────────
# Visualization
# ──────────────────────────────────────────────────────────────────────────────

def print_text_visualization(aut, prefix_trace, cycle_trace, cycle_entry, direction: str):
    """Print a human-readable step-by-step visualization of the run."""
    d = aut.get_dict()

    sep = "─" * 80
    print()
    print(sep)
    print("  VIOLATING RUN IN KOFOLA RAW COMPLEMENT")
    print(sep)

    direction_msg = {
        "over_approx": (
            "DIRECTION: kofola complement OVER-APPROXIMATES the correct language\n"
            "  → It accepts a word that belongs to L(original_aut).\n"
            "    A correct complement must NOT accept any word from L(aut)."
        ),
        "under_approx": (
            "DIRECTION: kofola complement UNDER-APPROXIMATES the correct language\n"
            "  → It rejects a word that does NOT belong to L(original_aut).\n"
            "    A correct complement MUST accept all words outside L(aut)."
        ),
    }
    print(direction_msg.get(direction, f"DIRECTION: {direction}"))
    print()

    col_w = 6  # column width for step index

    def fmt_step(idx, kind, src, bdd_cond, acc):
        label   = condition_to_str(bdd_cond, d)
        msname  = get_state_name(aut, src)
        acc_str = str(acc) if acc else "{}"
        return (
            f"  [{idx:{col_w}}] {kind:10s}  "
            f"state={src:<4d}  "
            f"on=[ {label:<20s} ]  "
            f"acc={acc_str:<12s}  "
            f"macrostate=[ {msname} ]"
        )

    step = 0

    if prefix_trace:
        print("  PREFIX:")
        for (src, bdd_cond, acc) in prefix_trace:
            print(fmt_step(step, "prefix", src, bdd_cond, acc))
            step += 1
        print()

    print(f"  CYCLE  (re-enters state {cycle_entry} = {get_state_name(aut, cycle_entry)!r}):")
    for (src, bdd_cond, acc) in cycle_trace:
        print(fmt_step(step, "cycle", src, bdd_cond, acc))
        step += 1

    print()
    print(f"  → cycle loops back to state {cycle_entry}  [{get_state_name(aut, cycle_entry)}]")
    print(sep)


def generate_dot(aut, prefix_trace, cycle_trace, cycle_entry) -> str:
    """
    Return a DOT string containing only the violating run (prefix + cycle loop).

    The output is a minimal automaton built from the run steps:
      - prefix states: filled light-blue
      - cycle states: filled light-green
    Each edge carries the transition label and acceptance marks.
    The cycle back-edge points to the cycle_entry state.
    """
    d = aut.get_dict()
    names = aut.get_state_names()

    def sname(s):
        if names and s < len(names) and names[s]:
            return names[s]
        return str(s)

    def acc_label(acc):
        return f" {{{', '.join(str(c) for c in acc.sets())}}}" if acc.count() else ""

    lines = []
    lines.append('digraph run {')
    lines.append('  rankdir=LR;')
    lines.append('  node [shape=rectangle, style=filled, fontname="monospace", fontsize=9];')
    lines.append('  edge [fontname="monospace", fontsize=9];')
    lines.append('  I [label="", style=invis, width=0];')

    # Collect ordered distinct states in appearance order
    all_steps = prefix_trace + cycle_trace
    seen_states = []
    seen_set = set()
    for (src, _, _) in all_steps:
        if src not in seen_set:
            seen_states.append(src)
            seen_set.add(src)
    # Also add cycle_entry if not already there
    if cycle_entry not in seen_set:
        seen_states.append(cycle_entry)
        seen_set.add(cycle_entry)

    cycle_states  = {src for (src, _, _) in cycle_trace}

    # Emit state nodes
    for s in seen_states:
        color = "lightgreen" if s in cycle_states else "lightblue"
        label = f"s{s}\\n{sname(s)}"
        # Escape backslashes and quotes for DOT
        label = label.replace('"', '\\"')
        lines.append(f'  s{s} [label="{label}", fillcolor={color}];')

    # Initial arrow
    init = aut.get_init_state_number()
    lines.append(f'  I -> s{init};')

    # Emit prefix transitions.
    # The trace order already determines destinations, which keeps projected
    # product runs stable even when the kofola automaton is nondeterministic.
    for i, (src, bdd_cond, acc) in enumerate(prefix_trace):
        if i + 1 < len(prefix_trace):
            dst = prefix_trace[i + 1][0]
        else:
            dst = cycle_trace[0][0] if cycle_trace else cycle_entry
        cond_str = condition_to_str(bdd_cond, d).replace('"', '\\"')
        edge_lbl = f"{cond_str}{acc_label(acc)}"
        lines.append(f'  s{src} -> s{dst} [label="{edge_lbl}", color=blue];')

    # Emit cycle transitions
    for i, (src, bdd_cond, acc) in enumerate(cycle_trace):
        if i + 1 < len(cycle_trace):
            dst = cycle_trace[i + 1][0]
        else:
            dst = cycle_entry  # loop back
        cond_str = condition_to_str(bdd_cond, d).replace('"', '\\"')
        edge_lbl = f"{cond_str}{acc_label(acc)}"
        lines.append(f'  s{src} -> s{dst} [label="{edge_lbl}", color=red];')

    lines.append('}')
    return '\n'.join(lines)


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Debug kofola complement: find and visualize violating accepted runs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Input automaton in HOA format"
    )
    parser.add_argument(
        "--kofola", "-k", default="build/src/kofola",
        help="Path to the kofola binary (default: build/src/kofola)"
    )
    parser.add_argument(
        "--params", "-p",
        default="tela=yes;tela_det_alg=inductive;sd_ind_or_opt=yes",
        help=(
            "kofola --params string. 'raw=yes' is appended automatically so "
            "that macrostate labels are preserved in the HOA output. "
            "(default: tela=yes;tela_det_alg=inductive;sd_ind_or_opt=yes)"
        )
    )
    parser.add_argument(
        "--dot", "-d", default=None, metavar="FILE",
        help="Write a DOT file of the kofola automaton with the run highlighted"
    )
    args = parser.parse_args()

    # ── Ensure raw=yes is present so macrostate labels are kept ───────────────
    params = args.params.rstrip(";")
    if "raw=yes" not in params:
        params += ";raw=yes"

    print(f"Input file : {args.input}")
    print(f"Kofola     : {args.kofola}")
    print(f"Params     : {params}")

    # ── Run kofola ─────────────────────────────────────────────────────────────
    print("\nRunning kofola (raw mode)…")
    try:
        kofola_hoa = run_kofola(args.kofola, args.input, params)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    # ── Load automata ──────────────────────────────────────────────────────────
    try:
        original_aut = spot.automaton(args.input)
    except Exception as e:
        print(f"ERROR loading input automaton: {e}", file=sys.stderr)
        return 2

    try:
        kofola_comp = spot.automaton(kofola_hoa)
    except Exception as e:
        print(f"ERROR parsing kofola HOA output: {e}", file=sys.stderr)
        print("kofola output was:\n", kofola_hoa, file=sys.stderr)
        return 2

    print(f"\nKofola complement : {kofola_comp.num_states()} states, "
          f"acc={kofola_comp.get_acceptance()}")

    # Check whether macrostate labels are present
    names = kofola_comp.get_state_names()
    has_labels = bool(names and any(n for n in names))
    if not has_labels:
        print(
            "WARNING: kofola output has no macrostate state-name labels.\n"
            "         Rebuild kofola or pass raw=yes explicitly."
        )

    # ── Spot complement ────────────────────────────────────────────────────────
    print("Computing Spot complement…")
    spot_comp = spot.complement(original_aut)
    print(f"Spot complement   : {spot_comp.num_states()} states")

    # ── Equivalence check ──────────────────────────────────────────────────────
    print("Checking language equivalence…")
    if spot.are_equivalent(kofola_comp, spot_comp):
        print("\n✓ SUCCESS: kofola complement is language-equivalent to Spot's complement.")
        return 0

    print("\n✗ FAILURE: kofola complement is NOT language-equivalent to Spot's complement!")

    # ── Find violation ─────────────────────────────────────────────────────────
    print("\nSearching for a violating word…")
    try:
        direction, word, projected_trace = find_violation(kofola_comp, spot_comp, original_aut)
    except Exception as e:
        print(f"ERROR finding violation: {e}", file=sys.stderr)
        return 2

    if word is None:
        print("ERROR: Could not construct a violating word (unexpected).", file=sys.stderr)
        return 2

    print(f"Violating word : {word}")

    # ── Replay on kofola automaton ─────────────────────────────────────────────
    if direction == "over_approx" and projected_trace is not None:
        print("Projecting accepting product run to kofola raw automaton…")
        prefix_trace, cycle_trace, cycle_entry = projected_trace
    else:
        print("Replaying word on kofola raw automaton…")
        try:
            prefix_trace, cycle_trace, cycle_entry = replay_word(kofola_comp, word)
        except RuntimeError as e:
            print(f"ERROR replaying word: {e}", file=sys.stderr)
            print(
                "This usually means the kofola complement did not actually accept "
                "the word—check that the direction detection is correct.",
                file=sys.stderr,
            )
            # Fall back: try finding run directly in kofola_comp
            print("Attempting to find any accepting run in kofola_comp as fallback…")
            try:
                run = kofola_comp.accepting_run()
                if run is None:
                    print("No accepting run found (automaton may be empty).", file=sys.stderr)
                    return 2
                # Build artificial word from this run
                print("Visualizing a generic accepting run in kofola complement:")
                prefix_trace = [
                    (kofola_comp.state_number(s.s), s.label, s.acc)
                    for s in run.prefix
                ]
                cycle_trace = [
                    (kofola_comp.state_number(s.s), s.label, s.acc)
                    for s in run.cycle
                ]
                cycle_entry = kofola_comp.state_number(list(run.cycle)[0].s) if run.cycle else 0

                # In fallback mode the BDD columns are already proper BDDs from the run
                print_text_visualization(kofola_comp, prefix_trace, cycle_trace, cycle_entry, direction)
            except Exception as e2:
                print(f"Fallback also failed: {e2}", file=sys.stderr)
            return 2

    # ── Visualize ──────────────────────────────────────────────────────────────
    print_text_visualization(kofola_comp, prefix_trace, cycle_trace, cycle_entry, direction)

    # ── DOT output ────────────────────────────────────────────────────────────
    if args.dot:
        dot_str = generate_dot(kofola_comp, prefix_trace, cycle_trace, cycle_entry)
        with open(args.dot, "w") as fh:
            fh.write(dot_str)
        print(f"\nDOT file written to: {args.dot}")
        print(f"Render with: dot -Tsvg {args.dot} -o run.svg")

    return 1


if __name__ == "__main__":
    sys.exit(main())
