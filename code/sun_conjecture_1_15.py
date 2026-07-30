#!/usr/bin/env python3
"""
Search for the least counterexample above 10^13 to Zhi-Wei Sun's
Conjecture 1.15 from:
    Zhi-Wei Sun, "Mixed Sums of Primes and Other Terms" (2009),
    arXiv:0901.3075.
Conjecture:
    Every integer n > 4 can be written as
        n = p + L_s + C_t,
    where p is an odd prime, L_s is a Lucas number, C_t is a
    Catalan number, and s,t >= 0.
The paper reports verification through 10^13, so the default search
starts at 10^13 + 1.
Dependencies:
    Python 3.10+
    sympy
No external databases or network access are used.
"""
from __future__ import annotations
import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path
from typing import Optional, Sequence
from sympy import isprime
KNOWN_VERIFIED = 10_000_000_000_000
DEFAULT_START = KNOWN_VERIFIED + 1
def lucas_values(limit: int) -> list[int]:
    """
    Return the distinct Lucas numbers <= limit in increasing order.
    L_0 = 2
    L_1 = 1
    L_{n+1} = L_n + L_{n-1}
    Thus the increasing list is:
        1, 2, 3, 4, 7, 11, 18, ...
    """
    if limit < 1:
        return []
    values = [1]
    if limit >= 2:
        values.append(2)
    # L_2=3, L_3=4, L_4=7, ...
    a, b = 3, 4
    while a <= limit:
        values.append(a)
        a, b = b, a + b
    return values
def catalan_values(limit: int) -> list[int]:
    """
    Return the distinct Catalan numbers <= limit in increasing order.
        C_t = binom(2t,t)/(t+1)
    C_0=C_1=1, so duplicate numerical values are removed.
    """
    if limit < 1:
        return []
    values: list[int] = []
    c = 1
    t = 0
    while c <= limit:
        if not values or values[-1] != c:
            values.append(c)
        # Exact integer recurrence:
        #
        # C_{t+1} = C_t * 2(2t+1)/(t+2)
        c = c * 2 * (2 * t + 1) // (t + 2)
        t += 1
    return values
def pair_sums(
    limit: int,
) -> tuple[list[int], list[int], int, int]:
    """
    Construct all distinct offsets L_s + C_t <= limit.
    Returns:
        even_sums
        odd_sums
        number of Lucas values
        number of Catalan values
    Separating by parity avoids calling isprime on even numbers.
    """
    lucas = lucas_values(limit)
    catalans = catalan_values(limit)
    sums = {
        lucas_value + catalan_value
        for lucas_value in lucas
        for catalan_value in catalans
        if lucas_value + catalan_value <= limit
    }
    even_sums = sorted(value for value in sums if value % 2 == 0)
    odd_sums = sorted(value for value in sums if value % 2 == 1)
    return (
        even_sums,
        odd_sums,
        len(lucas),
        len(catalans),
    )
def first_failure_in_chunk(
    task: tuple[int, int, Sequence[int], Sequence[int]],
) -> Optional[int]:
    """
    Search the half-open interval [lo, hi).
    Return the first n in the interval having no representation, or None.
    """
    lo, hi, even_sums, odd_sums = task
    for n in range(lo, hi):
        # p = n - offset must be odd.
        #
        # If n is odd, offset must be even.
        # If n is even, offset must be odd.
        offsets = even_sums if n & 1 else odd_sums
        represented = False
        for offset in offsets:
            p = n - offset
            if p < 3:
                # Offsets are sorted increasingly, so all later p
                # values will be still smaller.
                break
            if isprime(p):
                represented = True
                break
        if not represented:
            return n
    return None
def find_witness(n: int) -> Optional[tuple[int, int, int]]:
    """
    Independently search actual Lucas/Catalan pairs for a representation.
    This deliberately does not reuse the precomputed pair-sum lists.
    It is used to recheck any alleged counterexample.
    Returns:
        (p, L_s, C_t), if a representation exists;
        None otherwise.
    """
    lucas = lucas_values(n - 4)
    catalans = catalan_values(n - 4)
    for lucas_value in lucas:
        for catalan_value in catalans:
            p = n - lucas_value - catalan_value
            if p < 3:
                continue
            if (p & 1) and isprime(p):
                return p, lucas_value, catalan_value
    return None
def split_interval(
    lo: int,
    hi: int,
    pieces: int,
) -> list[tuple[int, int]]:
    """
    Split [lo,hi) into contiguous nonempty intervals.
    """
    length = hi - lo
    pieces = max(1, min(pieces, length))
    base, extra = divmod(length, pieces)
    result: list[tuple[int, int]] = []
    cursor = lo
    for i in range(pieces):
        width = base + (1 if i < extra else 0)
        result.append((cursor, cursor + width))
        cursor += width
    return result
def save_checkpoint(
    path: Path,
    original_start: int,
    next_start: int,
) -> None:
    """
    Atomically save the next untested integer.
    """
    payload = {
        "conjecture": (
            "Sun Conjecture 1.15: "
            "n = odd prime + Lucas number + Catalan number"
        ),
        "original_start": original_start,
        "last_completed": next_start - 1,
        "next_start": next_start,
        "saved_unix_time": time.time(),
    }
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    temporary.replace(path)
def load_checkpoint(path: Path) -> int:
    """
    Read and validate a checkpoint.
    """
    payload = json.loads(path.read_text(encoding="utf-8"))
    next_start = int(payload["next_start"])
    if next_start < 5:
        raise ValueError("Invalid checkpoint: next_start is below 5")
    return next_start
def search(
    start: int,
    end: Optional[int],
    block_size: int,
    workers: int,
    checkpoint: Optional[Path],
) -> Optional[int]:
    """
    Search integers in strictly increasing blocks.
    Arguments:
        start:
            First integer to test.
        end:
            Last integer to test, inclusive.
            None means continue until interrupted or a counterexample is found.
        block_size:
            Number of integers per checkpoint block.
        workers:
            Number of worker processes.
        checkpoint:
            Checkpoint path, or None to disable checkpoints.
    Returns:
        The least counterexample encountered, or None if the finite
        interval ends without one.
    """
    start = max(5, start)
    if end is not None and end < start:
        raise ValueError("--end must be greater than or equal to --start")
    if block_size <= 0:
        raise ValueError("--block-size must be positive")
    if workers <= 0:
        raise ValueError("--workers must be positive")
    original_start = start
    lo = start
    total_started = time.perf_counter()
    pool = None
    try:
        if workers > 1:
            pool = mp.Pool(processes=workers)
        while end is None or lo <= end:
            hi = lo + block_size
            if end is not None:
                hi = min(hi, end + 1)
            # For n <= hi-1 and p >= 3, a useful offset satisfies
            #
            #     offset = n-p <= (hi-1)-3 = hi-4.
            even_sums, odd_sums, lucas_count, catalan_count = (
                pair_sums(hi - 4)
            )
            intervals = split_interval(lo, hi, workers)
            tasks = [
                (
                    chunk_lo,
                    chunk_hi,
                    even_sums,
                    odd_sums,
                )
                for chunk_lo, chunk_hi in intervals
            ]
            block_started = time.perf_counter()
            if pool is None:
                results = [first_failure_in_chunk(tasks[0])]
            else:
                results = pool.map(first_failure_in_chunk, tasks)
            failures = [
                value
                for value in results
                if value is not None
            ]
            if failures:
                # Because the worker intervals are contiguous and each worker
                # returns its first failure, the minimum is the first failure
                # in the complete block.
                candidate = min(failures)
                # Independent pair-by-pair verification.
                witness = find_witness(candidate)
                if witness is not None:
                    p, lucas_value, catalan_value = witness
                    raise RuntimeError(
                        "Internal consistency failure: the alleged "
                        "counterexample actually has the representation "
                        f"{candidate} = {p} + {lucas_value} "
                        f"+ {catalan_value}"
                    )
                elapsed = time.perf_counter() - total_started
                print()
                print("*** COUNTEREXAMPLE FOUND ***")
                print(f"n = {candidate}")
                print(
                    "Exhaustively searched "
                    f"[{original_start:,}, {candidate:,}]"
                )
                print(f"Elapsed time: {elapsed:.2f} seconds")
                print(
                    "Independent pair-by-pair recheck found "
                    "no representation."
                )
                return candidate
            elapsed = time.perf_counter() - block_started
            rate = (hi - lo) / elapsed if elapsed else float("inf")
            print(
                f"covered [{lo:,}, {hi - 1:,}] | "
                f"Lucas={lucas_count}, "
                f"Catalan={catalan_count}, "
                f"distinct offsets="
                f"{len(even_sums) + len(odd_sums)} | "
                f"{elapsed:.2f}s | "
                f"{rate:,.0f} n/s",
                flush=True,
            )
            lo = hi
            if checkpoint is not None:
                save_checkpoint(
                    checkpoint,
                    original_start,
                    lo,
                )
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    assert end is not None
    print(
        f"No counterexample in "
        f"[{original_start:,}, {end:,}]."
    )
    return None
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Search for a counterexample to Sun's "
            "odd-prime + Lucas + Catalan conjecture."
        )
    )
    parser.add_argument(
        "--start",
        type=int,
        default=DEFAULT_START,
        help=(
            "first integer to test "
            f"(default: {DEFAULT_START})"
        ),
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help=(
            "last integer to test, inclusive; "
            "omit to continue indefinitely"
        ),
    )
    parser.add_argument(
        "--block-size",
        type=int,
        default=1_000_000,
        help=(
            "integers per checkpoint block "
            "(default: 1000000)"
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 1),
        help=(
            "worker processes "
            "(default: CPU count minus one)"
        ),
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=Path(
            "sun_lucas_catalan_checkpoint.json"
        ),
        help="checkpoint JSON path",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "resume from --checkpoint; "
            "this overrides --start"
        ),
    )
    parser.add_argument(
        "--no-checkpoint",
        action="store_true",
        help="disable checkpoint writes",
    )
    return parser.parse_args()
def main() -> None:
    mp.freeze_support()
    args = parse_args()
    checkpoint: Optional[Path]
    if args.no_checkpoint:
        checkpoint = None
    else:
        checkpoint = args.checkpoint
    start = args.start
    if args.resume:
        if checkpoint is None:
            raise SystemExit(
                "--resume cannot be combined with --no-checkpoint"
            )
        if not checkpoint.exists():
            raise SystemExit(
                f"Checkpoint not found: {checkpoint}"
            )
        start = load_checkpoint(checkpoint)
        print(
            f"Resuming at n={start:,} "
            f"from {checkpoint}"
        )
    search(
        start=start,
        end=args.end,
        block_size=args.block_size,
        workers=args.workers,
        checkpoint=checkpoint,
    )
if __name__ == "__main__":
    main()
