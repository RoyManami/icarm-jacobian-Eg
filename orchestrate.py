#!/usr/bin/env python3
"""
orchestrate.py -- target list + post-processing for the quotient harvest.

Usage:
    python3 orchestrate.py targets.json            # show the gap targets
    python3 orchestrate.py targets.json out1.txt out2.txt   # ingest Magma output

Each Magma output line is "<genus>\t<tag>\t<source...>".
This script collects the genera reported, intersects them with the official
gap list, and prints which Ekedahl-Serre gaps the harvest closed.
"""
import json, sys

def load_targets(path):
    t = json.load(open(path))
    gaps = set(t["gaps_56_240"]) | set(t["gaps_241_1297"])
    return t, gaps

def parse_harvest(path):
    out = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        try:
            g = int(parts[0])
        except ValueError:
            continue
        out[g] = parts[2] if len(parts) > 2 else ""
    return out

def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    t, gaps = load_targets(sys.argv[1])
    if len(sys.argv) == 2:
        print(f"Smallest unrealized genus: {t['smallest_gap']}")
        print(f"Gaps in [56,240]  ({len(t['gaps_56_240'])}): {t['gaps_56_240']}")
        print(f"Gaps in [241,1297]({len(t['gaps_241_1297'])}) -- first 30: "
              f"{t['gaps_241_1297'][:30]} ...")
        print(f"\nCompletely decomposable X_0(N) sources, N<=1200: "
              f"{len(t['good_levels_le1200'])} levels")
        print(f"High-level sources (N>=240, outside the N<240 census): "
              f"{ {k:v for k,v in t['high_sources_ge240'].items()} }")
        return
    filled = {}
    for f in sys.argv[2:]:
        for g, src in parse_harvest(f).items():
            if g in gaps:
                filled.setdefault(g, src)
    print(f"Gaps closed by harvest: {len(filled)}")
    for g in sorted(filled):
        print(f"  genus {g:4d}   <-   {filled[g]}")
    remaining = sorted(gaps - set(filled))
    print(f"\nStill open after harvest (first 20): {remaining[:20]} ...")

if __name__ == "__main__":
    main()
