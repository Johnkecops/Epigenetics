#!/usr/bin/env python3
"""
Module: test_deg_table
Purpose: Regression check for EpigenomicsPipeline.deg_table crashing on
         single_stage results (pandas DataFrame truthiness is ambiguous,
         so `df_or_none1 or df_or_none2` raises ValueError whenever the
         first operand is a non-None DataFrame).
Author: Dr. Arli Aditya Parikesit
Date: 2026
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from pipeline.pipeline import EpigenomicsPipeline


def test_deg_table_single_stage():
    """single_stage run populates 'degs' -> deg_table must return it, not crash."""
    p = EpigenomicsPipeline()
    degs = pd.DataFrame({"log2FoldChange": [1.0, -2.0]})
    p.results = {"degs": degs, "summary": {}}
    out = p.deg_table
    assert out is degs, f"expected degs table, got {out!r}"


def test_deg_table_dual_stage():
    """dual_stage run has no 'degs' key -> falls back to 'degs_stage1'."""
    p = EpigenomicsPipeline()
    degs1 = pd.DataFrame({"log2FoldChange": [3.0]})
    p.results = {"degs_stage1": degs1, "summary": {}}
    out = p.deg_table
    assert out is degs1, f"expected degs_stage1 table, got {out!r}"


def test_deg_table_no_results():
    p = EpigenomicsPipeline()
    assert p.deg_table is None


if __name__ == "__main__":
    test_deg_table_single_stage()
    test_deg_table_dual_stage()
    test_deg_table_no_results()
    print("OK - deg_table regression checks pass")
