#!/usr/bin/env python3
"""
TreeHarmonizer wrapper script.

This wrapper suppresses SyntaxWarnings from ete3 (an older library with
invalid escape sequences) before importing the main module.

Usage: python run_th.py [args...]
   or: ./run_th.py [args...]  (if made executable)
"""
import sys
import warnings

# Suppress SyntaxWarnings before any imports that might trigger them
warnings.filterwarnings("ignore", category=SyntaxWarning)

# Now import and run the main module
from th_main import main

if __name__ == "__main__":
    main()
