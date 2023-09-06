from __future__ import annotations
import rdkit.TestRunner
import typing
import os
import rdkit.RDConfig
import subprocess
import sys
import time

__all__ = [
    "BUILD_TYPE_ENVVAR",
    "OutputRedirectC",
    "RDConfig",
    "ReportResults",
    "RunScript",
    "RunTest",
    "TEST_FAILED",
    "TEST_PASSED",
    "isDebugBuild",
    "os",
    "redirect_stderr",
    "redirect_stdout",
    "subprocess",
    "sys",
    "time"
]


class OutputRedirectC():
    """
    Context manager which uses low-level file descriptors to suppress
      output to stdout/stderr, optionally redirecting to the named file(s).

      Suppress all output
      with Silence():
        <code>

      Redirect stdout to file
      with OutputRedirectC(stdout='output.txt', mode='w'):
        <code>

      Redirect stderr to file
      with OutputRedirectC(stderr='output.txt', mode='a'):
        <code>
      http://code.activestate.com/recipes/577564-context-manager-for-low-level-redirection-of-stdou/
      >>>

      
    """
    pass
class _RedirectStream():
    _stream = None
    pass
class redirect_stderr(_RedirectStream):
    """
    Context manager for temporarily redirecting stderr to another file.
    """
    _stream = 'stderr'
    pass
class redirect_stdout(_RedirectStream):
    """
    Context manager for temporarily redirecting stdout to another file.

            # How to send help() to stderr
            with redirect_stdout(sys.stderr):
                help(dir)

            # How to write help() to a file
            with open('help.txt', 'w') as f:
                with redirect_stdout(f):
                    help(pow)
        
    """
    _stream = 'stdout'
    pass
BUILD_TYPE_ENVVAR = 'RDKIT_BUILD_TYPE'
TEST_FAILED = -1
TEST_PASSED = 0
