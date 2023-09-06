from __future__ import annotations
import rdkit
import typing
import logging
import sys

__all__ = [
    "log_handler",
    "logger",
    "logging",
    "rdBase",
    "sys"
]


__version__ = '2023.09.1pre'
log_handler: logging.StreamHandler # value = <StreamHandler <stderr> (NOTSET)>
logger: logging.Logger # value = <Logger rdkit (WARNING)>
