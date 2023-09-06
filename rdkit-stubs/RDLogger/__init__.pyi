from __future__ import annotations
import rdkit.RDLogger
import typing
import sys
import traceback

__all__ = [
    "AttachFileToLog",
    "CRITICAL",
    "DEBUG",
    "DisableLog",
    "ERROR",
    "EnableLog",
    "INFO",
    "LogMessage",
    "WARNING",
    "logger",
    "sys",
    "traceback"
]


class logger():
    pass
def AttachFileToLog( spec: str, filename: str, delay: int = 100) -> None:
    """
    AttachFileToLog( spec: str, filename: str, delay: int = 100) -> None
        Causes the log to write to a file

        C++ signature :
            void AttachFileToLog(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=100])
    """
def DisableLog( arg1: str) -> None:
    """
    DisableLog( arg1: str) -> None

        C++ signature :
            void DisableLog(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def EnableLog( arg1: str) -> None:
    """
    EnableLog( arg1: str) -> None

        C++ signature :
            void EnableLog(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def LogMessage( arg1: str, arg2: str) -> None:
    """
    LogMessage( arg1: str, arg2: str) -> None
        Log a message to any rdApp.* log

        C++ signature :
            void LogMessage(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
CRITICAL = 4
DEBUG = 0
ERROR = 3
INFO = 1
WARNING = 2
_levels = ['rdApp.debug', 'rdApp.info', 'rdApp.warning', 'rdApp.error']
