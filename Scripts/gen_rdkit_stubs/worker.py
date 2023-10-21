import sys
import os
import pybind11_stubgen
from pybind11_stubgen.parser.mixins.parse import ExtractSignaturesFromPybind11Docstrings
import importlib


"""
worker script
1st param is the temporary directory
2nd param is the module name
"""
if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    gen_rdkit_stubs = importlib.import_module("gen_rdkit_stubs")

    parse_function_docstring_orig = ExtractSignaturesFromPybind11Docstrings.parse_function_docstring
    handle_property_orig = ExtractSignaturesFromPybind11Docstrings.handle_property

    def handle_property_patched(self, path, prop, **kwargs):
        if hasattr(prop, "fget") and hasattr(prop.fget, "__doc__") and prop.fget.__doc__ is None:
            prop = property(getattr(prop, "fget", None),
                            getattr(prop, "fset", None),
                            getattr(prop, "fdel", None),
                            "")
        return handle_property_orig(self, path, prop, **kwargs)

    def parse_function_docstring_patched(self, func_name, doc_lines, **kwargs):
        doc_lines = gen_rdkit_stubs.preprocess_doc_lines(doc_lines)
        i = 0
        for i, doc_line in enumerate(doc_lines):
            if doc_line:
                break
        for _ in range(i):
            doc_lines.pop(0)
        return parse_function_docstring_orig(self, func_name, doc_lines, **kwargs)

    ExtractSignaturesFromPybind11Docstrings.parse_function_docstring = parse_function_docstring_patched
    ExtractSignaturesFromPybind11Docstrings.handle_property = handle_property_patched

    tempdir = sys.argv[1]
    m = sys.argv[2]
    outer_dirs = sys.argv[3:]
    stored_argv = list(sys.argv)
    try:
        sys.argv = ["",
                    "--root-suffix",
                    "",
                    "--ignore-all-errors",
                    "-o",
                    tempdir,
                    m]
        pybind11_stubgen.main()
    except Exception as e:
        if isinstance(e, AssertionError):
            raise
        else:
            print(str(e))
    finally:
        sys.argv = stored_argv
    gen_rdkit_stubs.copy_stubs(os.path.join(tempdir, *m.split(".")), outer_dirs)
