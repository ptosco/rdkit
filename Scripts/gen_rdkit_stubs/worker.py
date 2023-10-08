import sys
import os
import pybind11_stubgen
import importlib


"""
worker script
1st param is the temporary directory
2nd param is the module name
"""
if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    gen_rdkit_stubs = importlib.import_module("gen_rdkit_stubs")
    pybind11_stubgen.function_docstring_preprocessing_hooks.append(gen_rdkit_stubs.preprocess_docstring)
    tempdir = sys.argv[1]
    m = sys.argv[2]
    outer_dirs = sys.argv[3:]
    try:
        pybind11_stubgen.main(args=("--root-module-suffix",
                                    "",
                                    "--no-setup-py",
                                    "--ignore-invalid",
                                    "all",
                                    "-o",
                                    tempdir,
                                    m))
    except Exception as e:
        if isinstance(e, AssertionError):
            raise
        else:
            print(str(e))
    gen_rdkit_stubs.copy_stubs(os.path.join(tempdir, *m.split(".")), outer_dirs)
