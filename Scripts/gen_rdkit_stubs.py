#!/usr/bin/env python

"""
Script to generate Python stubs for RDKit.

This script is invoked as part of the build process
by setting the CMake switch RDK_INSTALL_PYTHON_STUBS=ON.
If you decide to run this script outside the build process,
make sure that the RDKit Python modules for which stubs are
to be generated are the *first* RDKit modules available in
sys.path; otherwise, stubs will not be generated for the
intended RDKit version.

Usage:
./Scripts/gen_rdkit_stubs.py [output_dirs; defaults to $PWD]

Usage example:
$ cd $RDBASE
$ ./Scripts/gen_rdkit_stubs.py
$ cp -R rdkit-stubs $CONDA_PREFIX/lib/python3.*/site-packages

The scripts creates an rdkit-stubs directory in each
directory in output_dirs.
Warnings printed to console can be safely ignored.
"""

import sys
import os
import re
import tempfile
import shutil
from pathlib import Path
import pybind11_stubgen

PY_SIGNATURE_ARG_REGEX = re.compile(r"^\((\S+)\)([^\=]+)\=?(.*)?$")
DEF_REGEX = re.compile(r"^([^(]+)(\s*\(\s*).*(\s*\)\s*->\s*)[^:]+:(\s*)$")
PURGE_CPP_OBJECT_ANGLE_BRACKETS = re.compile(r"^(.*)<(\S*\.)?(\S+)\s*object\s*at\s*\S+\s*>(.*)$")
PURGE_OPEN_SQUARE_BRACKET = re.compile(r"\[(?!\])")
PURGE_CLOSE_SQUARE_BRACKET = re.compile(r"(?<!\[)\]")
IMPORT_MODULES = re.compile(r"^\s*import\s+(.*)$")
FROM_IMPORT_MODULES = re.compile(r"^\s*from\s+(\S+)\s+import\s+(.*)$")
RDKIT_MODULE_NAME = "rdkit"
RDKIT_RDCONFIG = "RDConfig.py"
PROTECTIONS = {
    "[": "__OPEN_SQUARE_BRACKET_TAG__",
    "]": "__CLOSE_SQUARE_BRACKET_TAG__",
    "=": "__EQUALS_TAG__"
}
INIT_PY = "__init__.py"

def parse_modules_to_set(modules):
    return {m.strip() for m in modules.split(",")}

def concat_parent_child_module(parent_module, child_module):
    if child_module != "*":
        parent_module += "." + child_module
    return parent_module

def clear_stubs(outer_dir):
    for entry in os.listdir(outer_dir):
        if entry == "CMakeLists.txt":
            continue
        entry = os.path.join(outer_dir, entry)
        if os.path.isdir(entry):
            shutil.rmtree(entry)
        else:
            os.remove(entry)
            
def copy_stubs(src_dir, dst_dir):
    for entry in os.listdir(src_dir):
        src_entry = os.path.join(src_dir, entry)
        dst_entry = os.path.join(dst_dir, entry)
        if os.path.isdir(src_entry):
            shutil.copytree(src_entry, dst_entry)
        else:
            shutil.copy(src_entry, dst_entry)

def generate_stubs(site_packages_path, output_dirs):
    modules = set()
    exclusions_cache = {}
    for p in (
        # order is important here
        sorted(site_packages_path.joinpath(RDKIT_MODULE_NAME).rglob("*.pyd")) +
        sorted(site_packages_path.joinpath(RDKIT_MODULE_NAME).rglob("*.so")) +
        sorted(site_packages_path.joinpath(RDKIT_MODULE_NAME).rglob("*.py"))
    ):
        module_file = p.relative_to(site_packages_path)
        d = p.parent
        if module_file.name == INIT_PY:
            m = str(d.relative_to(site_packages_path)).replace("/", ".")
            modules.add(m)
        else:
            noext, ext = os.path.splitext(str(module_file))
            noext = noext.replace("/", ".")
            init_py_file = str(d.joinpath(INIT_PY))
            if os.path.exists(init_py_file):
                exclusions = exclusions_cache.get(init_py_file, None)
                if exclusions is None:
                    exclusions = set()
                    with open(init_py_file, "r") as hnd:
                        for line in hnd:
                            import_modules_match = IMPORT_MODULES.match(line)
                            if import_modules_match:
                                comma_sep_modules = import_modules_match.group(1)
                                exclusions.update(parse_modules_to_set(comma_sep_modules))
                                continue
                            from_import_modules_match = FROM_IMPORT_MODULES.match(line)
                            if from_import_modules_match:
                                parent_module = from_import_modules_match.group(1)
                                comma_sep_modules = from_import_modules_match.group(2)
                                exclusions.update({concat_parent_child_module(parent_module, child_module)
                                                   for child_module in parse_modules_to_set(comma_sep_modules)})
                    exclusions_cache[init_py_file] = exclusions
            else:
                exclusions = set()
            if (ext == ".py") or (noext not in exclusions):
                modules.add(noext)
    outer_dirs = []
    for output_dir in output_dirs:
        outer_dir = os.path.join(output_dir, "rdkit-stubs")
        outer_dirs.append(outer_dir)
        if os.path.isdir(outer_dir):
            clear_stubs(outer_dir)
        elif os.path.isfile(outer_dir):
            os.remove(outer_dir)
        if not os.path.isdir(outer_dir):
            os.makedirs(outer_dir)
    with tempfile.TemporaryDirectory() as tempdir:
        for m in sorted(modules):
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
                print(str(e))
                if isinstance(e, AssertionError):
                    raise
        src_dir = os.path.join(tempdir, "rdkit")
        for outer_dir in outer_dirs:
            copy_stubs(src_dir, outer_dir)

def protect_quoted_square_brackets_and_equals(arg):
    open_quote = False
    protected_arg = ""
    for c in arg:
        if c == "'":
            open_quote ^= True
        elif open_quote:
            replacement = PROTECTIONS.get(c, None)
            if replacement is not None:
                c = replacement
        protected_arg += c
    return protected_arg

def deprotect_quoted_square_brackets_and_equals(arg):
    for k, v in PROTECTIONS.items():
        arg = arg.replace(v, k)
    return arg

def process_py_signature_arg(arg):
    arg = PURGE_CPP_OBJECT_ANGLE_BRACKETS.sub(r"\1\3()\4", arg)
    arg = protect_quoted_square_brackets_and_equals(arg)
    arg = PURGE_OPEN_SQUARE_BRACKET.sub("", arg)
    arg = PURGE_CLOSE_SQUARE_BRACKET.sub("", arg)
    arg = arg.replace(" ", "")
    py_signature_arg_match = PY_SIGNATURE_ARG_REGEX.match(arg)
    assert py_signature_arg_match
    py_signature_arg_type = py_signature_arg_match.group(1)
    py_signature_arg_name = py_signature_arg_match.group(2)
    py_signature_arg_default = deprotect_quoted_square_brackets_and_equals(py_signature_arg_match.group(3))
    res = f"{py_signature_arg_name}: {py_signature_arg_type}"
    if py_signature_arg_default:
        res += f" = {py_signature_arg_default}"
    return res

def process_src_line(src_line):
    def_match = DEF_REGEX.match(src_line)
    if def_match:
        func_name = def_match.group(1)
        func_open_bracket = def_match.group(2)
        func_end_bracket_and_arrow = def_match.group(3)
        func_colon_to_end = def_match.group(4)
        py_signature_regex = re.compile(r"^\s*" + func_name + r"\s*\((.*)\)\s*->\s*(\S+)\s*:\s*$")
        py_signature_match = py_signature_regex.match(src_line)
        if py_signature_match:
            py_signature_args = py_signature_match.group(1)
            if py_signature_args:
                py_signature_args = py_signature_args.split(", ")
            else:
                py_signature_args = []
            py_signature_ret = py_signature_match.group(2)
            processed_args = ", ".join(process_py_signature_arg(py_signature_arg) for py_signature_arg in py_signature_args)
            src_line = f"{func_name}{func_open_bracket}{processed_args}{func_end_bracket_and_arrow}{py_signature_ret}{func_colon_to_end}"
    return src_line

def preprocess_docstring(docstring):
    src_lines = docstring.split("\n")
    return "\n".join(map(process_src_line, src_lines))


if __name__ == "__main__":
    pybind11_stubgen.function_docstring_preprocessing_hooks.append(preprocess_docstring)
    site_packages_path = None
    for path_entry in sys.path:
        if not path_entry:
            continue
        rdkit_path = os.path.join(path_entry, RDKIT_MODULE_NAME)
        if os.path.isdir(rdkit_path) and os.path.isfile(os.path.join(rdkit_path, RDKIT_RDCONFIG)):
            site_packages_path = path_entry
            break
    if site_packages_path is None:
        raise ValueError("Failed to find rdkit in PYTHONPATH")
    site_packages_path = Path(site_packages_path)
    try:
        output_dirs = sys.argv[1:]
    except IndexError:
        output_dirs = [os.getcwd()]
    generate_stubs(site_packages_path, output_dirs)
