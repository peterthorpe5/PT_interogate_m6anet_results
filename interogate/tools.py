#!/usr/bin/env python
#
# tools.py


import gzip
import os
import gzip
import sys
import subprocess

class NotExecutableError(Exception):
    """Exception raised when expected executable is not executable"""
    def __init__(self, message):
        self.message = message


def is_exe(filename):
    """Returns True if path is to an executable file"""
    if os.path.isfile(filename) and os.access(filename, os.X_OK):
        return True
    else:
        try:
            exefile = check_output(["which", filename]).strip()
        except CalledProcessError:
            raise NotExecutableError("{0} does not exist".format(filename))
    return os.path.isfile(exefile) and os.access(exefile, os.X_OK)


def return_real_line(line):
    """function to return only true lines,
    no comments, or blank lines."""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):  # dont want comment lines
        return False
    if "\t" in line:
        cluster_line = line.rstrip("\n").split("\t")
    else:
        # different clustering program?
        line = line.rstrip("\n").split()
    return line


