"""Run scripts from within Python.

from runscript import runscript

(o, e, status) = runscript("echo $moo $1", 12, moo='mu')

o -> 'mu 12\\n', e -> '', status -> 0

from runscript import runcmd

f = runcmd("ls -1")     # simple interface returns a list of lines
for l in f:
    ...
"""

from __future__ import print_function

import sys, os
from subprocess import Popen, PIPE

def runscript(script, *args, **kwargs):
    """The script is run by a shell, and the output is returned.
    Positional arguments are set from args, and the kwargs are
    put into the environment to become shell variables. The args
    and kwargs are converted to strings before being passed to
    the script. The result is a pair of strings (stdout, stderr)."""

    shcmd = ['/bin/sh', '-s']
    shcmd.extend(map(str, args))

    for k, v in kwargs.items():
        os.environ[k] = str(v)

    f = Popen(shcmd, stdin = PIPE, stdout = PIPE, stderr = PIPE, universal_newlines = True)
    o, e = f.communicate(script)

    return o, e, f.wait()

def runcmd(*a, **k):
    """Same as runscript() plus default processing of results."""

    o, e, code = runscript(*a, **k)
    if e:
        print(e, file = sys.stderr)
    if code:
        sys.exit(code)

    return o.split('\n') # split the output into lines

