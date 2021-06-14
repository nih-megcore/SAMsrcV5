# Some useful stuff.

import sys, os, getopt
try:
    import subprocess
except:
    import popen2

__scriptname = os.path.basename(sys.argv[0])

def usage(s):
    global __usage
    __usage = s

def msg(s):
    sys.stderr.write(s)

def printerror(s):
    msg("%s: %s\n" % (__scriptname, s))

def printusage():
    msg("usage: %s %s\n" % (__scriptname, __usage))

def parseargs(*args):
    """usage: parseargs(opt, [lopt]) where opt is a string
    that looks like "a:bc:" and lopt (if present) is a list of
    long options without the --"""

    try:
        optlist, args = getopt.getopt(sys.argv[1:], *args)
    except Exception as msg:
        printerror(msg)
        printusage()
        sys.exit(1)
    return optlist, args

def run(cmd, raw = False):
    try:
        p = subprocess.Popen(cmd, shell = True, close_fds = True,
                             stdout = subprocess.PIPE)
        pr = p.stdout
    except:
        (pr, pw) = popen2.popen2(cmd)
        pw.close()
    if raw:
        r = pr.read()
    else:
        r = pr.readlines()
    pr.close()
    return r

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

    f = subprocess.Popen(shcmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, universal_newlines = True)
    o, e = f.communicate(script)

    return o, e, f.wait()

