# Read AFNI .HEAD files. The blocks are returned as a dict, using the
# block's name as a key.

from __future__ import print_function

import re, os.path, sys

def afni_open(filename, which):
    name = os.path.splitext(filename)[0]
    try:
        f = open("%s.%s" % (name, which), 'rb')
    except IOError:
        try:
            f = open("%s+orig.%s" % (name, which), 'rb')
        except IOError:
            print("can't open %s" % filename, file = sys.stderr)
            sys.exit(1)
    return f

# Reads are in binary mode to combat binary garbage that doesn't decode
# into utf-8. Such strings sometimes occur in .HEAD files. So we need to
# convert bytes back into strings; using the latin-1 encoding does this.

def bytes2str(l):
    return l.decode('latin-1')

def getline(f):
    return bytes2str(f.readline())

class scanner:
    def __init__(self):
        pass

    def __call__(self, p, l):
        self.match = p.search(l)
        return self.match

    def group(self, x):
        return self.match.group(x)

def afni_header_read(filename):
    f = afni_open(filename, 'HEAD')
    search = scanner()

    p1 = re.compile(r"type  *= (?P<type>.*)")
    p2 = re.compile(r"name  *= (?P<name>.*)")
    p3 = re.compile(r"count  *= (?P<count>.*)")

    d = {}
    l = getline(f)
    while len(l) != 0:
        if search(p1, l):
            thd_type = search.group('type')
        if search(p2, l):
            thd_name = search.group('name')
        if search(p3, l):
            thd_count = int(search.group('count'))
            d[thd_name] = _decode_block(thd_type, thd_count, f)
        l = getline(f)
    return d

def _decode_block(thd_type, thd_count, f):
    if thd_type == 'string-attribute':
        s = bytes2str(f.read(thd_count + 1))
        if s[0] != "'" or len(s) != thd_count + 1:
            raise RuntimeError("afni header string-attribute")
        s = s[1:].strip()
        s = s.replace('\\n', '\n').split('~')[:-1]
        if len(s) == 1: s = s[0]
        return s

    if thd_type == 'float-attribute':
        l = []
        while len(l) != thd_count:
            x = f.readline().strip().split()
            x = map(float, x)
            l.extend(x)
        return l

    if thd_type == 'integer-attribute':
        l = []
        while len(l) != thd_count:
            x = f.readline().strip().split()
            x = map(int, x)
            l.extend(x)
        return l
