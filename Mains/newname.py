#! /usr/bin/env python

# Change the names of the programs we copy to $(BDIR)

import sys

prefix = ['SAM', 'MUT', 'ROI']

if len(sys.argv) < 2:
    print('usage: %s name' % sys.argv[0])

def newname(name):

    for p in prefix:
        if name.startswith(p):
            n = len(p)
            s = p.lower()
            return "%s_%s" % (s, name[n:])

    # Not found, keep the old name

    return name

for name in sys.argv[1:]:
    print(newname(name))
