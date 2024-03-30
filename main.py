import sys
import hy
from systems_sandbox import *
from ltiza import *
hy.macros.require("systems_sandbox",None,"ALL")
hy.macros.require("ltiza",None,"ALL")

if __name__=='__main__':
    if len(sys.argv) < 2:
        hy.REPL(locals=locals()).run()
    else:
        hy.importer.runhy.run_path(sys.argv[1])
