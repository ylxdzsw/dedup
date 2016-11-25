from os import listdir
from cffi import FFI

modules = set(x[:-2] for x in listdir('.') if x[-2:] == '.c')
print modules
for module in modules:
    C = FFI()
    C.cdef(open(module+'.h').read())
    C.set_source(module, open(module+'.c').read())
    C.compile(verbose=True)
