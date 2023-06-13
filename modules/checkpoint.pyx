# distutils: language = c++
# cython: always_allow_keywords=True

import cython
from libcpp cimport bool

@cython.locals(keys=cython.dict)
cdef public void print_imports(const char * filename, dict keys):
  import inspect
  with open(filename,"w") as fp:
    print("from MASTER import MSR",file=fp)
    print("from argparse import Namespace",file=fp)
    for k,v in keys.items():
      if k[0] == '_':
        pass
      elif k in ['MSR','Namespace']:
        pass
      elif inspect.ismodule(v):
        name=inspect.getmodule(v).__name__
        if k == name:
          print(f'import {name}',file=fp)
        else:
          print(f'import {name} as {k}',file=fp)
      elif inspect.isclass(v):
        module = inspect.getmodule(v)
        if module and module.__name__ != '__main__':
          if hasattr(v,'__name__') and v.__name__ != k:
            print(f'from {module.__name__} import {v.__name__} as {k}',file=fp)
          else:
            print(f'from {module.__name__} import {k}',file=fp)

@cython.locals(keys=cython.dict)
cdef public bool is_PKDGRAV_imported(dict keys):
  import inspect
  for k,v in keys.items():
    if inspect.ismodule(v):
      name=inspect.getmodule(v).__name__
      if name == "PKDGRAV":
        return True
  return False