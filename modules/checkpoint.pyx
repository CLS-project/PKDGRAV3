# distutils: language = c++
# cython: always_allow_keywords=True

import cython

@cython.locals(keys=cython.dict)
cdef public void print_imports(const char * filename, dict keys):
  import inspect
  with open(filename,"w") as fp:
    for k,v in keys.items():
      if k[0] == '_':
        pass
      elif inspect.ismodule(v):
        name=inspect.getmodule(v).__name__
        if k == name:
          print(f'import {name}',file=fp)
        else:
          print(f'import {name} as {k}',file=fp)
      else:
        module = inspect.getmodule(v)
        if module and module.__name__ != '__main__':
          if hasattr(v,'__name__') and v.__name__ != k:
            print(f'from {module.__name__} import {v.__name__} as {k}',file=fp)
          else:
            print(f'from {module.__name__} import {k}',file=fp)