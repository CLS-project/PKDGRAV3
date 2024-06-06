# distutils: language = c++
# cython: always_allow_keywords=True

import cython
from libcpp cimport bool

@cython.locals(keys=cython.dict)
cdef public bool is_PKDGRAV_imported(dict keys):
  import inspect
  for k,v in keys.items():
    if inspect.ismodule(v):
      name=inspect.getmodule(v).__name__
      if name == "PKDGRAV":
        return True
  return False