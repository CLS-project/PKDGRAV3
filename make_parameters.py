from sys import argv,exit,version_info
if version_info >= (3,11):
  import tomllib as tl
else:
  import tomli as tl

def emit_str(section,parameters_h,parameters_pxi):
  for k,v in section:
    if isinstance(v,(str,int,float,bool)): continue
    if 'default' in v:
      print(f'    static constexpr auto str_{k:{w}} = "{k}";',file=parameters_h)
      print(f'        const char *str_{k}',file=parameters_pxi)
    else:
      emit_str(v.items(),parameters_h,parameters_pxi)

def get_types(default):
  if isinstance(default,float):
    return 'double','double'
  elif isinstance(default,bool):
    return 'bool','bool'
  elif isinstance(default,int):
    return 'int64_t','int64_t'
  elif isinstance(default,str):
    return 'std::string_view','string'
  else:
    return None, None

def emit_proto(section,parameters_h,parameters_pxi):
  for k,v in section:
    if isinstance(v,(str,int,float,bool)): continue
    if 'default' in v:
      default = v['default']
      help = v['help']
      size = v['size'] if 'size' in v else None
      print(f'    /// {help}',file=parameters_h)
      name = f'str_{k}'
      c_type,i_type = get_types(default)
      if c_type is None:
        if (len(default) == 0):
          c_type = 'PyObject*'
          i_type = 'object'
        else:
          item_type = get_types(default[0])[0]
          # name = f'get<PyObject*>(str_{k})'
          c_type = f'{item_type},{len(default)}'
          i_type = f'TinyVector[{item_type},BLITZ{len(default)}]'
      print(f'    auto get_{k:{w}}() const {{ return get<{c_type}>({name}); }}',file=parameters_h)
      print(f'        {i_type:<7} get_{k}()',file=parameters_pxi)
      print(f'    bool has_{k:{w}}() const {{ return has(str_{k}); }}',file=parameters_h)
      print(f'        bool    has_{k}()',file=parameters_pxi)
      print(f'    void set_{k:{w}}({c_type} value) {{ set<{c_type}>({name},value); }}',file=parameters_h)
      print(f'        void    set_{k}({i_type} value)',file=parameters_pxi)
    else:
      emit_proto(v.items(),parameters_h,parameters_pxi)

if len(argv) <= 3: exit('Usage: {} toml parameters.h parameters.pxi'.format(argv[0]))

with open(argv[1],"rb") as fp:
  f=tl.load(fp)


with open(argv[2], 'w') as parameters_h,open(argv[3], 'w') as parameters_pxi:
  print('''cdef extern from "pkd_parameters.h":
    cdef cppclass pkd_parameters:
        void set(const char *name,double value)
        void set(const char *name,uint64_t value)
        void set(const char *name,int64_t value)
        void prm2ppy(prmContext *prm)
        bool ppy2prm(prmContext *prm)
        bool    update(object kwobj,bool bIgnoreUnknown)
        object  arguments()
        object  specified()''',file=parameters_pxi)
  print('''#ifndef PKD_PARAMETERS_H
#define PKD_PARAMETERS_H 1
#include "pyrameters.h"
class pkd_parameters : public pyrameters {
public:
    pkd_parameters() = default;
    pkd_parameters(PyObject *arguments, PyObject *specified) : pyrameters(arguments,specified) {}''',file=parameters_h)

  w = 0
  for sk,sv in f.items():
    w = max(w,max(map(len, sv)))

  for sk,sv in f.items():
    emit_str(sv.items(),parameters_h,parameters_pxi)
  for sk,sv in f.items():
    emit_proto(sv.items(),parameters_h,parameters_pxi)

  print('''};
#endif
''',file=parameters_h)
