from sys import argv,exit,version_info
from os.path import basename
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

def emit_proto(section,parameters_h,parameters_pxi,enumerations_pxi,enumerations_h):
  have_enum = False
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
        elif (len(default) == 1):
          item_type = get_types(default[0])[0]
          c_type = f'std::vector<{item_type}>'
          i_type = f'vector[{item_type}]'
        else:
          item_type = get_types(default[0])[0]
          c_type = f'blitz::TinyVector<{item_type},{len(default)}>'
          i_type = f'TinyVector[{item_type},BLITZ{len(default)}]'
      type_name = v['name'] if 'name' in v else k.upper()
      if 'enum' in v:
        print(f'    auto get_{k:{w}}() const {{ return {type_name}(get<{c_type}>({name})); }}',file=parameters_h)
      else:
        print(f'    auto get_{k:{w}}() const {{ return get<{c_type}>({name}); }}',file=parameters_h)
      print(f'        {i_type:<7} get_{k}()',file=parameters_pxi)
      print(f'    bool has_{k:{w}}() const {{ return has(str_{k}); }}',file=parameters_h)
      print(f'        bool    has_{k}()',file=parameters_pxi)
      if 'enum' in v:
        print(f'    auto set_{k:{w}}(const {type_name} &value) {{ set<{c_type}>({name},{c_type}(value)); return value; }}',file=parameters_h)
      else:
        print(f'    auto set_{k:{w}}(const {c_type} &value) {{ set<{c_type}>({name},value); return value; }}',file=parameters_h)
      print(f'        {i_type:<7} set_{k}({i_type} value)',file=parameters_pxi)
      if 'enum' in v:
        have_enum = True
        print('  cpdef enum class {}(int):'.format(type_name),file=enumerations_pxi)
        print('enum class {} : int {{'.format(type_name),file=enumerations_h)
        for enum_key, enum_value in v['enum'].items():
          print('    {}'.format(enum_key), file=enumerations_pxi)
          print('  {} = {},'.format(enum_key, enum_value), file=enumerations_h)
        print('};',file=enumerations_h)
    else:
      if emit_proto(v.items(),parameters_h,parameters_pxi,enumerations_pxi,enumerations_h):
        have_enum = True
  return have_enum
if len(argv) <= 5: exit('Usage: {} toml parameters.h parameters.pxi pkd_enumerations.pxi pkd_enumerations.h'.format(argv[0]))

with open(argv[1],"rb") as fp:
  f=tl.load(fp)

enumerations_h_name = argv[5]
with open(argv[2], 'w') as parameters_h,open(argv[3], 'w') as parameters_pxi,open(argv[4], 'w') as enumerations_pxi,open(enumerations_h_name, 'w') as enumerations_h:
  print('''cdef extern from "pkd_parameters.h":
    cdef cppclass pkd_parameters:
        void set(const char *name,double value)
        void set(const char *name,uint64_t value)
        void set(const char *name,int64_t value)
        bool    update(object kwobj,bool bIgnoreUnknown)
        object  arguments()
        object  specified()''',file=parameters_pxi)
  print('''#ifndef PKD_PARAMETERS_H
#define PKD_PARAMETERS_H 1
#include "pyrameters.h"''',file=parameters_h)
  print(f'#include "{basename(enumerations_h_name)}"',file=parameters_h)
  print('''class pkd_parameters : public pyrameters {
public:
    pkd_parameters() = default;
    pkd_parameters(PyObject *arguments, PyObject *specified) : pyrameters(arguments,specified) {}''',file=parameters_h)
  print('cdef extern from "{}":'.format(basename(enumerations_h_name)), file=enumerations_pxi)
  print('''#ifndef PKD_ENUMERATIONS_H
#define PKD_ENUMERATIONS_H 1
''',file=enumerations_h)
  w = 0
  for sk,sv in f.items():
    w = max(w,max(map(len, sv)))

  have_enum = False
  for sk,sv in f.items():
    emit_str(sv.items(),parameters_h,parameters_pxi)
  for sk,sv in f.items():
    if emit_proto(sv.items(),parameters_h,parameters_pxi,enumerations_pxi,enumerations_h):
      have_enum = True
  if not have_enum: # if we don't have any enums, we don't need the enum header
    print('  pass',file=enumerations_pxi)
  print('''};
#endif
''',file=parameters_h)
  print('''
#endif
''',file=enumerations_h)
