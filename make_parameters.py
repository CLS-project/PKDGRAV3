from sys import argv,exit,version_info
if version_info >= (3,11):
  import tomllib as tl
else:
  import tomli as tl

if len(argv) <= 2: exit('Usage: {} toml parameters.h'.format(ARGV[0]))

with open(argv[1],"rb") as fp:
  f=tl.load(fp)


with open(argv[2], 'w') as parameters_h:
  print('''#include "pyrameters.h"
class pkd_parameters : public pyrameters {
public:
    pkd_parameters() = default;
    pkd_parameters(PyObject *arguments, PyObject *specified) : pyrameters(arguments,specified) {}''',file=parameters_h)

  w = 0
  for sk,sv in f.items():
    w = max(w,max(map(len, sv)))
  for sk,sv in f.items():
    for k,v in sv.items():
      default = v['default']
      help = v['help']
      size = v['size'] if 'size' in v else None
      print(f'    /// {help}',file=parameters_h)
      if isinstance(default,float):
        print(f'    auto get_{k:{w}}() {{ return get<double>("{k}"); }}',file=parameters_h)
      elif isinstance(default,bool):
        print(f'    auto get_{k:{w}}() {{ return get<bool>("{k}"); }}',file=parameters_h)
      elif isinstance(default,int):
        print(f'    auto get_{k:{w}}() {{ return get<int64_t>("{k}"); }}',file=parameters_h)
      elif isinstance(default,str):
        print(f'    auto get_{k:{w}}() {{ return get<std::string>("{k}"); }}',file=parameters_h)
      else:
        print(f'    auto get_{k:{w}}() {{ return get<PyObject*>("{k}"); }}',file=parameters_h)
      print(f'    bool has_{k:{w}}() {{ return has("{k}"); }}',file=parameters_h)

      # else:
      #   print(type(default))
  # print('    struct {',file=parameters_h)
  # for sk,sv in f.items():
  #   for k,v in sv.items():
  #       print(f'        bool {k:{w}}() {{ return has("{k}"); }}',file=parameters_h)

  # print('    } specified;',file=parameters_h)

  print('};',file=parameters_h)
