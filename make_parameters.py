from sys import argv,exit,version_info
if version_info >= (3,11):
  import tomllib as tl
else:
  import tomli as tl

if len(argv) <= 3: exit('Usage: {} toml parameters.h parameters.pxi'.format(argv[0]))

with open(argv[1],"rb") as fp:
  f=tl.load(fp)


with open(argv[2], 'w') as parameters_h,open(argv[3], 'w') as parameters_pxi:
  print('''cdef extern from "pkd_parameters.h":
    cdef cppclass pkd_parameters:
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
    for k,v in sv.items():
      print(f'    static constexpr auto str_{k:{w}} = "{k}";',file=parameters_h)
  for sk,sv in f.items():
    for k,v in sv.items():
      default = v['default']
      help = v['help']
      size = v['size'] if 'size' in v else None
      print(f'    /// {help}',file=parameters_h)
      if isinstance(default,float):
        print(f'    auto get_{k:{w}}() {{ return get<double>(str_{k}); }}',file=parameters_h)
        print(f'        float   get_{k}()',file=parameters_pxi)
      elif isinstance(default,bool):
        print(f'    auto get_{k:{w}}() {{ return get<bool>(str_{k}); }}',file=parameters_h)
        print(f'        bool    get_{k}()',file=parameters_pxi)
      elif isinstance(default,int):
        print(f'    auto get_{k:{w}}() {{ return get<int64_t>(str_{k}); }}',file=parameters_h)
        print(f'        int64_t get_{k}()',file=parameters_pxi)
      elif isinstance(default,str):
        print(f'    auto get_{k:{w}}() {{ return get<std::string>(str_{k}); }}',file=parameters_h)
        print(f'        string  get_{k}()',file=parameters_pxi)
      else:
        print(f'    auto get_{k:{w}}() {{ return get<PyObject*>(str_{k}); }}',file=parameters_h)
        print(f'        object  get_{k}()',file=parameters_pxi)
      print(f'    bool has_{k:{w}}() {{ return has(str_{k}); }}',file=parameters_h)
      print(f'        bool    has_{k}()',file=parameters_pxi)

  print('''};
#endif
''',file=parameters_h)
