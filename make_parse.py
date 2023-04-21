from sys import argv,exit,version_info
if version_info >= (3,11):
  import tomllib as tl
else:
  import tomli as tl

def group(a, *ns):
  for n in ns:
    a = [a[i:i+n] for i in range(0, len(a), n)]
  return a

def join(a, *cs):
  return [cs[0].join(join(t, *cs[1:])) for t in a] if cs else a

def hexdump(file,data):
  toHex = lambda c: '0x{:02X}'.format(c)
  make = lambda f, *cs: join(group(list(map(f, data)), 16), *cs)
  hs = make(toHex, ',')
  for i, (h) in enumerate(hs):
    file.write('{},\n'.format(h))

if len(argv) <= 3: exit('Usage: {} toml parse.py m_parse.h'.format(ARGV[0]))

with open(argv[1],"rb") as fp:
  f=tl.load(fp)

with open(argv[2], 'w') as parse_py:
  print('''def parse():
    from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter,Action,Namespace
    from sys import argv

    class PkdParser(ArgumentParser):
        def __init__(self, **kwargs):
            super(PkdParser, self).__init__(**kwargs,
                formatter_class=lambda prog: ArgumentDefaultsHelpFormatter(prog,max_help_position=50, width=132))
            self._specified = Namespace()
        def setSpecified(self,dest):
            setattr(self._specified,dest,True)
        def getSpecified(self):
            return self._specified

    def add_bool(self,name,default=None,**kwargs):
        class BoolAction(Action):
            def __init__(self, option_strings, dest, nargs=None, **kwargs):
                if nargs is not None: raise ValueError("nargs not allowed")
                super(BoolAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
            def __call__(self, parser, namespace, values, option_string=None):
                parser.setSpecified(self.dest)
                setattr(namespace, self.dest, True if option_string[0:1]=='+' else False)
        if default:
            return self.add_argument('-'+name,'+'+name,action=BoolAction,default=default,**kwargs)
        else:
            return self.add_argument('+'+name,'-'+name,action=BoolAction,default=default,**kwargs)

    def add_flag(self,name,default=None,**kwargs):
        class FlagAction(Action):
            def __init__(self, option_strings, dest, **kwargs):
                super(FlagAction, self).__init__(option_strings, dest, **kwargs)
            def __call__(self, parser, namespace, values, option_string=None):
                parser.setSpecified(self.dest)
                setattr(namespace, self.dest, values)
        return self.add_argument('-'+name,default=default,action=FlagAction,**kwargs)

    parser = PkdParser(description='PKDGRAV3 n-body code',prefix_chars='-+',add_help=False)
    parser.add_argument('--help',action='help',help='show this help message and exit')
''',file=parse_py)

  section=0
  for sk,sv in f.items():
    section+=1
    print(f"    s{section} = parser.add_argument_group('{sk}')",file=parse_py)
    for k,v in sv.items():
      default = v['default']
      flag = v['flag']
      help = v['help']
      size = v['size'] if 'size' in v else None
      if isinstance(default,float):
        print(f"    add_flag(s{section},'{flag}',default={default:.15g},dest='{k}',type=float,help='{help}')",file=parse_py)
      elif isinstance(default,bool):
        print(f"    add_bool(s{section},'{flag}',default={default},dest='{k}',help='{help}')",file=parse_py)
      elif isinstance(default,int):
        print(f"    add_flag(s{section},'{flag}',default={default},dest='{k}',type=int,help='{help}')",file=parse_py)
      elif isinstance(default,str):
        if len(default)>0:
          print(f"    add_flag(s{section},'{flag}',default='{default}',dest='{k}',help='{help}')",file=parse_py)
        else:
          print(f"    add_flag(s{section},'{flag}', dest='{k}',help='{help}')",file=parse_py)
      elif isinstance(default,list):
        print(f"    s{section}.add_argument('-{flag}',dest='{k}',type=int,action='append',help='{help}')",file=parse_py)
      else:
        print(type(default),file=parse_py)
    print('',file=parse_py)

  print('''
    parser.add_argument('script',nargs='?',default=None,help='File containing parameters or analysis script')

    (params,extra) = parser.parse_known_args()
    spec = parser.getSpecified()
    for k in vars(params):
        if not k in vars(spec): setattr(spec,k,False)
    argv[1:] = extra # Consume the parameters we parsed out
    if params.script is not None: argv[0]=params.script
    return (params,spec)

def update(pars,args,spec):
    for key,value in pars.items():
        if key in vars(args) and not getattr(spec,key):
            setattr(args,key,value)
            setattr(spec,key,True)
        # else: this is a rogue variable?''',file=parse_py)

with open(argv[2], 'rb') as parse_py:
  with open(argv[3], 'w') as m_parse:
    m_parse.write('char parse_py[] = {\n')
    hexdump(m_parse,parse_py.read())
    m_parse.write('0x00 };\n')
