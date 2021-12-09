from sys import argv,exit
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

if len(argv) <= 2: exit('Usage: {} input output'.format(ARGV[0]))

with open(argv[1], 'rb') as parse_py:
  with open(argv[2], 'w') as m_parse:
    m_parse.write('char parse_py[] = {\n')
    hexdump(m_parse,parse_py.read())
    m_parse.write('0x00 };\n')
