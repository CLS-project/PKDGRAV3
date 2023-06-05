from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.statemachine import ViewList
from sphinx.util.nodes import nested_parse_with_titles
from sys import version_info
if version_info >= (3,11):
    import tomllib as tl
else:
    import tomli as tl
import re

class make_parameters(Directive):
    required_arguments = 1
    final_argument_whitespace = False

    def emit_section(self,sk,sv):
      sect = nodes.section(ids=[sk])
      sect += nodes.title(text=sk)
      if 'docs' in sv:
        node = nodes.section()
        vl = ViewList(sv['docs'].split('\n'),source='')
        nested_parse_with_titles(self.state, vl, node)
        # items += node
        sect += node.children
      items=[]
      subs=[]
      for key,v in sv.items():
        if isinstance(v,(str,int,float,bool)): continue
        if 'default' in v:
          if 'private' in v and v['private']:
            continue
          # The default value with empty strings as none
          if 'omitted' in v:
            default=v['omitted']
          else:
            default = v['default']
            if isinstance(default,str):
              default = f'"{default}"' if len(default)>0 else "none"
          # Prefer docs over help and reparse as restructed text.
          text = v['docs'] if 'docs' in v else v['help']
          node = nodes.section()
          text = f'.. index:: single: parameters ; {key}\n\n{text}'
          vl = ViewList(text.split('\n'),source='')
          nested_parse_with_titles(self.state, vl, node)
          term=nodes.section()
          nested_parse_with_titles(self.state, nodes.paragraph(text=f'{key} (default {default})'), term)
          items += nodes.definition_list_item('',
                  nodes.term('','',*term.children),
                  nodes.definition('',*node.children))
          # dli.insert(0,nodes.target('', '', ids=[f'{key}'.lower()]))
        else:
          subs += [self.emit_section(key,v)]
      if (len(items)>0):
        sect += nodes.definition_list('',*items)
      sect += subs
      return sect

    def run(self):
        # We are given a TOML file with all of the parameters.
        # If the file changes (or this script) then we need to rebuild,
        # so we mark the dependencies.
        file = self.arguments[0]
        self.state.document.settings.env.note_dependency(file)
        self.state.document.settings.env.note_dependency(__file__)
        # Read the TOML file
        with open(file,"rb") as fp:
          f=tl.load(fp)

        docs=[]
        # Loop over each document section
        for sk,sv in f.items():
          docs.append(self.emit_section(sk,sv))
        return docs

def setup(app):
    app.add_directive('make_parameters', make_parameters)
