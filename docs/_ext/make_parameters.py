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

# Global dictionary to store data from files
parameters_toml = None
parameters_dict = {}

class ParametersDirective(Directive):
    def format_parameter(self,key,v,main_entry=True):
      # The default value with empty strings as none
      if 'omitted' in v:
        default=v['omitted']
      else:
        default = v['default']
        if isinstance(default,str):
          default = f'"{default}"' if len(default)>0 else "none"
        if isinstance(default,list) and len(default) <= 1:
          default = "none"
      if isinstance(v['default'],float):
        classifier = 'float'
      elif isinstance(v['default'],bool):
        classifier = 'Boolean'
      elif isinstance(v['default'],int):
        classifier = 'integer'
      elif isinstance(v['default'],str):
        classifier = 'string'
      else:
        classifier = None

      # Prefer docs over help and reparse as restructed text.
      text = v['docs'] if 'docs' in v else v['help']
      node = nodes.section()
      vl = ViewList(text.split('\n'),source='')
      nested_parse_with_titles(self.state, vl, node)
      text = f'.. index:: {"! " if main_entry else ""}single: parameters ; {key}\n'\
             f'{f"  :name: {key.lower()}" if main_entry else ""}\n'\
             f'\n{key} (default {default})'
      term=nodes.section()
      vl = ViewList(text.split('\n'),source='')
      nested_parse_with_titles(self.state, vl, term)
      dli = nodes.definition_list_item()
      dli += nodes.term('','',*term.children[0:-1],*term.children[-1].children)
      if classifier is not None:
        dli += nodes.classifier('default',classifier)
      dli += nodes.definition('',*node.children)
      return dli

class load_parameters(Directive):
    # This directive takes one argument: the filename
    required_arguments = 1

    def record_section(self,sk,sv):
      for key,v in sv.items():
        if isinstance(v,(str,int,float,bool)): continue
        if 'default' in v:
           parameters_dict[key] = v
        else:
          self.record_section(key,v)

    def run(self):
        # We are given a TOML file with all of the parameters.
        # If the file changes (or this script) then we need to rebuild,
        # so we mark the dependencies.
        file = self.arguments[0]
        self.state.document.settings.env.note_dependency(file)
        self.state.document.settings.env.note_dependency(__file__)
        # Read the TOML file
        with open(file,"rb") as fp:
          parameters_toml = tl.load(fp)

        # Loop over each document section
        parameters_dict = {}
        for sk,sv in parameters_toml.items():
          self.record_section(sk,sv)

        # This directive does not modify the document
        return []

class show_parameters(ParametersDirective):
    # This directive takes at least one argument: the parameter
    required_arguments = 1
    optional_arguments = 1000

    def run(self):
      items = []
      for key in self.arguments:
        v = parameters_dict[key]
        items += self.format_parameter(key,v,False)
      return [nodes.definition_list('',*items)]

class make_parameters(ParametersDirective):
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
          items += self.format_parameter(key,v)
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
    app.add_directive('load_parameters', load_parameters)
    app.add_directive('make_parameters', make_parameters)
    app.add_directive('show_parameters', show_parameters)
