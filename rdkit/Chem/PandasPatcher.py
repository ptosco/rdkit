import logging
import re
import importlib
from xml.dom import minidom
from xml.parsers.expat import ExpatError
from rdkit.Chem import Mol

log = logging.getLogger(__name__)
RDK_MOLS_AS_IMAGE_ATTR = "__rdkitMolAsImage"
InteractiveRenderer = None
PrintAsBase64PNGString = None
molJustify = None
pandas_frame = None

pandas_formats = None
for pandas_formats_name in ("pandas.io", "pandas"):
  try:
    pandas_formats = importlib.import_module(f"{pandas_formats_name}.formats")
  except ModuleNotFoundError:
    continue
  break
if pandas_formats is None:
  log.warning("Failed to import pandas formats module")
  raise ModuleNotFoundError
to_html_class = None
for to_html_class_name in ("DataFrameRenderer", "DataFrameFormatter"):
  if hasattr(pandas_formats, "format") and hasattr(pandas_formats.format, to_html_class_name):
    to_html_class = getattr(pandas_formats.format, to_html_class_name)
    if hasattr(to_html_class, "to_html"):
      break
    else:
      to_html_class = None
if to_html_class is None:
  log.warning("Failed to find the pandas to_html method to patch")
  raise AttributeError
orig_get_adjustment = None
for get_adjustment_name in ("get_adjustment", "_get_adjustment"):
  if hasattr(pandas_formats.format, get_adjustment_name):
    orig_get_adjustment = getattr(pandas_formats.format, get_adjustment_name)
    break
if orig_get_adjustment is None:
  log.warning("Failed to find the pandas get_adjustment() function to patch")
  raise AttributeError
adj = orig_get_adjustment()
if not hasattr(adj, "len"):
  log.warning(f"Failed to find the pandas {adj.__class.__name__}.len() method to patch")
  raise AttributeError
write_cell_class = None
for html_formatter_module_name in ("format", "html"):
  try:
    write_cell_class = importlib.import_module(f"{pandas_formats.__name__}.{html_formatter_module_name}")
  except ModuleNotFoundError:
    continue
  if hasattr(write_cell_class, "HTMLFormatter"):
    write_cell_class = getattr(write_cell_class, "HTMLFormatter")
    break
  else:
    write_cell_class = None
if write_cell_class is None:
  log.warning("Failed to find the pandas HTMLFormatter class to patch")
  raise AttributeError
orig_write_cell = None
if not hasattr(write_cell_class , "_write_cell"):
  log.warning("Failed to find the HTMLFormatter._write_cell() method to patch")
  raise AttributeError
orig_write_cell = getattr(write_cell_class , "_write_cell")
if not (hasattr(pandas_formats, "printing") and hasattr(pandas_formats.printing, "pprint_thing")):
  log.warning("Failed to find the pprint_thing function")
  raise AttributeError
try:
  import pandas as pd
except ImportError:
  log.warning("Failed to import pandas")
  raise


orig_to_html = getattr(to_html_class, "to_html")
pprint_thing = pandas_formats.printing.pprint_thing

def is_molecule_image(s):
  result = False
  try:
    # is text valid XML / HTML?
    xml = minidom.parseString(s)
    root_node = xml.firstChild
    # check data-content attribute
    if (root_node.nodeName in ['svg', 'img', 'div'] and
        'data-content' in root_node.attributes.keys() and
        root_node.attributes['data-content'].value == 'rdkit/molecule'):
      result = True
  except ExpatError:
    pass  # parsing xml failed and text is not a molecule image
  return result

styleRegex = re.compile("^(.*style=[\"'][^\"^']*)([\"'].*)$")

def is_mol(x):
  return isinstance(x, Mol)

def default_formatter(x):
  return pprint_thing(x, escape_chars=("\t", "\r", "\n"))

def mol_formatter(x):
  return is_mol(x) and PrintAsBase64PNGString(x) or default_formatter(x)

def get_mol_formatters(df):
  return {col: mol_formatter for col in df.select_dtypes("object").applymap(is_mol).any().keys()}

def check_rdk_attr(frame, attr):
  return hasattr(frame, attr) and getattr(frame, attr)

def set_rdk_attr(frame, attr):
  setattr(frame, attr, True)

def patched_to_html(self, *args, **kwargs):
  """A patched version of the to_html method
     that allows rendering molecule images in data frames.
  """
  frame = None
  if self.__class__.__name__ == "DataFrameRenderer":
    fmt = self.fmt
    frame = fmt.frame
  elif self.__class__.__name__ == "DataFrameFormatter":
    fmt = self
    frame = fmt.frame
  else:
    raise ValueError(f"patched_to_html: unexpected class {self.__class__.__name__}")
  if not check_rdk_attr(frame, RDK_MOLS_AS_IMAGE_ATTR):
    return orig_to_html(self, *args, **kwargs)
  orig_formatters = fmt.formatters
  try:
    formatters = orig_formatters or {}
    if not isinstance(formatters, dict):
      formatters = {col: formatters[i] for i, col in enumerate(self.columns)}
    else:
      formatters = dict(formatters)
    formatters.update(get_mol_formatters(frame))
    fmt.formatters = formatters
    res = orig_to_html(self, *args, **kwargs)
    should_inject = InteractiveRenderer and InteractiveRenderer.isEnabled()
    return InteractiveRenderer.injectHTMLHeaderBeforeTable(res) if should_inject else res
  finally:
    fmt.formatters = orig_formatters

def patched_write_cell(self, s, *args, **kwargs):
  """ Disable escaping of HTML in order to render img / svg tags """
  styleTags = f"text-align: {molJustify};"
  def_escape = self.escape
  try:
    if hasattr(self.frame, RDK_MOLS_AS_IMAGE_ATTR) and is_molecule_image(s):
      self.escape = False
      kind = kwargs.get('kind', None)
      if kind == 'td':
        tags = kwargs.get('tags', None) or ''
        match = styleRegex.match(tags)
        if match:
          tags = styleRegex.sub(f'\\1 {styleTags}\\2', tags)
        else:
          if tags:
            tags += ' '
          tags += f'style="{styleTags}"'
        kwargs['tags'] = tags
    return orig_write_cell(self, s, *args, **kwargs)
  finally:
    self.escape = def_escape

def patched_get_adjustment():
  """ Avoid truncation of data frame values that contain HTML content """
  adj = orig_get_adjustment()
  orig_len = adj.len
  adj.len = lambda text, *args, **kwargs: (
    0 if is_molecule_image(text) else orig_len(text, *args, **kwargs)
  )
  return adj

def renderImagesInAllDataFrames(images=True):
  if images:
    set_rdk_attr(pd.core.frame.DataFrame, RDK_MOLS_AS_IMAGE_ATTR)
  elif hasattr(pd.core.frame.DataFrame, RDK_MOLS_AS_IMAGE_ATTR):
    delattr(pd.core.frame.DataFrame, RDK_MOLS_AS_IMAGE_ATTR)

def changeMoleculeRendering(frame, renderer='image'):
  if not renderer.lower().startswith('str'):
    set_rdk_attr(frame, RDK_MOLS_AS_IMAGE_ATTR)
  elif hasattr(frame, RDK_MOLS_AS_IMAGE_ATTR):
    delattr(frame, RDK_MOLS_AS_IMAGE_ATTR)

def patchPandas():
  if getattr(to_html_class, "to_html") != patched_to_html:
    setattr(to_html_class, "to_html", patched_to_html)
  if getattr(write_cell_class, "_write_cell") != patched_write_cell:
    setattr(write_cell_class, "_write_cell", patched_write_cell)
  if getattr(pandas_formats.format, get_adjustment_name) != patched_get_adjustment:
    setattr(pandas_formats.format, get_adjustment_name, patched_get_adjustment)

def unpatchPandas():
  if orig_to_html:
    setattr(to_html_class, "to_html", orig_to_html)
  if orig_write_cell:
    setattr(write_cell_class, "_write_cell", orig_write_cell)
  if orig_get_adjustment:
    setattr(pandas_formats.format, get_adjustment_name, orig_get_adjustment)
