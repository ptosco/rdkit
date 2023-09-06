from __future__ import annotations
import rdkit.VLib.Output
import typing
from rdkit.VLib.Node import VLibNode
import rdkit.VLib.Node

__all__ = [
    "OutputNode",
    "VLibNode"
]


class OutputNode(rdkit.VLib.Node.VLibNode):
    """
    base class for nodes which dump output

       Assumptions:

         - destination supports a write() method

         - strFunc, if provided, returns a string representation of
           the input

         - inputs (parents) can be stepped through in lockstep


       Usage Example:
       
         >>> from rdkit.VLib.Supply import SupplyNode
         >>> supplier = SupplyNode(contents=[1,2,3])
         >>> from io import StringIO
         >>> sio = StringIO()
         >>> node = OutputNode(dest=sio,strFunc=lambda x:'%s '%(str(x)))
         >>> node.AddParent(supplier)
         >>> node.next()
         1
         >>> sio.getvalue()
         '1 '
         >>> node.next()
         2
         >>> sio.getvalue()
         '1 2 '

       
    """
    pass
