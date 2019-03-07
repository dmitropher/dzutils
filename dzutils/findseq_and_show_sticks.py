import pymol
def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                             'PDB Loader Service',
                             label = 'PDB Loader Service',
                             command = lambda s=self : fetchPDBDialog(s))

def findseq_and_show_sticks (needle, haystack="", selname="" , het="", firstOnly="", radius=5):
    """
    Finds the sequence given and shows sticks of the selection and surrounding radius (default 5 A)
    
    Colors selection orange and the rest green.
    """
    objects = pymol.cmd.get_names("objects")
    for name in objects:
        selname = f"""{name}_seq"""
        pymol.cmd.findseq (needle, name, selName=selName, het=het, firstOnly=firstOnly)
        pymol.cmd.show ("sticks"  , f"""{name} & ({selname} + ({selname} around {radius})) & ! (e. h + name c+n+o)""")
        pymol.util.cbag(name)
        pymol.util.cbao(selname)

pymol.cmd.extend('findseq_and_show_sticks', findseq_and_show_sticks)
