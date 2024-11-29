import graphviz
import rich

"""'
graph = {
    AGCCG: 
    GCCGA:


}
"""


def create_graph():
    """
    STEPS:
    1) Take contiguos pairs of k-mers. Vi, Vj, Eij
    2) check if Vi or Vj is in the graph
        - if both: add 1 to Eij.
        - if one of them: add either Vi or Vj to G.
        - if none: add Vi, Vj, Eij to G.
    """
