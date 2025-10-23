class Simplex:
    def __init__(self):
        self.vertices = frozenset() #set of integers refering to the vertexes of the simplex
        self.filtration= 0.0 #float
        self.dim = 0 #integer

    def add_vertex(self,vertex):
        self.vertices = self.vertices | {vertex}

    def add_filtration(self,filtration):
        self.filtration = filtration

    def add_dim(self,dim):
        self.dim = dim
    

def simplex_key(simplex):
    """
    allows to first sort by filtration value, then by dimension, and finally by lexicographic order.
    (tuples are sorted by lexicographic order by python)
"""
    return (simplex.filtration, simplex.dim, tuple(sorted(simplex.vertices)))
