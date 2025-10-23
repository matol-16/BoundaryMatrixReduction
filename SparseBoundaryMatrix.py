from Simplex import Simplex, simplex_key

class SparseBoundaryMatrix:
    
  def __init__(self, n_columns=0, n_rows=0):
    self.n_columns = n_columns
    self.n_rows=n_rows
    self.columns= dict() # conserves the ordering of insertion. each column index is associated to a set of vertices
    self.reduced= False
    self.sorted_simplices = [] # array of sorted simplices. Allows to recover filtration values
    self.sorted_simplices_indexes = dict() # allows to retrieve in O(1) the index of any sub-simplex


  def set_coef(self,i,j):
    if j not in self.columns:
        self.columns[j] = set()
    self.columns[j].add(i)

  def get_coef(self,i,j):
    return 1 if j in self.columns and i in self.columns[j] else 0

  def low(self, i: int):
    """
    Implements the low function, which yields the lowest non- zero coefficient of
    column i in matrix self.

    Storing the low values may be more efficient.
    """
    if i not in self.columns or not self.columns[i]:
        return -1
    return max(self.columns[i])

  def read_filtration_file(self, filename):
    """
    Reads the AMSCII file under the following format : f(sigma) dim(sigma) v_0 ... v_{dim(sigma)}
    output: array of simplex objects
    """
    simplices=[]
    f = open(filename,"r")
    lines = f.readlines()
    for line in lines:
      simplex = Simplex()
      line = line.split()
      filtration = float(line[0])
      simplex.add_filtration(filtration)
      dim = int(line[1])
      simplex.add_dim(dim)

      for i in range(dim+1):
        vertex = int(line[i+2])
        simplex.add_vertex(vertex)

      simplices.append(simplex)

    return simplices

  def store_sub_simplex_index(self):
    """
        Stores the index of each sub_simplex given its vertices.
        This allows to retrieve in O(1) the index of any sub-simplex
    """
    if self.sorted_simplices == []:
      print("simplices are not sorted yet. Please sort them first")
      return

    for simplex_index, simplex in enumerate(self.sorted_simplices):
      self.sorted_simplices_indexes[simplex.vertices] = simplex_index #here, we use frozen sets to have hashable keys


  def from_simplices(self, filename):
    """
        computes the boundary matrix from the vector of simplices F.

        input: sorted array of simplices F. each simplex is a Simplex Object
        output: boundary matrix B
    """
    simplices=self.read_filtration_file(filename)

    self.n_columns=len(simplices)
    self.n_rows=len(simplices)

    self.sorted_simplices = sorted(simplices, key=simplex_key)

    self.store_sub_simplex_index()

    for simplex_index, simplex in enumerate(self.sorted_simplices):
      for vertex_drop in simplex.vertices: #each face as the simplex formed by dim-1 points
        sub_simplex_vertices = simplex.vertices - {vertex_drop}
        if sub_simplex_vertices:
          sub_simplex_index= self.sorted_simplices_indexes[sub_simplex_vertices] #retrieval in O(1)
          self.set_coef(sub_simplex_index,simplex_index)

    return self

  def gaussian_elimination(self):
      """
      input: the sparse matrix obtained from a filtration.

      output: None, in place reduction
      """
      lows=dict() #store low values to optimize computations

      for j in self.columns.keys():
        equal_lows=True
        while equal_lows:
          low_j=self.low(j)
          if low_j==-1:
            break
          elif low_j in lows.keys():
            #add the column which has the same pivot to column j
            self.columns[j]^=self.columns[lows[low_j]]  #here, We use symmetric difference as we work in Z/2Z
          else:
            lows[low_j]=j
            equal_lows=False

      self.reduced=True

  def get_barcode(self, output_filename = None):
    if not self.reduced:
      print("matrix not reduced, first conducting gaussian elimination")
      self.gaussian_elimination()

    barcodes=[] #format: (k b d) for interval [b,d) in dimension k
    paired = set()

    for column_index, column_set in self.columns.items():
      pivot_row_index = self.low(column_index)
      if pivot_row_index == -1:
        continue
      if pivot_row_index < column_index: # should always be the case for a reduced matrix
        barcodes.append((self.sorted_simplices[pivot_row_index].filtration,self.sorted_simplices[column_index].filtration,self.sorted_simplices[pivot_row_index].dim))
        paired.add(column_index)
        paired.add(pivot_row_index)

    for unpaired_simplex_index in range(self.n_columns):
          if unpaired_simplex_index not in paired:
            barcodes.append((self.sorted_simplices[unpaired_simplex_index].filtration,float('inf'),self.sorted_simplices[unpaired_simplex_index].dim))

    barcodes.sort()

    if output_filename:
      with open(output_filename, 'w') as f:
          for (birth, death, dim) in barcodes:
              death_str = str(death)
              f.write(f"{dim} {birth} {death_str}\n")
      f.close()

    return barcodes