from SparseBoundaryMatrix import SparseBoundaryMatrix


## This scripts provides a simple example of usage of the code


# Create boundary matrix from filtration
matrix = SparseBoundaryMatrix()
matrix.from_simplices("data/test_file.txt")

# Reduce the matrix
matrix.gaussian_elimination()

barcodes = matrix.get_barcode(output_filename="output_barcode.txt")

print("Computed barcode:")
for (birth, death, dim) in barcodes:
  death_str = 'inf' if death == float('inf') else str(death)
  print(f"{dim} {birth} {death_str}")