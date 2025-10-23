import time
import matplotlib.pyplot as plt
from SparseBoundaryMatrix import SparseBoundaryMatrix

### calculate computation times of our methods for the given filtrations
### Compare them with their theoretical complexities


def evaluate_time_performance(dataset_letters = ["B","C","A","D"]):
  """
  Evaluates the time performance of the algorithms in SparseBoundaryMatrix on the datasets under the format filtration_X.txt, for X a letter.
  output: a dictionnary {number of simplices: computation time}t
  """
  #for each filtration_X.txt file for X from A to D, run our SparseBoundaryMatrix pipeline
  #record computation time and the number m of simplices
  performance= dict() #keys are the number of simplices, values the corresponding computation time

  for x in dataset_letters:
    print(f"------------start computation for {x}--------------")
    start_time = time.time()
    matrix = SparseBoundaryMatrix()
    print("get boundary from filtration")
    matrix.from_simplices(f"data/filtration_{x}.txt")
    print("--time ="+str(time.time()-start_time))
    print("conduct gaussian elimination")
    matrix.gaussian_elimination()
    print("--time ="+str(time.time()-start_time))
    print("calculate barcode")
    barcode = matrix.get_barcode()
    print("--time ="+str(time.time()-start_time))
    end_time = time.time()
    performance[len(matrix.sorted_simplices)] = end_time - start_time
    print(f"end computation for {x}, time: {end_time-start_time}; m={len(matrix.sorted_simplices)}")
  return performance

if __name__ == "__main__":
    performance = evaluate_time_performance() # This may take up to 10 minutes
    #plot the performance dictionnary
    plt.plot(list(performance.keys()), list(performance.values()))
    plt.xlabel('Number of simplices')
    plt.ylabel('Computation time (s)')
    plt.savefig("plots/performance_measure.png")
    plt.show()
    
