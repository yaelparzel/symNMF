import sys
import numpy as np
import symnmfmodule

np.random.seed(1234)

def print_matrix(mat):
    """
    Print matrix with 4 decimal places, comma separated
    Input: mat - a list of lists representing the matrix
    """
    for row in mat:
        print(",".join([f"{val:.4f}" for val in row]))

def parse_args():
    """ 
    Extract and validate command line arguments
    Output: k, goal, file_name if valid, otherwise None, None, None
    """
    if len(sys.argv) != 4:
        return None, None, None
    try:
        k = int(sys.argv[1])
    except ValueError:
        return None, None, None
    
    goal = sys.argv[2]
    file_name = sys.argv[3]
    return k, goal, file_name

def load_data(file_name):
    """
    Read the text file into a Python list of lists
    Input: file_name - path to the text file
    Output: a list of lists representing the matrix or None if an error occurred
    """
    try:
        X = np.loadtxt(file_name, delimiter=',', ndmin=2)
        return X.tolist()
    except Exception:
        return None

def execute_goal(k, goal, X):
    """
    Route to the appropriate C function based on the requested goal
    Input: k - number of clusters, goal - the goal to achieve, X - the input matrix
    Output: the result of the computation or None if an error occurred
    """
    if goal == 'sym':
        return symnmfmodule.sym(X)
    elif goal == 'ddg':
        return symnmfmodule.ddg(X)
    elif goal == 'norm':
        return symnmfmodule.norm(X)
    elif goal == 'symnmf':
        W = symnmfmodule.norm(X)
        n = len(W)
        m = sum(sum(row) for row in W) / (n * n)
        high = 2 * np.sqrt(m / k)
        
        H_init = np.random.uniform(0, high, (n, k)).tolist()
        return symnmfmodule.symnmf(H_init, W)
    return None

def main():
    k, goal, file_name = parse_args()
    if k is None:
        print("An Error Has Occurred")
        sys.exit(1)
        
    X = load_data(file_name)
    if X is None or len(X) == 0:
        print("An Error Has Occurred")
        sys.exit(1)
    
    if not 1 <= k < len(X):
        print("An Error Has Occurred")
        sys.exit(1)

    try:
        result = execute_goal(k, goal, X)
    except Exception:
        result = None

    if result is None:
        print("An Error Has Occurred")
        sys.exit(1)
        
    print_matrix(result)

if __name__ == "__main__":
    main()