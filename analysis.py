import sys
import numpy as np
import symnmfmodule
import kmeans  # Importing HW1 kmeans.py
from sklearn.metrics import silhouette_score, adjusted_rand_score

MAX_ITER = 300
EPS = 1e-4

np.random.seed(1234)

def parse_args():
    """
    Extract and validate command line arguments
    Output: k (int), file_name (str) or None, None if invalid
    """
    if len(sys.argv) != 3:
        return None, None
    try:
        k = int(sys.argv[1])
        file_name = sys.argv[2]
        return k, file_name
    except ValueError:
        return None, None

def load_data(file_name):
    """
    Read the text file into a Python list of lists
    Input: file_name (str)
    Output: X (list of lists) or None if error
    """
    try:
        X = np.loadtxt(file_name, delimiter=',', ndmin=2)
        return X.tolist()
    except Exception:
        return None

def get_symnmf_labels(X, k):
    """
    Run SymNMF and return cluster assignments
    Input: X (list of lists), k (int)
    Output: labels (list of int)
    """
    W = symnmfmodule.norm(X)
    n = len(W)
    m = sum(sum(row) for row in W) / (n * n)
    high = 2 * np.sqrt(m / k)
    
    H_init = np.random.uniform(0, high, (n, k)).tolist()
    H_final = symnmfmodule.symnmf(H_init, W)
    
    return np.argmax(np.array(H_final), axis=1)

def get_kmeans_labels(X, k):
    """
    Run K-means and return cluster assignments
    Input: X (list of lists), k (int)
    Output: labels (list of int)
    """
    centroids = kmeans.kmeans(X, k, max_iter=MAX_ITER, eps=EPS)
    return [kmeans.closest_centroid(point, centroids) for point in X]

def print_metrics(X, symnmf_labels, kmeans_labels):
    """
    Calculate and print scores
    Input: X (list of lists), symnmf_labels (list of int), kmeans_labels (list of int)
    """
    X_np = np.array(X)
    nmf_score = silhouette_score(X_np, symnmf_labels)
    kmeans_score = silhouette_score(X_np, kmeans_labels)
    ari_score = adjusted_rand_score(kmeans_labels, symnmf_labels)
    
    print(f"nmf: {nmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")
    print(f"ari: {ari_score:.4f}")

def main():
    k, file_name = parse_args()
    if k is None:
        print("An Error Has Occurred")
        sys.exit(1)
        
    X = load_data(file_name)
    if X is None:
        print("An Error Has Occurred")
        sys.exit(1)

    symnmf_labels = get_symnmf_labels(X, k)
    kmeans_labels = get_kmeans_labels(X, k)
    
    print_metrics(X, symnmf_labels, kmeans_labels)

if __name__ == "__main__":
    main()