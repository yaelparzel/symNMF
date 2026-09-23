# Original file from HW1 (Without main function + changing defualt contants)

import sys

EPS = 1e-4
DEFAULT_ITER = 300

def is_natural(potential_natural):
    """
    Check if the input is a natural number (positive integer)
    Input: potential_natural - the input to check
    Output: True if it is a natural number, False otherwise
    """
    # 3 and 3.0 are ok 3.5/abc/-3/3.00001 are invalid
    if not isinstance(potential_natural, str):
        return False
    potential_natural = potential_natural.strip()
    if not potential_natural:
        return False
    if potential_natural[0] == "-":
        return False
    if potential_natural[0] == "+":
        potential_natural = potential_natural[1:]
        if not potential_natural:
            return False
    if "." in potential_natural:
        if potential_natural.count(".") > 1:
            return False
        whole, frac = potential_natural.split(".")
        if not whole.isdigit():
            return False
        if not all(c == "0" for c in frac):
            return False
    else:
        if not potential_natural.isdigit():
            return False
    return True

def read_points(stream):
    """
    Read the input stream into a list of points
    Input: stream - the file or text stream containing comma-separated coordinates
    Output: a list of points, where each point is a list of floats
    """
    points = []
    for line in stream:
        line = line.strip()
        if line == "":
            continue
        coords = [float(x) for x in line.split(",")]
        points.append(coords)
    return points


def distance(p, q):
    """
    Calculate the Euclidean distance between two points
    Input: p, q - two points with the same dimension
    Output: the Euclidean distance between them
    """
    total = sum((p[i] - q[i]) ** 2 for i in range(len(p)))
    return total ** 0.5

def closest_centroid(point, centroids):
    """
    Find the centroid closest to a given point
    Input: point - a data point, centroids - a list of centroid points
    Output: the index of the nearest centroid
    """
    best_index = 0
    best_dist = distance(point, centroids[0])
    for i in range(1, len(centroids)):
        d = distance(point, centroids[i])
        if d < best_dist:
            best_dist = d
            best_index = i
    return best_index

# start with first k points as centroids and loop: update centroids until nothing changes more than eps or hit max_iter
def kmeans(points, k, max_iter, eps=EPS):
    """
    Run the k-means clustering algorithm until convergence
    Input: points - data points, k - number of clusters, max_iter - maximum iterations, eps - convergence threshold
    Output: the final centroids for the clusters
    """
    dim = len(points[0])
    centroids = [list(points[i]) for i in range(k)]

    for _ in range(max_iter):
        sums = [[0.0] * dim for _ in range(k)]
        counts = [0] * k
 
        for point in points:
            idx = closest_centroid(point, centroids)
            counts[idx] += 1
            bucket = sums[idx]
            for j in range(dim):
                bucket[j] += point[j]

        max_delta = 0.0
        new_centroids = []
        for i in range(k):
            if counts[i] > 0:
                new_center = [sums[i][j] / counts[i] for j in range(dim)]
            else:
                new_center = centroids[i]
            delta = distance(new_center, centroids[i])
            if delta > max_delta:
                max_delta = delta
            new_centroids.append(new_center)

        centroids = new_centroids

        if max_delta < eps:
            break
    
    return centroids


def format_centroids(centroids):
    """
    Format centroids as comma-separated values with 4 decimal places
    Input: centroids - a list of centroid points
    Output: a string containing the formatted centroid rows
    """
    lines = []
    for c in centroids:
        lines.append(",".join("%.4f" % value for value in c))
    return "\n".join(lines)
