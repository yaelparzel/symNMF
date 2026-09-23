# Original file from HW1

import sys

EPS = 0.001
DEFAULT_ITER = 400

def is_natural(potential_natural):
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
    points = []
    for line in stream:
        line = line.strip()
        if line == "":
            continue
        coords = [float(x) for x in line.split(",")]
        points.append(coords)
    return points


def distance(p, q):
    total = sum((p[i] - q[i]) ** 2 for i in range(len(p)))
    return total ** 0.5

def closest_centroid(point, centroids):
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
    lines = []
    for c in centroids:
        lines.append(",".join("%.4f" % value for value in c))
    return "\n".join(lines)


def main():
    args = sys.argv[1:]

    if len(args) == 1:
        k_arg, iter_arg = args[0], None
    elif len(args) == 2:
        k_arg, iter_arg = args[0], args[1]
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    if not is_natural(k_arg):
        print("Incorrect number of clusters!")
        sys.exit(1)
    k = int(float(k_arg))

    if iter_arg is None:
        max_iter = DEFAULT_ITER
    else:
        if not is_natural(iter_arg):
            print("Incorrect maximum iteration!")
            sys.exit(1)
        max_iter = int(float(iter_arg))
        if not (1 < max_iter < 800):
            print("Incorrect maximum iteration!")
            sys.exit(1)

    try:
        points = read_points(sys.stdin)
    except Exception:
        print("An Error Has Occurred")
        sys.exit(1)

    number_of_points = len(points)

    if number_of_points == 0:
        print("An Error Has Occurred")
        sys.exit(1)

    if not (1 < k < number_of_points):
        print("Incorrect number of clusters!")
        sys.exit(1)

    try:
        centroids = kmeans(points, k, max_iter)
    except Exception:
        print("An Error Has Occurred")
        sys.exit(1)

    print(format_centroids(centroids))


if __name__ == "__main__":
    main()