import numpy as np


class Kmeans:
    """Kmeans clustering.

    Example usage:
        > clf = Kmeans()
        > clf.fit(patient_ID, data)
        > clf.predict(x)
    """
    def __init__(self, k=5, max_iter=10000, eps=1e-5):
        """
        Args:
            step_size: Step size for iterative solvers only.
            max_iter: Maximum number of iterations for the solver.
            eps: Threshold for determining convergence.
            theta_0: Initial guess for theta. If None, use the zero vector.
            verbose: Print loss values during training.
        """
        self.k = k
        self.tol = eps
        self.max_iter = max_iter
        self.centroids = {}
        self.clusters = {}
        self.patient_clusters = {}


    def fit(self, patient_ID, data):
        """
        Args:
            x: Training example inputs. Shape (n_examples, dim).
        """
        # *** START CODE HERE ***
        for i in range(self.k):
            self.centroids[i] = data[i]

        self.clusters = {}
        for i in range(self.max_iter):
            for j in range(self.k):
                self.clusters[j] = []
                self.patient_clusters[j] = []

            for k in range(data.shape[0]):
                x = data[k]
                distances = [np.linalg.norm(x - self.centroids[centroid]) for centroid in self.centroids]
                closest_class = distances.index(min(distances))
                self.clusters[closest_class].append(x)
                self.patient_clusters[closest_class].append(patient_ID[k])

            prev_centroids = self.centroids

            optimized = True

            for cluster in self.clusters:
                self.centroids[cluster] = np.average(self.clusters[cluster], axis=0)
                if np.sum((self.centroids[cluster] - prev_centroids[cluster]) / self.centroids[cluster] * 100.0) > self.tol:
                    optimized = False

            if optimized:
                break
        # *** END CODE HERE ***

    def predict(self, x):
        """
        Args:
            x: Inputs of shape (n_examples, dim).

        Returns:
            Outputs of shape (n_examples,).
        """
        # *** START CODE HERE ***
        # *** END CODE HERE

def main():



if __name__ == '__main__':
    main()