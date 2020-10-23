import numpy as np


class Kmeans:
    """Kmeans clustering.

    Example usage:
        > clf = Kmeans()
        > clf.fit(patient_ID, data)
        > clf.predict(x)
    """
    def __init__(self, k=5, predictor_count=30, max_iter=10000, eps=1e-5):
        """
        Args:
            step_size: Step size for iterative solvers only.
            max_iter: Maximum number of iterations for the solver.
            eps: Threshold for determining convergence.
            theta_0: Initial guess for theta. If None, use the zero vector.
            verbose: Print loss values during training.
        """
        self.k = k
        self.predictor_count = predictor_count
        self.tol = eps
        self.max_iter = max_iter
        self.centroids = {}
        self.clusters = {}
        self.patient_clusters = {}

    def read_data(self, csv_path):
        # Load features
        headers = np.loadtxt(csv_path, delimiter=',', dtype="S", max_rows=1)

        x_cols = np.asarray((np.arange(1, len(headers))))
        patient_cols = [0]

        inputs = np.loadtxt(csv_path, delimiter=',', skiprows=1, usecols=x_cols)
        patient_ID = np.loadtxt(csv_path, dtype="S", delimiter=',', skiprows=1, usecols=patient_cols)
        feature_names = headers[1:]

        frequency = np.sum(inputs, axis=0)
        ind = np.argpartition(frequency, -1 * self.predictor_count)[-1 * self.predictor_count:]
        ind = ind[np.argsort(frequency[ind])]

        final_features = inputs[:, ind]
        feature_names = feature_names[ind]

        return patient_ID, final_features, feature_names


    def fit(self, patient_ID, data):
        """
        Args:
            x: Training example inputs. Shape (n_examples, dim).
        """
        # *** START CODE HERE ***
        for i in range(self.k):
            self.centroids[i] = data[i]

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
    model = Kmeans()

    patient_ID, features, feature_names = model.read_data("/Users/juliapark/GitHub/julia-evolution-data/matrix-chrY.csv")
    model.fit(patient_ID, features)
    print(model.patient_clusters[0])
    print(model.patient_clusters[1])
    print(model.patient_clusters[2])
    print(model.patient_clusters[3])
    print(model.patient_clusters[4])





if __name__ == '__main__':
    main()