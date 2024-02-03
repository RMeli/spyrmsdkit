"""
Utility functions for analysis.
"""

import numpy as np


def adjacency_matrix(ag):
    """
    Compute adjacency matrix for selection, based on bonds.

    The adjacency matrix is computed for the molecular graph consisting of
    atoms (vertices) and bonds (edges).

    Parameters
    ----------
    atoms: AtomGroup
        Selection

    Returns
    -------
    np.array
        Adjacency matrix
    """

    b = ag.get_connections("bonds", outside=False).to_indices()

    # Map bond indices to selection adjacency matrix
    # NumPy magic to re-index
    _, indices_flat = np.unique(b, return_inverse=True)
    indices = indices_flat.reshape(b.shape)

    n_atoms = len(ag)
    A = np.zeros((n_atoms, n_atoms), dtype=int)
    A[indices[:, 0], indices[:, 1]] = 1

    return A + A.T
