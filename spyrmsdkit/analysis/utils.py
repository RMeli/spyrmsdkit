"""
Utility functions for analysis.
"""

import numpy as np


def adjacency_matrix(atoms):
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
    n_atoms = len(atoms)

    # Allocate adjacency matrix
    A = np.zeros((n_atoms, n_atoms), dtype=int)

    # Map bond indices to selection adjacency matrix
    b = atoms.bonds.to_indices()
    _, indices_flat = np.unique(b, return_inverse=True)
    indices = indices_flat.reshape(b.shape)

    A = np.zeros((n_atoms, n_atoms), dtype=int)
    A[indices[:, 0], indices[:, 1]] = 1

    return A + A.T
