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
    n_atoms = len(ag)

    # Allocate adjacency matrix
    A = np.zeros((n_atoms, n_atoms), dtype=int)

    # Map bond indices to selection adjacency matrix
    b = ag.bonds.to_indices()

    # Remove dangling bonds
    # An AtomGroup can contain bonds to atoms that are not part of the AtomGroup
    mask = np.isin(b[:, 0], ag.atoms.indices) & np.isin(
        b[:, 1], ag.atoms.indices
    )
    b = b[mask, :]

    # NumPy magic to re-index
    _, indices_flat = np.unique(b, return_inverse=True)
    indices = indices_flat.reshape(b.shape)

    A = np.zeros((n_atoms, n_atoms), dtype=int)
    A[indices[:, 0], indices[:, 1]] = 1

    return A + A.T
