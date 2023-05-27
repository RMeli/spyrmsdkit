import MDAnalysis as mda
from spyrmsdkit.analysis import utils
import numpy as np

def test_adjacency_matrix():
    """
    Test construction of the adjacency matrix.
    """
    bonds = np.array([
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 5],
        [5, 0],
    ])

    benzene = mda.Universe.empty(
        6, atom_resindex=[0] * 6, residue_segindex=[0], trajectory=True
    )
    benzene.add_TopologyAttr('bonds', bonds)

    A = utils.adjacency_matrix(benzene.atoms)

    expected = np.array([
        [0, 1, 0, 0, 0, 1],
        [1, 0, 1, 0, 0, 0],
        [0, 1, 0, 1, 0, 0],
        [0, 0, 1, 0, 1, 0],
        [0, 0, 0, 1, 0, 1],
        [1, 0, 0, 0, 1, 0],
    ])

    np.testing.assert_array_equal(A, expected)

