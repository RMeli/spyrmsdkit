import pytest
from numpy.testing import assert_allclose

import MDAnalysis as mda

from spyrmsdkit.analysis.spyrmsd import SPyRMSD

import numpy as np
import copy


class TestSPyRMSD(object):
    @pytest.fixture()
    def benzene(self):
        from collections import namedtuple

        Molecule = namedtuple("Molecule", ["coords", "bonds"])

        coords = np.array(
            [
                [+0.00000, +1.40272, 0.00000],
                [-1.21479, +0.70136, 0.00000],
                [-1.21479, -0.70136, 0.00000],
                [+0.00000, -1.40272, 0.00000],
                [+1.21479, -0.70136, 0.00000],
                [+1.21479, +0.70136, 0.00000],
            ]
        )

        bonds = np.array(
            [
                [0, 1],
                [1, 2],
                [2, 3],
                [3, 4],
                [4, 5],
                [5, 0],
            ]
        )

        return Molecule(coords, bonds)

    @pytest.fixture
    def benzene_with_hydrogens(self):
        from collections import namedtuple

        Molecule = namedtuple("Molecule", ["coords", "bonds"])

        coords = np.array(
            [
                [+0.00000, +1.40272, 0.00000],
                [0.00000, 2.49029, 0.00000],
                [-1.21479, +0.70136, 0.00000],
                [-2.15666, 1.24515, 0.00000],
                [-1.21479, -0.70136, 0.00000],
                [-2.15666, -1.24515, 0.00000],
                [+0.00000, -1.40272, 0.00000],
                [0.00000, -2.49029, 0.00000],
                [+1.21479, -0.70136, 0.00000],
                [2.15666, -1.24515, 0.00000],
                [+1.21479, +0.70136, 0.00000],
                [2.15666, 1.24515, 0.00000],
            ]
        )

        bonds = np.array(
            [
                [0, 1],
                [0, 2],
                [0, 10],
                [2, 3],
                [2, 4],
                [4, 5],
                [4, 6],
                [6, 7],
                [6, 8],
                [8, 9],
                [8, 10],
                [10, 11],
            ]
        )

        return Molecule(coords, bonds)

    @pytest.fixture()
    def benzene_traj_still(self, benzene):
        """
        Build fictitious trajectory for benzene, which remains still at every frame.
        Returns atom group of the whole universe (benzene).
        """
        B = mda.Universe.empty(
            6, atom_resindex=[0] * 6, residue_segindex=[0], trajectory=True
        )
        B.add_TopologyAttr("bonds", benzene.bonds)
        B.add_TopologyAttr("type", ["C"] * 6)
        B.add_TopologyAttr("resname", ["BNZ"])

        coords = np.zeros((7, B.atoms.n_atoms, 3))
        coords[:, :, :] = benzene.coords

        B.load_new(np.array(coords), order="fac")

        return B.atoms

    @pytest.fixture()
    def benzene_with_hydrogens_traj_still(self, benzene_with_hydrogens):
        """
        Build fictitious trajectory for benzene, which remains still at every frame.
        Returns atom group of the whole universe (benzene).
        """
        B = mda.Universe.empty(
            12, atom_resindex=[0] * 12, residue_segindex=[0], trajectory=True
        )
        B.add_TopologyAttr("bonds", benzene_with_hydrogens.bonds)
        B.add_TopologyAttr("type", ["C", "H"] * 6)
        B.add_TopologyAttr("resname", ["BNZ"])

        coords = np.zeros((7, B.atoms.n_atoms, 3))
        coords[:, :, :] = benzene_with_hydrogens.coords

        B.load_new(np.array(coords), order="fac")

        return B.atoms

    @pytest.fixture()
    def benzene_traj_rotating_random(self, benzene):
        """
        Build fictitious trajectory for benzene, rotating of 60 degrees, with additional
        random coordinates.
        Returns atom group of the whole universe (benzene and random atoms).
        """
        from MDAnalysis.lib.transformations import rotation_matrix

        angles = np.array([0, 60, 120, 180, 240, 300, 360]) * np.pi / 180
        n_frames = len(angles)

        B = mda.Universe.empty(
            20,
            n_residues=2,
            atom_resindex=[0] * 6 + [1] * 14,
            residue_segindex=[0, 0],
            trajectory=True,
        )
        B.add_TopologyAttr("bonds", benzene.bonds)
        B.add_TopologyAttr("type", ["C"] * 6 + ["X"] * 14)
        B.add_TopologyAttr("resname", ["BNZ"] + ["UNK"])

        coords = np.zeros((n_frames, 20, 3))
        for i, angle in enumerate(angles):
            cog = np.mean(benzene.coords, axis=0)
            assert np.allclose(cog, [0, 0, 0])

            R = rotation_matrix(angle, [0, 0, 1], cog)[:3, :3]
            coords[i, :6, :] = (
                benzene.coords @ R.T
            )  # Apply rotation to all coordinates
            coords[i, 6:, :] = np.random.random((14, 3))

        assert np.allclose(coords[0, :6], benzene.coords)
        assert np.allclose(coords[-1, :6], benzene.coords)

        B.load_new(np.array(coords), order="fac")

        return B.atoms

    @pytest.fixture()
    def benzene_traj_rotating(self, benzene_traj_rotating_random):
        """
        Build fictitious trajectory for benzene, rotating of 60 degrees.
        Returns atom group of benzene only.
        """
        return benzene_traj_rotating_random.select_atoms("resname BNZ")

    def test_rmsd_benzene_self(self, benzene_traj_rotating):
        assert len(benzene_traj_rotating.universe.trajectory) == 7

        R = SPyRMSD(benzene_traj_rotating)
        R.run()

        assert np.allclose(R.results.rmsd[:, -1], 0.0, atol=1e-5)

    def test_rmsd_benzene(self, benzene_traj_still, benzene_traj_rotating):
        assert len(benzene_traj_still.universe.trajectory) == 7
        assert len(benzene_traj_rotating.universe.trajectory) == 7

        R = SPyRMSD(
            benzene_traj_rotating,
            reference=benzene_traj_still,
            ref_frame=-1,
        )
        R.run()

        assert np.allclose(R.results.rmsd[:, -1], 0.0, atol=1e-5)

    def test_hydrogen_stripping(self, benzene_with_hydrogens_traj_still):
        assert len(benzene_with_hydrogens_traj_still.universe.trajectory) == 7
        assert benzene_with_hydrogens_traj_still.n_atoms == 12

        R = SPyRMSD(
            benzene_with_hydrogens_traj_still.select_atoms(
                "resname BNZ and (not type H)"
            ),
        )
        R.run()

        assert np.allclose(R.results.rmsd[:, -1], 0.0, atol=1e-5)

    def test_rmsd_benzene_rand(
        self, benzene_traj_rotating, benzene_traj_rotating_random
    ):
        assert len(benzene_traj_rotating.universe.trajectory) == 7
        assert len(benzene_traj_rotating_random.universe.trajectory) == 7

        R = SPyRMSD(
            benzene_traj_rotating_random.select_atoms("resname BNZ"),
            reference=benzene_traj_rotating,
        )
        R.run()

        assert np.allclose(R.results.rmsd[:, -1], 0.0, atol=1e-5)

    def test_isomorphisms(self, benzene_traj_still):
        """
        Test isomorphisms between two selections.
        """
        R = SPyRMSD(benzene_traj_still.atoms)
        R.run()

        assert R.isomorphisms is not None
        assert len(R.isomorphisms) == 12

        # Check identity
        assert ([0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5]) in R.isomorphisms
