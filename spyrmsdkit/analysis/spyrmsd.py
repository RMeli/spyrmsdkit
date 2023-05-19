"""
SPyRMSD --- :mod:`spyrmsdkit.analysis.SPyRMSD`
===========================================================

This module contains the :class:`SPyRMSD` class.

"""
from MDAnalysis.analysis.base import AnalysisBase

from typing import Union, TYPE_CHECKING
import logging
import numpy as np
# TODO: Add due

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe, AtomGroup

logger = logging.getLogger("MDanalysis.analysis.spyrmsd")

class SPyRMSD(AnalysisBase):
    """SPyRMSD class.

    This class is used to perform analysis on a trajectory.

    Parameters
    ----------
    mobile : AtomGroup
        Group of atoms for which the RMSD is calculated. If a trajectory is
        associated with the atoms then the computation iterates over the
        trajectory.
    reference : AtomGroup (optional)
        Group of reference atoms; if ``None`` then the current frame of
        `atomgroup` is used.
    ref_frame : int (optional)
        frame index to select frame from `reference`

    Attributes
    ----------
    universe: :class:`~MDAnalysis.core.universe.Universe`
        The universe to which this analysis is applied
    atomgroup: :class:`~MDAnalysis.core.groups.AtomGroup`
        The atoms to which this analysis is applied
    results: :class:`~MDAnalysis.analysis.base.Results`
        results of calculation are stored here, after calling
        :meth:`SPyRMSD.run`
    start: Optional[int]
        The first frame of the trajectory used to compute the analysis
    stop: Optional[int]
        The frame to stop at for the analysis
    step: Optional[int]
        Number of frames to skip between each analyzed frame
    n_frames: int
        Number of frames analysed in the trajectory
    times: numpy.ndarray
        array of Timestep times. Only exists after calling
        :meth:`SPyRMSD.run`
    frames: numpy.ndarray
        array of Timestep frame indices. Only exists after calling
        :meth:`SPyRMSD.run`
    
    Raises
    ------
    SelectionError
        If the selections from `atomgroup` and `reference` do not match.
    """

    def __init__(
        self,
        mobile: "AtomGroup",
        reference: "AtomGroup" = None,
        reference_frame: int = 0,
        **kwargs
    ):
        super().__init__(mobile.universe.trajectory, **kwargs)
        
        self.mobile_atoms = mobile
        self.ref_atoms = reference if reference is not None else self.mobile_atoms
        self.reference_frame = reference_frame

        if len(self.ref_atoms) != len(self.mobile_atoms):
            err = ("Reference and trajectory atom selections do "
                   "not contain the same number of atoms: "
                   f"N_ref={self.ref_atoms.n_atoms:d}, "
                   f"N_traj={self.mobile_atoms.n_atoms:d}")

            logger.exception(err)
            raise SelectionError(err)


        self.universe = universe_or_atomgroup.universe
        self.atomgroup = universe_or_atomgroup.select_atoms(select)

    def _prepare(self):
        """
        Prepare SPyRMSD analysis.

        SPyRMSD needs adjacency matrices for the molecular graphs of the
        selections and needs to perform graph isomorphism matching. The
        molecular graph is assumed to be constant, therefore graph
        isomorphisms can be computed here. 
        """
        current_frame = self.ref_atoms.universe.trajectory.ts.frame

        # Columns: frame, time, rmsd
        self.results.rmsd = np.zeros((self.n_frames, 3))

        self.ref_adj = adjacency_matrix(self.ref_atoms)
        self.ref_aprops = self.ref_atoms.types.copy()

        self.mobile_adj = adjacency_matrix(self.mobile_atoms)
        self.mobile_aprops = self.mobile_atoms.types.copy()

        # Compute isomorphisms at the first iteration
        self.isomorphisms = None

        # TODO: Check consistency of masses after graph isomorphism?

        try:
            self.ref_atoms.universe.trajectory[self.ref_frame]
            self._ref_coordinates64 = \
                self.ref_atoms.positions.copy().astype(np.float64)
        finally:
            # Move back to the original frame
            self.ref_atoms.universe.trajectory[current_frame]

        # Pre-allocate memory for the mobile coordinates
        self._mobile_coordinates64 = \
            self.mobile_atoms.positions.copy().astype(np.float64)

        # Pre-allocate memory for the results
        self.results.rmsd = np.zeros((self.n_frames, 3))

    def _single_frame(self):
        """Calculate data from a single frame of trajectory"""
        """
        Raises
        ------
        ImportError
             If the spyrmsd package is not installed
        """
        try:
            import spyrmsd.rmsd
        except ImportError:
            raise ImportError("""ERROR --- spyrmsd was not found!
                spyrmsd is required to compute symmetry-corrected RMSDs.
                try installing it using pip eg:
                    pip install spyrmsd
                or conda eg:
                    conda install spyrmsd -c conda-forge
                """)

        # Get current coordinates
        self._mobile_coordinates64[:] = self.mobile_atoms.positions

        # Get frame number and time for current timestep
        self.results.rmsd[self._frame_index, 0] = self._ts.frame
        self.results.rmsd[self._frame_index, 1] = self._trajectory.time

        # Compute minimum RMSD from graph isomorphisms
        min_rmsd, self.isomorphisms = spyrmsd.rmsd._rmsd_isomorphic_core(
            self._mobile_coordinates64,
            self._ref_coordinates64,
            self.mobile_aprops,
            self.ref_aprops,
            self.mobile_adj,
            self.ref_adj,
            center=False,
            minimize=False,
            isomorphisms=self.isomorphisms,
        )

        self.results.rmsd[self._frame_index, -1] = min_rmsd