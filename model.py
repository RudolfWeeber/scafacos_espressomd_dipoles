"""This is the base class for a simple ferrogel model. It contains the ESPResSo system instance as a member.
"""


import espressomd
from espressomd.interactions import HarmonicBond
from espressomd.magnetostatics import DipolarDirectSumCpu, Scafacos

from time import time
import sys
import pickle
import numpy as np
from random import shuffle


class GelModel:

    # Lennard Jones interaction parameters
    lj_sigma = 1.
    lj_cut = 2**(1. / 6) * lj_sigma
    lj_eps = 1.

    # Harmonic bond parameters
    bond_strength = 200

    # Langevin thermostat parameters
    kT = 1
    gamma = 1

    # Dipole parameters
    dip_lambda = 2.
    mu = np.sqrt(lj_sigma**3 * dip_lambda)

    # The gel network is constructed by first placing beads as nodes on a
    # regular cubic lattice and then connecting them by chains of beads.
    # To produce a roughly spherical sample, nodes further than a given radius
    # from the center are not placed.

    # Diameter of gel sphere in units of node-node spacings
    node_grid_diameter = 10

    # Length of the bead chains between the nodes
    chain_length = 14
    node_spacing = (1 + chain_length) * lj_cut

    # Through diffusion and shape change, the gel can move outside the simulation box.
    # To prevent this, the box is resized and/or the gel is moved if there is
    # not at least a space of lower_margin on all sides.
    # After a resize, a space of margin will be around the gel at all sides.
    margin = 4
    lower_margin = 2

    # Configuration for the initial relaxation of the gel network
    warmup_steps_without_dipolar = 20000

    # Integrator: Time step
    dt = 0.015

    # Initial configuration of the domain decomposition and verlet lists
    skin = 0.4

    # Setup espresso system
    s = espressomd.System(box_l=[1,1,1]) # Resizing occurs later

    def node_positions(self):
        """Returns the positions of the node particles as n-by-3 array"""
        lim = int(self.node_grid_diameter / 2)
        positions = []
        for i in range(-lim, lim + 1, 1):
            for j in range(-lim, lim + 1, 1):
                for k in range(-lim, lim + 1, 1):
                    if i**2 + j**2 + k**2 >= (self.node_grid_diameter / 2)**2:
                        continue
                    positions.append(np.array((i, j, k)) * self.node_spacing)
        positions = np.array(positions)

        return positions

    def connect_nodes(self, n1, n2):
        """Adds the chain between nodes with particle ids n1 and n2"""
        # Normalized distance between the nodes
        d = self.s.distance_vec(self.s.part[n1], self.s.part[
                                n2]) / self.s.distance(self.s.part[n1], self.s.part[n2])

        # Store particle of previous iteration
        previous_p = None

        # Loop for the particles in the chain
        for k in range(1, self.chain_length + 1, 1):

            # Place particle
            p = self.s.part.add(pos=self.s.part[n1].pos + d * self.lj_cut * k)

            if k == 1:
                # Bind first particle to node n1
                p.add_bond((self.bond, n1))
            else:
                # Bind particle to previous particle in the chain
                p.add_bond((self.bond, previous_p))

            # Remember handle of particle placed in this round
            previous_p = p

            if k == self.chain_length:
                # Bind last particle in the chain to node n2
                p.add_bond((self.bond, n2))

    def verify_initial_state(self):
        """Verifies various aspects of the initial state"""
        # Short-hand for the Espresso system handle
        s = self.s

        # Number of particles
        assert len(s.part) == self.n_nodes + self.n_chains * self.chain_length

        # Energy close to 0
        assert abs(s.analysis.energy()["total"]) < 1e-3

        # Dipole moment per particle
        assert all(np.abs(np.sum(s.part[:].dip**2, 1) - self.mu**2) < 1E-5)

        # No bonds on nodes
        for p in s.part[:self.n_nodes]:
            assert p.bonds == ()

        # Bonds on chains
        for i in range(self.n_chains):
            # First particle bound to node
            p0 = s.part[self.n_nodes + i * self.chain_length]
            assert len(p0.bonds) == 1
            assert p0.bonds[0][1] < self.n_nodes

            # Middle of chain has bond to previous chain particle
            for j in range(1, self.chain_length - 1, 1):
                p = s.part[self.n_nodes + i * self.chain_length + j]
                assert len(p.bonds) == 1
                assert p.bonds[0][1] == p.id - 1

            # LaSt particle in the chain has bond to previous chain particle
            # and next node
            p = s.part[self.n_nodes + (i + 1) * self.chain_length - 1]
            assert len(p.bonds) == 2
            assert p.bonds[0][1] == p.id - 1
            assert p.bonds[1][1] < self.n_nodes
            assert p.bonds[1][1] != p0.bonds[0][1]

    def adjust_box_and_shift_system(self, force=False):
        """Adjust box and shift particles

           If there is less than self.margin free space at any side of the
           gel, shift it.
           If the difference between the gel extensions and the box length is
           less than 2*self.margin ore more than 3*self.margin, resize the
           simulation box.
           If force=True, the resize will be done in any case.

        """
        # Short-hand for Espresso system handle
        s = self.s

        # Bounds of the gel
        lower_left_front = np.amin(s.part[:].pos, 0)
        upper_right_back = np.amax(s.part[:].pos, 0)

        # Cache particle positions
        pos = s.part[:].pos

        # Decide whether to move
        move = False
        if np.any(pos < self.lower_margin):
            move = True
        if np.any((s.box_l - pos) < self.lower_margin):
            move = True

        # Extent of the gel
        v_diff = upper_right_back - lower_left_front

        # Decide whether to resize the box
        resize = False
        if force:
            move = True
            resize = True
            print("Forced box resize")
        if any(s.box_l - v_diff < 2 * self.lower_margin):
            resize = True
            move = True
            print("Box too small. Resizing")
        if any(s.box_l - v_diff > 3 * self.margin):
            resize = True
            move = True
            print("Box too big. Resizing")

        # Execute resize and move if needed
        if resize:
            s.box_l = v_diff + 2 * self.margin
            print("box", s.box_l)
        if move:
            shift = s.box_l / 2 - lower_left_front - v_diff / 2
            print("Moving particles by", shift)
            s.part[:].pos = s.part[:].pos + shift

    def setup_system(self):
        """Setup the Espresso system parameters."""

        # Shorthand for Espresso system handle
        s = self.s

        # Positions and number of nodes of the gel network
        self.node_pos = self.node_positions()
        self.n_nodes = len(self.node_pos)

        # Box size is highest coordinate of any node +1/2 node_spacing
        l = np.amax(self.node_pos) + self.node_spacing / 2
        s.box_l = l, l, l

        # Open boundary conditions
        s.periodicity = 0, 0, 0

        # Thermostat
        s.thermostat.set_langevin(kT=self.kT, gamma=self.gamma)

        # Integrator
        s.time_step = self.dt
        s.cell_system.skin = self.skin

        # Lennard Jones interaction
        s.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=self.lj_sigma, epsilon=self.lj_eps, cutoff=self.lj_cut, shift="auto")

        # Harmonic bond
        self.bond = HarmonicBond(r_0=self.lj_cut, k=self.bond_strength)
        s.bonded_inter.add(self.bond)

    def setup_particles(self, warmup=True):
        """Places the beads of the gel network and cross-links them."""
        # Short-hand for Espresso system handle
        s = self.s

        # Place the node particles
        s.part.add(
            id=np.arange(self.node_pos.shape[0], dtype=int), pos=self.node_pos)

        # Cache references to some class attributes to speed up the nested loop
        dist_func = self.s.distance
        n_nodes = self.n_nodes
        node_spacing = self.node_spacing
        part = self.s.part

        # Iterate over all pairs of node particles and connect those that are
        # neighbors (one node_spacing apart)
        self.n_chains = 0
        for i in range(n_nodes):
            for j in range(i + 1, n_nodes, 1):
                if abs(dist_func(part[i], part[j]) - node_spacing) < 1E-7:
                    self.connect_nodes(i, j)
                    self.n_chains += 1

        # Randomize dipole moment orientation
        cos_theta = 2 * np.random.random(len(s.part)) - 1
        sin_theta = np.sin(np.arccos(cos_theta))
        phi = 2 * np.pi * np.random.random(len(s.part))
        dip = np.zeros((len(s.part), 3))
        dip[:, 0] = self.mu * sin_theta * np.sin(phi)
        dip[:, 1] = self.mu * sin_theta * np.cos(phi)
        dip[:, 2] = self.mu * cos_theta
        s.part[:].dip = dip

        self.verify_initial_state()
