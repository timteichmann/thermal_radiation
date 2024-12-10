# (c) 2024 Tim Teichmann

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

class Setup:
    def __init__(self, axisymmetric=False):
        self.points = []
        self.lines = []
        self.groups = []
        self.group2line = []
        self.group2temp = []
        self.group2eps = []
        self.species2temp = []
        self.group2species = []

        self.axisymmetric = axisymmetric

    def add_point(self, point):
        point = (round(point[0], 9), round(point[1], 9))
        self.points.append(np.array(point))
        return len(self.points) - 1

    def add_line(self, line, group):
        group_id = self.get_group(group)
        self.lines.append(line)
        line_id = len(self.lines) - 1
        self.group2line[group_id].append(line_id)
        return line_id

    def connect_point(self, point, group):
        point_id = self.add_point(point)
        self.add_line((point_id - 1, point_id), group)
        return point_id

    def add_group(self, group, temp, eps):
        if group in self.groups:
            raise ValueError("group " + str(group) + " already present")
        self.groups.append(group)
        self.group2line.append([])
        self.group2temp.append(temp)
        self.group2eps.append(eps)
        if not temp in self.species2temp:
            self.species2temp.append(temp)
        self.group2species.append(self.species2temp.index(temp))

    def get_group(self, group):
        if not group in self.groups:
            raise ValueError("group " + str(group) + " not found")
        return self.groups.index(group)

    def plot(self, annotate=False):
        for i, group in enumerate(self.groups):
            for j, line_id in enumerate(self.group2line[i]):
                line = self.lines[line_id]
                x1, y1 = self.points[line[0]]
                x2, y2 = self.points[line[1]]

                if j==0:
                    first_line = plt.plot((x1, x2), (y1, y2))
                else:
                    plt.plot((x1, x2), (y1, y2), color=first_line[0].get_color())

                if annotate:
                    plt.annotate(str(line_id + 1) + " - " + str(group), (x1 + (x2 - x1)/2, y1 + (y2 - y1)/2))

        plt.gca().set_aspect("equal")
        plt.show()

    def setup_case(self, input_file, surf_file, species_file, fnum, dt, n_steady, n_ave):
        # bounding box
        xmin = ymin = np.inf
        xmax = ymax = -np.inf
        for line in self.lines:
            x1, y1 = self.points[line[0]]
            x2, y2 = self.points[line[1]]

            xmin = min(xmin, x1, x2)
            xmax = max(xmax, x1, x2)
            ymin = min(ymin, y1, y2)
            ymax = max(ymax, y1, y2)

        with open(surf_file, "w") as f:
            def nl(text):
                f.write(text + "\n")
            nl("# Geometry definition:")
            nl("")
            nl(str(len(self.points)) + " points")
            nl(str(len(self.lines)) + " lines")
            nl("")
            nl("Points")
            nl("")
            for i, point in enumerate(self.points):
                nl(str(i + 1) + "\t" + str(point[0]) + "\t" + str(point[1]))
            nl("")
            nl("Lines")
            nl("")
            for i, line in enumerate(self.lines):
                nl(str(i + 1) + "\t" + str(line[0] + 1) + "\t" + str(line[1] + 1))

        with open(species_file, "w") as f:
            def nl(text):
                f.write(text + "\n")
            for specie, temp in enumerate(self.species2temp):
                nl("T" + str(specie + 1) + " -1 1e-26 0 0 0 0 0 1 0")

        with open(input_file, "w") as f:
            def nl(text):
                f.write(text + "\n")

            nl("# Simulation setup:")
            nl("variable            ns equal \"1e20\"")
            nl("variable            Ts equal \"300\"")
            nl("variable            ne equal \"1/4 * v_ns * sqrt(8*1.380649e-23*v_Ts/PI/1e-26)\"")
            nl("")

            nl("# Variables:")
            for specie, temp in enumerate(self.species2temp):
                nl("variable            T" + str(specie + 1) + " equal \"" + str(temp) + "\"")
            nl("")

            for group, temp, eps in zip(self.groups, self.group2temp, self.group2eps):
                nl("variable            e_" + str(group) + " equal \"" + str(eps) + "\"")
                nl("variable            n_" + str(group) + " equal \"v_e_" + str(group) + " * v_ns\"")
                nl("")

            nl("# Geometry:")
            nl("seed                12345") # reproducible for now
            nl("dimension           2")
            nl("global              gridcut 0.0 comm/sort yes surfmax 10000 splitmax 1000")
            nl("")
            if self.axisymmetric:
                nl("boundary            oo ao pp")
            else:
                nl("boundary            oo oo pp")
            nl("create_box          " + str(xmin) + " " + str(xmax) + " " + str(ymin) + " " + str(ymax) + " -0.5 0.5")
            nl("create_grid         1 1 1")
            nl("")
            nl("# Surfaces:")
            nl("read_surf           " + str(surf_file))
            for group, lines in zip(self.groups, self.group2line):
                surf_ids = ""
                for line in lines:
                    surf_ids += " " + str(line + 1)
                nl("group               " + str(group) + " surf id" + surf_ids)

            nl("")
            nl("surf_collide        all adiabatic")
            for group, eps in zip(self.groups, self.group2eps):
                nl("surf_react          " + str(group) + " global ${e_" + str(group) + "} 0.0")
                nl("surf_modify         " + str(group) + " collide all react " + str(group))

            nl("")
            nl("# Emits:")
            species = ""
            for specie in range(len(self.species2temp)):
                species += " T" + str(specie + 1)
            nl("species             " + str(species_file) + species)

            nl("")
            for group, specie, eps in zip(self.groups, self.group2species, self.group2eps):
                # only add mixture and fix emit for groups with eps > 0
                if not eps == 0.0:
                    nl("mixture             " + str(group) + " T" + str(specie + 1) + " nrho ${n_" + str(group) + "} temp ${T" + str(specie + 1) + "}")
                    nl("fix                 " + str(group) + " emit/surf " + str(group) + " " + str(group))

            nl("")
            nl("# Steady-state run:")
            nl("global              fnum " + str(fnum))
            nl("timestep            " + str(dt))
            nl("stats               1")
            nl("run                 " + str(int(n_steady)))
            nl("reset_timestep      0")

            nl("")
            nl("# Computes:")
            nl("variable            SB equal \"" + str(5.670374419e-8) + "\"")
            for group in self.groups:
                nl("")
                nl("compute             N_" + str(group) + "_ surf " + str(group) + " species nflux norm flow")
                nl("fix                 N_" + str(group) + "_ ave/surf " + str(group) + " 1 1 1 c_N_" + str(group) + "_[*] ave running")
                nl("compute             N_" + str(group) + " reduce sum f_N_" + str(group) + "_[*]")
                sum_term = ""
                for specie in range(len(self.species2temp)):
                    if not sum_term == "":
                        sum_term += " + "
                    sum_term += "c_N_" + str(group) + "[" + str(specie + 1) + "]*v_T" + str(specie + 1) + "^4"

                nl("variable            Q_" + str(group) + " equal \"v_SB * (" + sum_term + ")/v_ne\"")

            nl("")
            nl("# Averaging run:")
            heat_flows = ""
            for group in self.groups:
                heat_flows += " v_Q_" + str(group)
            nl("stats_style         step np" + heat_flows)
            nl("run                 " + str(n_ave))
