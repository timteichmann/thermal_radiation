# (c) 2024 Tim Teichmann

R = 0.773
H = 0.892
r = R/2
h = 0.01784
t = 1e-9

from Setup import *
S = Setup(axisymmetric=True)

# create the groups
S.add_group("P1", 300, 0.9)
S.add_group("P2", 300, 0.9)
S.add_group("P3R", 4.5, 0.1)
S.add_group("P3L", 4.5, 0.1)
S.add_group("P4R", 4.5, 0.1)
S.add_group("P4L", 4.5, 0.1)
S.add_group("P34T", 4.5, 0.0)
S.add_group("R1", 300, 0.9)

S.add_point((-h-t, 0.0))
S.connect_point((-h-t, r), "P3L")
S.connect_point((-h+t, r), "P34T")
S.connect_point((-h+t, 0.0), "P3R")

S.add_point((h-t, 0.0))
S.connect_point((h-t, r), "P4L")
S.connect_point((h+t, r), "P34T")
S.connect_point((h+t, 0.0), "P4R")

S.add_point((H, 0.0))
S.connect_point((H, R), "P2")
S.connect_point((-H, R), "R1")
S.connect_point((-H, 0.0), "P1")

S.setup_case("in.case2", "geometry.surf", "photon.species", 1e16, 1e-2, 100, 1000)
