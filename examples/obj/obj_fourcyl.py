from vfvchem_geom2d import Geom2D

# Channel with four cylinders in a 2 x 2 arrangement
# The cylinders cut holes in the channel and are also domains themselves
#
# Boundaries and Interfaces
# bnd_i - inlet
# bnd_o - outlet
# bnd_w - bottom and top walls of channel
# bnd_h - cylinder walls (same curves as itf_h; e.g., no-slip for momentum transfer)
# itf_h - cylinder walls (between dom_c and dom_h)
#
# Domains
# dom_c - channel
# dom_h - cylinders
#
# Group names must begin with dom_, bnd_, or itf_
# The ids shown before grouping are only used to select curves and domains
#

# define geometry
# the cylinders overlap the channel; the channel will be automatically cut around them
geom = Geom2D("obj_fourcyl")
rec1 = geom.add_rectangle(0.0, 0.0, 0.2, 0.05)  # arguments: x_min, y_min, x_max, y_max
ell1 = geom.add_ellipse(0.1125, 0.0125, 0.005, 0.005)  # arguments: x_center, y_center, x_radius, y_radius
ell2 = geom.add_ellipse(0.1375, 0.0125, 0.005, 0.005)  # arguments: x_center, y_center, x_radius, y_radius
ell3 = geom.add_ellipse(0.1125, 0.0375, 0.005, 0.005)  # arguments: x_center, y_center, x_radius, y_radius
ell4 = geom.add_ellipse(0.1375, 0.0375, 0.005, 0.005)  # arguments: x_center, y_center, x_radius, y_radius
geom.finalize_geometry(0.0025)  # arguments: size_curve (maximum segment length on curved edges)

# show geometry with domain, boundary, and interface ids
# this will show the automatically generated ids
geom.show_geometry()  # call to show geometry with ids

# name and/or join ids
# arguments: group name, list of ids to join
# an interface curve may also be in one boundary group
geom.group_dom("dom_c", [4])  # channel
geom.group_dom("dom_h", [0, 1, 2, 3])  # cylinders
geom.group_bnd("bnd_i", [1])  # inlet
geom.group_bnd("bnd_o", [2])  # outlet
geom.group_bnd("bnd_w", [0, 3])  # bottom and top walls
geom.group_bnd("bnd_h", [4, 5, 6, 7])  # cylinder walls
geom.group_itf("itf_h", [4, 5, 6, 7])  # cylinder walls
geom.update_group()  # call to check the groups

# show group names
geom.show_geometry()

# export geometry
# arguments: file name, h_min (optional; boundary node spacing used in vfvm2d)
# h_min is only used to check that segments are long enough for vfvm2d; it does not modify the geometry
geom.export_geometry("obj_fourcyl.obj", 0.001)
