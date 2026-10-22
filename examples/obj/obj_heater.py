from vfvchem_geom2d import Geom2D

# Channel with bump surrounded by solid regions.
#
# Boundaries and Interfaces
# bnd_i - inlet
# bnd_o - outlet
# bnd_b - bottom solid
# bnd_t - top solid
# bnd_bc - bottom channel wall (same curves as itf_b; e.g., no-slip for momentum transfer)
# bnd_tc - top channel wall (same curves as itf_t; e.g., no-slip for momentum transfer)
# itf_b - bottom solid-channel interface
# itf_t - top solid-channel interface
#
# Domains
# dom_c - channel
# dom_b - bottom solid
# dom_t - top solid
#
# Group names must begin with dom_, bnd_, or itf_
# The ids shown before grouping are only used to select curves and domains
#

# define geometry
geom = Geom2D("obj_heater")
reg1 = geom.add_polygon(  # channel
    [
        (0.00, 0.00),
        (1.50, 0.00),
        (1.50, 0.50),
        (2.00, 0.50),
        (2.00, 0.00),
        (3.50, 0.00),
        (3.50, 0.50),
        (2.50, 0.50),
        (2.50, 1.00),
        (1.00, 1.00),
        (1.00, 0.50),
        (0.00, 0.50),
    ],
)
reg2 = geom.add_polygon(  # bottom
    [
        (0.00, 0.00),
        (1.50, 0.00),
        (1.50, 0.50),
        (2.00, 0.50),
        (2.00, 0.00),
        (3.50, 0.00),
        (3.50, -0.50),
        (0.00, -0.50),
    ],
)
reg3 = geom.add_polygon(  # top
    [
        (0.00, 0.50),
        (1.00, 0.50),
        (1.00, 1.00),
        (2.50, 1.00),
        (2.50, 0.50),
        (3.50, 0.50),
        (3.50, 1.50),
        (0.00, 1.50),
    ],
)
geom.finalize_geometry()  # call at the end of geometry definition

# show geometry with domain, boundary, and interface ids
# this will show the automatically generated ids
geom.show_geometry()  # call to show geometry with ids

# name and/or join ids
# arguments: group name, list of ids to join
# an interface curve may also be in one boundary group
geom.group_dom("dom_c", [0])  # channel
geom.group_dom("dom_b", [1])  # bottom solid
geom.group_dom("dom_t", [2])  # top solid
geom.group_bnd("bnd_i", [1])  # inlet
geom.group_bnd("bnd_o", [0])  # outlet
geom.group_bnd("bnd_b", [4, 3, 2])  # bottom solid
geom.group_bnd("bnd_t", [7, 6, 5])  # top solid
geom.group_bnd("bnd_bc", [8, 9, 10, 11, 12])  # bottom channel wall
geom.group_bnd("bnd_tc", [17, 16, 15, 14, 13])  # top channel wall
geom.group_itf("itf_b", [8, 9, 10, 11, 12])  # bottom solid-channel interface
geom.group_itf("itf_t", [17, 16, 15, 14, 13])  # top solid-channel interface
geom.update_group()  # call to check the groups

# show group names
geom.show_geometry()

# export geometry
# arguments: file name, h_min (optional; boundary node spacing used in vfvm2d)
# h_min is only used to check that segments are long enough for vfvm2d; it does not modify the geometry
geom.export_geometry("obj_heater.obj", 0.05)
