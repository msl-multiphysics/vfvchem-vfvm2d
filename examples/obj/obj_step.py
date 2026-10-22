from vfvchem_geom2d import Geom2D

# Channel with a backward-facing step
#
# Boundaries
# bnd_i - inlet
# bnd_o - outlet
# bnd_w - walls
#
# Domains
# dom_f - fluid
#
# Group names must begin with dom_, bnd_, or itf_
# The ids shown before grouping are only used to select curves and domains
#

# define geometry
geom = Geom2D("obj_step")
reg1 = geom.add_polygon(  # channel
    [
        (0.0, 0.2),
        (0.4, 0.2),
        (0.4, 0.0),
        (2.4, 0.0),
        (2.4, 0.4),
        (0.0, 0.4),
    ],
)
geom.finalize_geometry()  # call at the end of geometry definition

# show geometry with domain, boundary, and interface ids
# this will show the automatically generated ids
geom.show_geometry()  # call to show geometry with ids

# name and/or join ids
# arguments: group name, list of ids to join
geom.group_dom("dom_f", [0])  # fluid
geom.group_bnd("bnd_i", [5])  # inlet
geom.group_bnd("bnd_o", [3])  # outlet
geom.group_bnd("bnd_w", [0, 1, 2, 4])  # walls
geom.update_group()  # call to check the groups

# show group names
geom.show_geometry()

# export geometry
# arguments: file name, h_min (optional; boundary node spacing used in vfvm2d)
# h_min is only used to check that segments are long enough for vfvm2d; it does not modify the geometry
geom.export_geometry("obj_step.obj", 0.02)
