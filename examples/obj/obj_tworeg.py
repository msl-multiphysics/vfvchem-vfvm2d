from vfvchem_geom2d import Geom2D

# Square split along its diagonal into two triangles
#
# Boundaries and Interfaces
# bnd_b - bottom
# bnd_r - right
# bnd_t - top
# bnd_l - left
# itf_d - diagonal interface
#
# Domains
# dom_lo - lower triangle
# dom_up - upper triangle
#
# Group names must begin with dom_, bnd_, or itf_
# The ids shown before grouping are only used to select curves and domains
#

# define geometry
geom = Geom2D("obj_tworeg")
reg1 = geom.add_polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])  # lower triangle
reg2 = geom.add_polygon([(0.0, 0.0), (1.0, 1.0), (0.0, 1.0)])  # upper triangle
geom.finalize_geometry()  # call at the end of geometry definition

# show geometry with domain, boundary, and interface ids
# this will show the automatically generated ids
geom.show_geometry()  # call to show geometry with ids

# name and/or join ids
# arguments: group name, list of ids to join
geom.group_dom("dom_lo", [0])  # lower triangle
geom.group_dom("dom_up", [1])  # upper triangle
geom.group_bnd("bnd_b", [0])  # bottom
geom.group_bnd("bnd_r", [1])  # right
geom.group_bnd("bnd_t", [2])  # top
geom.group_bnd("bnd_l", [3])  # left
geom.group_itf("itf_d", [4])  # diagonal interface
geom.update_group()  # call to check the groups

# show group names
geom.show_geometry()

# export geometry
# arguments: file name, h_min (optional; boundary node spacing used in vfvm2d)
# h_min is only used to check that segments are long enough for vfvm2d; it does not modify the geometry
geom.export_geometry("obj_tworeg.obj", 0.02)
