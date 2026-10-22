from vfvchem_geom2d import Geom2D

# Square with bottom-left and top-right corner rectangles
#
# Boundaries and Interfaces
# bnd_bl - bottom-left outer (left + bottom)
# bnd_mb - middle outer (bottom + right)
# bnd_tr - top-right outer (right + top)
# bnd_mt - middle outer (top + left)
# itf_bl - bottom-left interface
# itf_tr - top-right interface
#
# Domains
# dom_m - middle
# dom_bl - bottom-left
# dom_tr - top-right
#
# Group names must begin with dom_, bnd_, or itf_
# The ids shown before grouping are only used to select curves and domains
#

# define geometry
# rec2 and rec3 overlap with rec1; rec1 will be automatically cut into three regions
geom = Geom2D("obj_threereg")
rec1 = geom.add_rectangle(0, 0, 1, 1)  # arguments: x_min, y_min, x_max, y_max
rec2 = geom.add_rectangle(0, 0, 0.75, 0.25)  # arguments: x_min, y_min, x_max, y_max
rec3 = geom.add_rectangle(0.25, 0.75, 1, 1)  # arguments: x_min, y_min, x_max, y_max
geom.finalize_geometry()  # call at the end of geometry definition

# show geometry with domain, boundary, and interface ids
# this will show the automatically generated ids
geom.show_geometry()  # call to show geometry with ids

# name and/or join ids
# arguments: group name, list of ids to join
geom.group_dom("dom_m", [2])  # middle
geom.group_dom("dom_bl", [0])  # bottom-left
geom.group_dom("dom_tr", [1])  # top-right
geom.group_bnd("bnd_bl", [4, 5])  # bottom-left outer (left + bottom)
geom.group_bnd("bnd_mb", [3, 2])  # middle outer (bottom + right)
geom.group_bnd("bnd_tr", [7, 6])  # top-right outer (right + top)
geom.group_bnd("bnd_mt", [1, 0])  # middle outer (top + left)
geom.group_itf("itf_bl", [8, 9])  # bottom-left interface
geom.group_itf("itf_tr", [10, 11])  # top-right interface
geom.update_group()  # call to check the groups

# show group names
geom.show_geometry()

# export geometry
# arguments: file name, h_min (optional; boundary node spacing used in vfvm2d)
# h_min is only used to check that segments are long enough for vfvm2d; it does not modify the geometry
geom.export_geometry("obj_threereg.obj", 0.02)
