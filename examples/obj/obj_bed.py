from vfvchem_geom2d import Geom2D

# Channel (0.1 x 0.02) with a packed bed in the lower half under a clear fluid layer
#
# Boundaries and Interfaces
# bnd_lb - bed inlet (left, y < 0.01)
# bnd_lf - fluid inlet (left, y > 0.01)
# bnd_rb - bed outlet (right, y < 0.01)
# bnd_rf - fluid outlet (right, y > 0.01)
# bnd_b - bottom wall
# bnd_t - top wall
# itf_bed - bed-fluid interface (y = 0.01)
#
# Domains
# dom_bed - packed bed
# dom_fluid - clear fluid
#
# Group names must begin with dom_, bnd_, or itf_
# The ids shown before grouping are only used to select curves and domains
#

# define geometry
geom = Geom2D("obj_bed")
reg1 = geom.add_polygon([(0.0, 0.0), (0.1, 0.0), (0.1, 0.01), (0.0, 0.01)])  # packed bed
reg2 = geom.add_polygon([(0.0, 0.01), (0.1, 0.01), (0.1, 0.02), (0.0, 0.02)])  # clear fluid
geom.finalize_geometry()  # call at the end of geometry definition

# show geometry with domain, boundary, and interface ids
# this will show the automatically generated ids
geom.show_geometry()  # call to show geometry with ids

# name and/or join ids
# arguments: group name, list of ids to join
geom.group_dom("dom_bed", [0])  # packed bed
geom.group_dom("dom_fluid", [1])  # clear fluid
geom.group_bnd("bnd_b", [0])  # bottom wall
geom.group_bnd("bnd_rb", [1])  # bed outlet
geom.group_bnd("bnd_lb", [2])  # bed inlet
geom.group_bnd("bnd_rf", [3])  # fluid outlet
geom.group_bnd("bnd_t", [4])  # top wall
geom.group_bnd("bnd_lf", [5])  # fluid inlet
geom.group_itf("itf_bed", [6])  # bed-fluid interface
geom.update_group()  # call to check the groups

# show group names
geom.show_geometry()

# export geometry
# arguments: file name, h_min (optional; boundary node spacing used in vfvm2d)
# h_min is only used to check that segments are long enough for vfvm2d; it does not modify the geometry
geom.export_geometry("obj_bed.obj", 5e-4)
