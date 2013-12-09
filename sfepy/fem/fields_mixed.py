"""
Notes
-----

Important attributes of continuous (order > 0) :class:`Field` and
:class:`SurfaceField` instances:

- `vertex_remap` : `econn[:, :n_vertex] = vertex_remap[conn]`
- `vertex_remap_i` : `conn = vertex_remap_i[econn[:, :n_vertex]]`

where `conn` is the mesh vertex connectivity, `econn` is the
region-local field connectivity.
"""
# import time
# import numpy as nm
# 
# from sfepy.base.base import output, assert_
from sfepy.base.base import Struct
import fea
# from sfepy.fem.utils import prepare_remap
# from sfepy.fem.dof_info import expand_nodes_to_dofs
# from sfepy.fem.global_interp import get_ref_coors
from sfepy.fem.facets import get_facet_dof_permutations
from sfepy.fem.fields_base import Field, VolumeField, SurfaceField
from sfepy.fem.fields_nodal import H1NodalMixin 
# from sfepy.fem.extmods.bases import evaluate_in_rc


class Hdiv_RaviartThomas(H1NodalMixin, VolumeField, Struct): #H1NodalMixin,
    family_name = 'volume_Hdiv_RaviartThomas'

    def _setup_approximations(self):
        self.aps = {}
        self.aps_by_name = {}
        for ig in self.igs:
            name = self.interp.name + '_%s_ig%d' % (self.region.name, ig)
            ap = fea.HdivApproximation(name, self.interp, self.region, ig)
            self.aps[ig] = ap
            self.aps_by_name[ap.name] = ap

    def setup_dof_conns(self, dof_conns, dpn, dc_type,
                            region, is_trace=False):
        """Setup dof connectivities of various kinds as needed by terms."""
        dpn = 1 # number of components enforced to 1
        VolumeField.setup_dof_conns(self, dof_conns, dpn, dc_type, 
                                    region, is_trace=False)
        print 'JV: field_mixed.Hdiv_RaviartThomas.setup_dof_conns'

    def _setup_facet_orientations(self):
        """
        This function is necessary here since the commit on 2013/10/11
        in sfepy/fem/facets.py
        """
        order = self.approx_order
        self.node_desc = self.interp.describe_nodes()

        edge_nodes = self.node_desc.edge_nodes
        if edge_nodes is not None:
            n_fp = self.gel.edges.shape[1]
            self.edge_dof_perms = get_facet_dof_permutations(n_fp, self.igs,
                                                             order+2)

        face_nodes = self.node_desc.face_nodes
        if face_nodes is not None:
            raise NotImplementedError('JV: face_dof_perms is not implemented \
                                      yet.!')
#             n_fp = self.gel.faces.shape[1]
#             self.face_dof_perms = get_facet_dof_permutations(n_fp, self.igs,
#                                                              order)
