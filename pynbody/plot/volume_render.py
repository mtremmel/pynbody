import pylab as p
import matplotlib
import numpy as np
from .. import sph, config
from .. import units as _units

class RenderVolume(object):
	def __init__(self, sim, resolution=500, width=None):
		sim.physical_units() # make sure things are in physical units!
		self.sim = sim
		self.resolution = resolution
		self.width = width

	def density(self, family=None, vmin=None, vmax=None, dynamic_range=4, log=True, color=None, colortable=None, create_figure=True):
		import mayavi
		from mayavi import mlab
		from tvtk.util.ctf import PiecewiseFunction, ColorTransferFunction
		import palettable

		if create_figure:
			fig = mlab.figure(size=(500, 500), bgcolor=(0, 0, 0))

		ss = self.sim
		if family == 'gas':
			ss = self.sim.g
		if family == 'dm':
			ss = self.sim.dm
		if family == 'star':
			ss = self.sim.s

		grid_data = sph.to_3d_grid(ss, qty='rho', nx=self.resolution,
		                           x2=None if self.width is None else self.width / 2)

		if log:
			grid_data = np.log10(grid_data)
			if vmin is None:
				vmin = grid_data.max() - dynamic_range
			if vmax is None:
				vmax = grid_data.max()

		else:
			if vmin is None:
				vmin = np.min(grid_data)
			if vmax is None:
				vmax = np.max(grid_data)

		grid_data[grid_data < vmin] = vmin
		grid_data[grid_data > vmax] = vmax

		otf = PiecewiseFunction()
		otf.add_point(vmin, 0.0)
		otf.add_point(vmax, 1.0)

		sf = mayavi.tools.pipeline.scalar_field(grid_data)
		V = mlab.pipeline.volume(sf, color=color, vmin=vmin, vmax=vmax)

		V.trait_get('volume_mapper')['volume_mapper'].blend_mode = 'maximum_intensity'

		if color is None:
			ctf = ColorTransferFunction()
			#unless color is specified, use cubehelix map
			if colortable is None:
				colortable = palettable.cubehelix.cubehelix1_16.colors
			colortable = np.array(colortable)

			vbins = np.arange(vmin,vmax,(vmax-vmin)/len(colortable))
			for i in range(len(colortable)):
				ctf.add_rgb_point(vbins[i],colortable[i][0]/256.,colortable[i][1]/256., colortable[i][2]/256.)

			V._volume_property.set_color(ctf)
			V._ctf = ctf
			V.update_ctf = True

		V._otf = otf
		V._volume_property.set_scalar_opacity(otf)

		return V

	def gas_density(self,**kwargs):
		return self.density(family='gas', **kwargs)

	def dm_density(self,**kwargs):
		return self.density(family='dm', **kwargs)

	def star_density(self, **kwargs):
		return self.density(family='star', **kwargs)
