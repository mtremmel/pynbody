import pylab as p
import matplotlib
import numpy as np
from .. import sph, config
from .. import units as _units
from .. import filt

class RenderVolume(object):
	def __init__(self, sim, resolution=500, width=None):
		sim.physical_units() # make sure things are in physical units!
		self.sim = sim
		self.resolution = resolution
		self.width = width
		self._loaded_data = {} #store calculated grids for faster plotting

	def _create_colormap(self, colortable, vbins):
		from tvtk.util.ctf import ColorTransferFunction

		ctf = ColorTransferFunction()

		for i in range(len(colortable)):
			ctf.add_rgb_point(vbins[i], colortable[i][0] / 256., colortable[i][1] / 256., colortable[i][2] / 256.)

		return ctf

	def density(self, family=None, vmin=None, vmax=None, dynamic_range=4,
	            log=True, color=None, colortable=None, create_figure=True):
		import mayavi
		from mayavi import mlab
		from tvtk.util.ctf import PiecewiseFunction, ColorTransferFunction
		import palettable

		if create_figure:
			fig = mlab.figure(size=(500, 500), bgcolor=(0, 0, 0))

		data_name = 'all_den'
		ss = self.sim
		if family == 'gas':
			data_name = 'gas_den'
			ss = self.sim.g
		if family == 'dm':
			data_name = 'dm_den'
			ss = self.sim.dm
		if family == 'star':
			data_name = 'star_den'
			ss = self.sim.s

		if data_name in self._loaded_data.keys():
			print("using previously calculated data grid")
			grid_data = self._loaded_data[data_name]
		else:
			grid_data = sph.to_3d_grid(ss, qty='rho', nx=self.resolution,
		                           x2=None if self.width is None else self.width / 2)
			self._loaded_data[data_name] = grid_data

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
			if colortable is None: #default colormap is cubehelix
				colortable = np.array(palettable.cubehelix.cubehelix1_16.colors)
			vbins = np.arange(vmin, vmax, (vmax - vmin) / len(colortable))
			ctf = self._create_colormap(colortable,vbins)
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

	def star_tform(self,min_age=None, max_age=None, log=False, maxstarsize=0.5,
	            color=None, colortable=None, create_figure=True, age_bins=None):

		import mayavi
		from mayavi import mlab
		from tvtk.util.ctf import PiecewiseFunction, ColorTransferFunction
		import palettable

		data_name = 'star_tform'

		if create_figure:
			fig = mlab.figure(size=(500, 500), bgcolor=(0, 0, 0))


		if starsize is not None:
			smf = filt.HighPass('smooth', str(starsize) + ' kpc')
			sim.s[smf]['smooth'] = array.SimArray(starsize, 'kpc', sim=self.sim)

		if data_name in self._loaded_data.keys():
			print("using previously calculated data grid")
			grid_data = self._loaded_data[data_name]
		else:
			grid_data = sph.to_3d_grid(self.sim.s, qty='tform', nx=self.resolution, snap_slice=filt.HighPass('tform',0),
		                           x2=None if self.width is None else self.width / 2)
			self._loaded_data[data_name] = grid_data

		grid_data = grid_data.in_units('Gyr')

		sim_time = self.sim.properties['time'].in_units('Gyr')


		if max_age is None:
			vmin=0
		else:
			vmin = sim_time-max_age

		if min_age is None:
			vmax = sim_time
		else:
			vmax = sim_time - min_age

		if log is True:
			if vmin <= 0:
				vmin = 0.001
			grid_data = np.log10(grid_data)
			vmax = np.log10(vmax)

		grid_data[(grid_data>vmax)] = vmax
		grid_data[(grid_data<vmin)] = vmin

		otf = PiecewiseFunction()
		otf.add_point(vmin, 0.0)
		otf.add_point(0.5, 1.0)

		sf = mayavi.tools.pipeline.scalar_field(grid_data)
		V = mlab.pipeline.volume(sf, color=color, vmin=vmin, vmax=vmax)

		V.trait_get('volume_mapper')['volume_mapper'].blend_mode = 'maximum_intensity'

		if color is None:
			if colortable is None: #default colormap is BlueOrange10
				colortable = np.array(palettable.lightbartlein.diverging.BlueOrange10_6.colors)[::-1]
			if age_bins is None:
				age_bins = np.array([0.01, 0.1, 1.0, 2.0, 4.0, 10.0])
			else:
				age_bins = np.array(age_bins) #make sure age_bins is an array
			vbins = sim_time - age_bins
			vbins = vbins[::-1]
			if log is True:
				vbins = np.log10(vbins)
			ctf = self._create_colormap(colortable,vbins)
			V._volume_property.set_color(ctf)
			V._ctf = ctf
			V.update_ctf = True

		V._otf = otf
		V._volume_property.set_scalar_opacity(otf)

		return V

