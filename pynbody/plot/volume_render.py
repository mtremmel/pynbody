import pylab as p
import matplotlib
import numpy as np
from .. import sph, config
from .. import units as _units
from .. import filt, array

class RenderVolume(object):
	def __init__(self, sim, resolution=500):
		sim.physical_units() # make sure things are in physical units!
		self.sim = sim
		self.resolution = resolution
		self._loaded_data = {} #store calculated grids for faster plotting
		self._starsize=0.25 #default starsize

	def _create_colormap(self, colortable, vbins):
		from tvtk.util.ctf import ColorTransferFunction

		ctf = ColorTransferFunction()

		for i in range(len(colortable)):
			ctf.add_rgb_point(vbins[i], colortable[i][0] / 256., colortable[i][1] / 256., colortable[i][2] / 256.)

		return ctf

	def _create_grid_data(self, qty, family=None, width=None, snap_slice=None, recalc=False):
		ss = self.sim
		data_name = 'all_'+qty
		if family in ['star','stars']:
			data_name = 'star_'+qty
			ss = self.sim.s
		if family=='gas':
			data_name = 'gas_' + qty
			ss = self.sim.g
		if family in ['dm','dark']:
			data_name = 'dm_' + qty
			ss = self.sim.dm

		if width:
			data_name = data_name+'_'+str(width)

		if data_name in self._loaded_data.keys() and not recalc:
			print("using previously calculated data grid")
			grid_data = self._loaded_data[data_name]
		else:
			grid_data = sph.to_3d_grid(ss, qty=qty, nx=self.resolution,snap_slice=snap_slice,
		                           x2=None if width is None else width / 2)
			self._loaded_data[data_name] = grid_data
		return grid_data

	def set_starsize(self, size):
		newsize = float(size) #make sure the input actually can be converted
		self._starsize = newsize

	def render(self, qty, family=None, width=None, vmin=None, vmax=None, dynamic_range=4,
	           log=True, color=None, colortable=None, create_figure=True, snap_slice=None,
	           recalc=False, bins=None, clear=True, cut='low', max_opacity=None):

		import mayavi
		from mayavi import mlab
		from tvtk.util.ctf import PiecewiseFunction, ColorTransferFunction
		import palettable

		if type(qty) != str:
			raise ValueError("qty must be a strong, e.g. 'rho', 'temp'")

		if family is not None:
			if family not in ['gas','star','stars','dm', 'dark']:
				raise ValueError("family must be one of these strings: 'gas','star','dm'")

		if family in ['star','stars']:
			if self._starsize:
				smf = filt.HighPass('smooth', str(self._starsize) + ' kpc')
				self.sim.s[smf]['smooth'] = array.SimArray(self._starsize, 'kpc', sim=self.sim)
			if snap_slice is None:
				#default stars to not include BHs
				snap_slice = filt.HighPass('tform',0)

		grid_data = self._create_grid_data(qty, family=family, width=width, snap_slice=snap_slice, recalc=recalc)

		if create_figure:
			fig = mlab.figure(size=(500, 500), bgcolor=(0, 0, 0))
		if clear:
			mlab.clf()

		if family=='star' and (qty=='tform' or qty=='age'):
			sim_time = self.sim.properties['time'].in_units('Gyr')
			if bins is None and vmin is None and vmax is None:
				print("using default bins for ", qty)
				if qty=='tform':
					bins = sim_time - np.array([0.0, 1.0, 3.0, 4.0, 6.0, 10.0, 14.0])
					bins = bins[::-1]
				if qty == 'age':
					bins = np.array([0.0, 1.0, 3.0, 4.0, 6.0, 10.0, 14.0])
				if log:
					bins = np.log10(bins)
		if bins is not None:
			vmin = bins.min()
			vmax = bins.max()

		if log:
			grid_data = np.log10(grid_data)
			if vmax is None:
				vmax = grid_data.max()
			if vmin is None:
				vmin = grid_data.max() - dynamic_range
		else:
			if vmin is None:
				vmin = np.min(grid_data)
			if vmax is None:
				vmax = np.max(grid_data)

		#grid_data[grid_data < vmin] = vmin
		#grid_data[grid_data > vmax] = vmax
		otf = PiecewiseFunction()
		if not max_opacity:
			max_opacity = 1
		if cut=='high':
			otf.add_point(vmax,0.0)
			otf.add_point(vmin,max_opacity)
		if cut=='low':
			otf.add_point(vmin,0)
			otf.add_point(vmax,max_opacity)
		if cut=='none':
			otf.add_point((vmax-vmin)/2.,max_opacity)
			otf.add_point(vmin,max_opacity*0.5)
			otf.add_point(vmax,max_opacity*0.5)


		sf = mayavi.tools.pipeline.scalar_field(grid_data)
		V = mlab.pipeline.volume(sf, color=color, vmin=vmin, vmax=vmax)

		V.trait_get('volume_mapper')['volume_mapper'].blend_mode = 'maximum_intensity'

		if color is None:
			if colortable is None: #default colormap is cubehelix unless looking at stellar age
				if family=='star' and qty in ['tform','age']:
					colortable = np.array(palettable.lightbartlein.diverging.BlueOrange10_6.colors)
					if qty == 'tform': colortable = colortable[::-1] #reverse blue/orange for tform values
				else:
					colortable = np.array(palettable.cubehelix.cubehelix1_16.colors)
			vbins = np.arange(vmin, vmax, (vmax - vmin) / len(colortable))
			ctf = self._create_colormap(colortable,vbins)
			V._volume_property.set_color(ctf)
			V._ctf = ctf
			V.update_ctf = True

		V._otf = otf
		V._volume_property.set_scalar_opacity(otf)

		return V

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

	def star_tform(self,min_age=None, max_age=None, log=False, starsize=0.25,
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
			self.sim.s[smf]['smooth'] = array.SimArray(starsize, 'kpc', sim=self.sim)

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
		otf.add_point(10.0, 1.0)

		sf = mayavi.tools.pipeline.scalar_field(grid_data)
		V = mlab.pipeline.volume(sf, color=color, vmin=vmin, vmax=vmax)

		V.trait_get('volume_mapper')['volume_mapper'].blend_mode = 'maximum_intensity'

		if color is None:
			if colortable is None: #default colormap is BlueOrange10
				colortable = np.array(palettable.lightbartlein.diverging.BlueOrange10_6.colors)[::-1]
			if age_bins is None:
				age_bins = np.array([1.0, 3.0, 4.0, 0.6,10.0,14.0])
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

