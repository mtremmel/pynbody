import pylab as p
import matplotlib
import numpy as np
from .. import sph, config
from .. import units as _units
from .. import filt, array

default_colors = {  'red': (1., 0., 0.),
					'blue': (0., 0., 1.),
					'yellow':  (1.,1., 0.),
					'green': (0.,1.,0.)
					}

def _create_colormap(colortable, vbins):
	from tvtk.util.ctf import ColorTransferFunction

	ctf = ColorTransferFunction()

	for i in range(len(colortable)):
		ctf.add_rgb_point(vbins[i], colortable[i][0] / 256., colortable[i][1] / 256., colortable[i][2] / 256.)

	return ctf

def _set_bins(vmin, vmax, nbins):
	binsize  = (vmax-vmin)/nbins
	bins = np.arange(vmin, vmax+binsize, binsize)
	return bins

def _get_opacities(vmin, vmax, max_opacity, cut):
	from tvtk.util.ctf import PiecewiseFunction
	otf = PiecewiseFunction()
	if not max_opacity:
		max_opacity = 1
	if cut == 'high':
		otf.add_point(vmax, 0.0)
		otf.add_point(vmin, max_opacity)
	if cut == 'low':
		otf.add_point(vmin, 0)
		otf.add_point(vmax, max_opacity)
	if cut == 'none':
		otf.add_point((vmax - vmin) / 2., max_opacity)
		otf.add_point(vmin, max_opacity * 0.5)
		otf.add_point(vmax, max_opacity * 0.5)
	if cut == 'middle':
		otf.add_point((vmax - vmin) / 2., 0)
		otf.add_point(vmin, max_opacity)
		otf.add_point(vmax, max_opacity)
	if cut == 'both':
		otf.add_point((vmax-vmin)/2.,max_opacity)
		otf.add_point(vmin, 0)
		otf.add_point(vmax, 0)
	return otf

class RenderVolume(object):
	def __init__(self, sim, resolution=500, load_file=None):
		if sim is None and load_file is None:
			raise RuntimeError("Cannot create a volume render with neither a simulation nor previously saved grids!")
		if sim:
			sim.physical_units() # make sure things are in physical units!
		self.sim = sim
		self.resolution = resolution
		self._loaded_data = {} #store calculated grids for faster plotting
		self._starsize=0.25 #default starsize
		self._load_file = load_file
		if load_file:
			self.load(load_file)

	def clear_cache(self):
		import gc
		del self._loaded_data
		gc.collect()


	def load(self, filename, overwrite=True):
		import pickle
		f = open(filename,'rb')
		saved_data = pickle.load(f)
		if type(saved_data)!=dict:
			raise ValueError("Error Loading datafile "+filename+" Expecting a pickle file with a dictionary of 3d grids!")
		for key in saved_data.keys():
			if key not in self._loaded_data.keys() or overwrite is False:
				print("loading in ", key, " from ", filename)
				self._loaded_data[key] = saved_data[key]

	def save(self, filename):
		import pickle
		if len(self._loaded_data.keys())==0:
			raise RuntimeError("no data currently loaded!")
		f = open(filename,'wb')
		pickle.dump(self._loaded_data,f)
		f.close()

	# def _create_binned_grid_data(self, ss, qty, binned_qty, bins, width, log):
	# 	output = []
	# 	for i in range(len(bins)-1):
	# 		if log:
	# 			ss_part = ss[filt.BandPass(binned_qty, 10**bins[i], 10**bins[i+1])]
	# 		else:
	# 			ss_part = ss[filt.BandPass(binned_qty, bins[i], bins[i + 1])]
	# 		grid_part = sph.to_3d_grid(ss_part, qty=qty, nx=self.resolution,
	# 		                           x2=None if width is None else width / 2)
	# 		output.append(grid_part)
	# 	return output

	# def _create_grid_data_old(self, qty, nbins, family=None, width=None, recalc=False, weight=None,
	#                       vmin=None, vmax=None, log=True, save=True, dynamic_range=4):
	# 	ss = self.sim
	# 	data_name = 'all_'+qty
	# 	if family in ['star','stars']:
	# 		data_name = 'star_'+qty
	# 		ss = self.sim.s[filt.HighPass('tform',0)] #remove black holes
	# 	if family=='gas':
	# 		data_name = 'gas_' + qty
	# 		ss = self.sim.g
	# 	if family in ['dm','dark']:
	# 		data_name = 'dm_' + qty
	# 		ss = self.sim.dm
	#
	# 	if width:
	# 		data_name = data_name+'_'+str(width)
	#
	# 	if weight is not None:
	# 		data_name = data_name+"_"+weight
	#
	# 	bins = self._set_bins(ss, qty, vmin, vmax, dynamic_range, log, nbins)
	#
	# 	if data_name in self._loaded_data.keys() and weight:
	# 		if nbins != len(self._loaded_data[data_name]):
	# 			if not recalc:
	# 				print("Warning! Provided colortable requires different binning... recalculating weighted 3d render")
	# 				recalc = True
	#
	# 	if data_name in self._loaded_data.keys() and not recalc:
	# 		print("using previously calculated data grid")
	# 		grid_data = np.copy(self._loaded_data[data_name])
	#
	# 	else:
	# 		if weight:
	# 			grid_data = self._create_binned_grid_data(ss, weight, qty, bins, width, log)
	# 		else:
	# 			grid_data = sph.to_3d_grid(ss, qty=qty, nx=self.resolution,
	# 		                           x2=None if width is None else width / 2)
	# 		if save:
	# 			self._loaded_data[data_name] = np.copy(grid_data)
	# 	return grid_data, bins

	def _create_grid_data(self, qty, filter=None, family=None, width=None, recalc=False, save=True):
		ss = self.sim
		data_name = 'all_'+qty
		if family in ['star','stars']:
			data_name = 'star_'+qty
			if ss is not None:
				ss = self.sim.s[filt.HighPass('tform',0)] #remove black holes
		if family=='gas':
			data_name = 'gas_' + qty
			if ss is not None:
				ss = self.sim.g
		if family in ['dm','dark']:
			data_name = 'dm_' + qty
			if ss is not None:
				ss = self.sim.dm

		if width:
			data_name = data_name+'_'+str(width)

		if filter:
			data_name += '_'+filter._descriptor
			if hasattr(filter,'_min'):
				data_name += '_'+str(filter._min)
			if hasattr(filter,'_max'):
				data_name += '_'+str(filter._max)
			ss = ss[filter]

		#bins = _set_bins(ss, qty, vmin, vmax, dynamic_range, log, nbins)

		# if data_name in self._loaded_data.keys() and weight:
		# 	if nbins != len(self._loaded_data[data_name]):
		# 		if not recalc:
		# 			print("Warning! Provided colortable requires different binning... recalculating weighted 3d render")
		# 			recalc = True

		if data_name in self._loaded_data.keys() and not recalc:
			print("using previously calculated data grid")
			grid_data = np.copy(self._loaded_data[data_name])

		else:
			if ss is None:
				raise ValueError("simulation has not been provided (None) but requested data does not already exist")
			if family in ['star', 'stars']:
				if self._starsize:
					smf = filt.HighPass('smooth', str(self._starsize) + ' kpc')
					ss[smf]['smooth'] = array.SimArray(self._starsize, 'kpc', sim=self.sim)

			grid_data = sph.to_3d_grid(ss, qty=qty, nx=self.resolution,
			                           x2=None if width is None else width / 2)
			if save:
				self._loaded_data[data_name] = np.copy(grid_data)
		return grid_data

	def set_starsize(self, size):
		newsize = float(size) #make sure the input actually can be converted
		self._starsize = newsize

	def render(self, qty, filter=None, family=None, width=None, vmin=None, vmax=None, dynamic_range=4,
	           log=True, color=None, colortable=None, create_figure=True,
	           recalc=False, clear=True, cut='low', max_opacity=None, weight=None):

		import mayavi
		from mayavi import mlab
		import palettable

		nbins = None

		if colortable is None and color is None:
			#default to something reasonable
			colortable = np.array(palettable.matplotlib.Viridis_16.colors)
			if qty in ['tform', 'age']:
				colortable = np.array(palettable.lightbartlein.diverging.BlueOrange10_6.colors)
			if qty == 'temp':
				colortable = np.array(palettable.lightbartlein.diverging.BlueDarkRed18_16.colors)
			if qty == 'rho':
				colortable = np.array(palettable.cubehelix.cubehelix1_16.colors)
		if color:
			if type(color)==str:
				if color not in default_colors.keys():
					raise ValueError("Provided color name does not match default colors", default_colors.keys())
				color = default_colors[color]

		if colortable is not None:
			nbins = len(colortable)

		if type(qty) != str:
			raise ValueError("qty must be a string, e.g. 'rho', 'temp'")

		if family is not None:
			if family not in ['gas','star','stars','dm', 'dark']:
				raise ValueError("family must be one of these strings: 'gas','star','dm'")

		grid_data = self._create_grid_data(qty, filter=filter, family=family, width=width, recalc=recalc)

		if create_figure:
			fig = mlab.figure(size=(500, 500), bgcolor=(0, 0, 0))
		if clear:
			mlab.clf()

		if dynamic_range:
			grid_data[(grid_data < np.max(grid_data) / 10 ** dynamic_range)] = np.max(
				grid_data) / 10 ** dynamic_range
		else:
			if log:
				if vmin:
					grid_data[(grid_data<10**vmin)] = 10**vmin
				if vmax:
					grid_data[(grid_data>10**vmax)] = 10**vmax

				grid_data = np.log10(grid_data)
			else:
				if vmax:
					grid_data[(grid_data >vmax)] = vmax
				if vmin:
					grid_data[(grid_data<vmin)] = vmin

		global_max = np.max(grid_data)
		global_min = np.min(grid_data)

		otf = _get_opacities(global_min, global_max, max_opacity, cut)
		sf = mayavi.tools.pipeline.scalar_field(grid_data)
		V = mlab.pipeline.volume(sf, color=color, vmin=global_min, vmax=global_max)
		if colortable is not None:
			bins = _set_bins(global_min, global_max, nbins)
			print("bins created:", bins)
			ctf = _create_colormap(colortable,bins)
			V._volume_property.set_color(ctf)
			V._ctf = ctf
			V.update_ctf = True
		V.trait_get('volume_mapper')['volume_mapper'].blend_mode = 'maximum_intensity'
		V._otf = otf
		V._volume_property.set_scalar_opacity(otf)

		# else:
		# 	V = []
		# 	for i in range(len(grid_data)):
		# 		if np.abs(grid_data[i].max()) == np.inf: #skip any bins with zero weight
		# 			continue
		# 		otf = _get_opacities(np.max(grid_data[i])-dynamic_range_weights, np.max(grid_data[i]), max_opacity, 'low')
		# 		sf = mayavi.tools.pipeline.scalar_field(grid_data[i])
		# 		V_part = mlab.pipeline.volume(sf, color=tuple(colortable[i]/255), vmin=np.max(grid_data[i])-dynamic_range_weights, vmax=np.max(grid_data[i]))
		# 		V_part.trait_get('volume_mapper')['volume_mapper'].blend_mode = 'maximum_intensity'
		# 		V_part._otf = otf
		# 		V_part._volume_property.set_scalar_opacity(otf)
		# 		V.append(V_part)


		return V

	# def density(self, family=None, vmin=None, vmax=None, dynamic_range=4,
	#             log=True, color=None, colortable=None, create_figure=True):
	# 	import mayavi
	# 	from mayavi import mlab
	# 	from tvtk.util.ctf import PiecewiseFunction, ColorTransferFunction
	# 	import palettable
	#
	# 	if create_figure:
	# 		fig = mlab.figure(size=(500, 500), bgcolor=(0, 0, 0))
	#
	# 	data_name = 'all_den'
	# 	ss = self.sim
	# 	if family == 'gas':
	# 		data_name = 'gas_den'
	# 		ss = self.sim.g
	# 	if family == 'dm':
	# 		data_name = 'dm_den'
	# 		ss = self.sim.dm
	# 	if family == 'star':
	# 		data_name = 'star_den'
	# 		ss = self.sim.s
	#
	# 	if data_name in self._loaded_data.keys():
	# 		print("using previously calculated data grid")
	# 		grid_data = self._loaded_data[data_name]
	# 	else:
	# 		grid_data = sph.to_3d_grid(ss, qty='rho', nx=self.resolution,
	# 	                           x2=None if self.width is None else self.width / 2)
	# 		self._loaded_data[data_name] = grid_data
	#
	# 	if log:
	# 		grid_data = np.log10(grid_data)
	# 		if vmin is None:
	# 			vmin = grid_data.max() - dynamic_range
	# 		if vmax is None:
	# 			vmax = grid_data.max()
	#
	# 	else:
	# 		if vmin is None:
	# 			vmin = np.min(grid_data)
	# 		if vmax is None:
	# 			vmax = np.max(grid_data)
	#
	# 	grid_data[grid_data < vmin] = vmin
	# 	grid_data[grid_data > vmax] = vmax
	#
	# 	otf = PiecewiseFunction()
	# 	otf.add_point(vmin, 0.0)
	# 	otf.add_point(vmax, 1.0)
	#
	# 	sf = mayavi.tools.pipeline.scalar_field(grid_data)
	# 	V = mlab.pipeline.volume(sf, color=color, vmin=vmin, vmax=vmax)
	#
	# 	V.trait_get('volume_mapper')['volume_mapper'].blend_mode = 'maximum_intensity'
	#
	# 	if color is None:
	# 		if colortable is None: #default colormap is cubehelix
	# 			colortable = np.array(palettable.cubehelix.cubehelix1_16.colors)
	# 		vbins = np.arange(vmin, vmax, (vmax - vmin) / len(colortable))
	# 		ctf = self._create_colormap(colortable,vbins)
	# 		V._volume_property.set_color(ctf)
	# 		V._ctf = ctf
	# 		V.update_ctf = True
	#
	# 	V._otf = otf
	# 	V._volume_property.set_scalar_opacity(otf)
	#
	# 	return V
	#
	# def gas_density(self,**kwargs):
	# 	return self.density(family='gas', **kwargs)
	#
	# def dm_density(self,**kwargs):
	# 	return self.density(family='dm', **kwargs)
	#
	# def star_density(self, **kwargs):
	# 	return self.density(family='star', **kwargs)
	#
	# def star_tform(self,min_age=None, max_age=None, log=False, starsize=0.25,
	#             color=None, colortable=None, create_figure=True, age_bins=None):
	#
	# 	import mayavi
	# 	from mayavi import mlab
	# 	from tvtk.util.ctf import PiecewiseFunction, ColorTransferFunction
	# 	import palettable
	#
	# 	data_name = 'star_tform'
	#
	# 	if create_figure:
	# 		fig = mlab.figure(size=(500, 500), bgcolor=(0, 0, 0))
	#
	#
	# 	if starsize is not None:
	# 		smf = filt.HighPass('smooth', str(starsize) + ' kpc')
	# 		self.sim.s[smf]['smooth'] = array.SimArray(starsize, 'kpc', sim=self.sim)
	#
	# 	if data_name in self._loaded_data.keys():
	# 		print("using previously calculated data grid")
	# 		grid_data = self._loaded_data[data_name]
	# 	else:
	# 		grid_data = sph.to_3d_grid(self.sim.s, qty='tform', nx=self.resolution, snap_slice=filt.HighPass('tform',0),
	# 	                           x2=None if self.width is None else self.width / 2)
	# 		self._loaded_data[data_name] = grid_data
	#
	# 	grid_data = grid_data.in_units('Gyr')
	#
	# 	sim_time = self.sim.properties['time'].in_units('Gyr')
	#
	#
	# 	if max_age is None:
	# 		vmin=0
	# 	else:
	# 		vmin = sim_time-max_age
	#
	# 	if min_age is None:
	# 		vmax = sim_time
	# 	else:
	# 		vmax = sim_time - min_age
	#
	# 	if log is True:
	# 		if vmin <= 0:
	# 			vmin = 0.001
	# 		grid_data = np.log10(grid_data)
	# 		vmax = np.log10(vmax)
	#
	# 	grid_data[(grid_data>vmax)] = vmax
	# 	grid_data[(grid_data<vmin)] = vmin
	#
	# 	otf = PiecewiseFunction()
	# 	otf.add_point(vmin, 0.0)
	# 	otf.add_point(10.0, 1.0)
	#
	# 	sf = mayavi.tools.pipeline.scalar_field(grid_data)
	# 	V = mlab.pipeline.volume(sf, color=color, vmin=vmin, vmax=vmax)
	#
	# 	V.trait_get('volume_mapper')['volume_mapper'].blend_mode = 'maximum_intensity'
	#
	# 	if color is None:
	# 		if colortable is None: #default colormap is BlueOrange10
	# 			colortable = np.array(palettable.lightbartlein.diverging.BlueOrange10_6.colors)[::-1]
	# 		if age_bins is None:
	# 			age_bins = np.array([1.0, 3.0, 4.0, 0.6,10.0,14.0])
	# 		else:
	# 			age_bins = np.array(age_bins) #make sure age_bins is an array
	# 		vbins = sim_time - age_bins
	# 		vbins = vbins[::-1]
	# 		if log is True:
	# 			vbins = np.log10(vbins)
	# 		ctf = self._create_colormap(colortable,vbins)
	# 		V._volume_property.set_color(ctf)
	# 		V._ctf = ctf
	# 		V.update_ctf = True
	#
	# 	V._otf = otf
	# 	V._volume_property.set_scalar_opacity(otf)
	#
	# 	return V

