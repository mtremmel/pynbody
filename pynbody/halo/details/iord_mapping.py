from __future__ import annotations

import abc
import warnings

import numpy as np

from pynbody.util import binary_search, is_sorted


class IordToOffsetFromFile(abc.ABC):
    @abc.abstractmethod
    def map_preserving_order(self, i: np.ndarray | int) -> np.ndarray | int:
        """Given an array of iord values, return the corresponding fpos values.

        This can be used ONLY for simulations where the particle order within the file preserves the order of Iord."""
        pass

class IordToOffsetFromFileAscending(IordToOffsetFromFile):
    def __init__(self,iord_info):
        """
        The input is a list of dictionary items denoting blocks of iord values with the following format
        blockN = {'ind_start': starting indices of block within a SimSnap object,
            'ind_end':ending indice of block,
            'first_iord':first iord,
            'last_iord':last iord,
            'missing_iords':a list of missing iords for the block(e.g. corresponding to deleted particles)}
        """
        self._iord_info = iord_info
        for i in range(len(self._iord_info)): #make sure missing iords are also in ascending order
            osort = np.argsort(self._iord_info[i]['missing_iords'])
            self._iord_info[i]['missing_iords'] = self._iord_info[i]['missing_iords'][osort]

    def do_map(self,i):
        offset_array = np.ones(len(i)).astype(np.int64)*-1 #initialize an array of index values
        cnt = 1
        for block in self._iord_info:
            mask_block = (i>=block['first_iord'])&(i<=block['last_iord']) #get all iords within the block
            block_ind = np.where(mask_block)[0]
            offset_array[block_ind] = i[block_ind] - block['first_iord'] + block['ind_start'] #produce an initial offset
            osort_block = np.argsort(i[block_ind])
            ind_miss = np.searchsorted(i[block_ind],block['missing_iords'],sorter=osort_block)
            uind_miss, cnt_miss = np.unique(ind_miss,return_counts=True)
            for ind,cnt in zip(uind_miss,cnt_miss):
                offset_array[block_ind[osort_block[ind:]]]-=cnt
            if offset_array[mask_block][osort_block].min()<block['ind_start']:
                warnings.warn("indices found for halo particles go below the provided starting indices for block", cnt, "provided in file")
            if offset_array[mask_block][osort_block].max()>block['ind_end']:
                warnings.warn("indices found for halo particles go above the provided starting indices for block", cnt, "provided in file")
            cnt += 1
        good = offset_array>=0 #we only care about those elements that were able to be mapped to indices.
        #if the iord information file does not encompass all iords in the halo, those missing will be ignored.
        if len(offset_array[good])<len(offset_array) and len(offset_array[good])>0:
            warnings.warn("Not all iords in the halo were mapped according to _iordInfo.npy file. Ignoring those that failed to receive a mapping",
                              RuntimeWarning)
        if len(offset_array[good])==0:
            raise RuntimeError("Halo has no particles that were able to be mapped according to their iords")
        return offset_array[good]
    
    def map_preserving_order(self,i):
        offset_array = self.do_map(i)
        return offset_array
    
class IordToOffsetFromFileDescending(IordToOffsetFromFileAscending): #same as Ascending, but reversed.
    def __init__(self,iord_info):
        self._iord_info = iord_info
        for i in range(len(self._iord_info)): #make sure missing iords are also in ascending order
            self._iord_info[i]['missing_iords'] *= -1 #if we multiply the descending iords by -1, then the order becomes ascending
            self._iord_info[i]['last_iord'] *= -1
            self._iord_info[i]['first_iord'] *= -1
            osort = np.argsort(self._iord_info[i]['missing_iords'])
            self._iord_info[i]['missing_iords'] = self._iord_info[i]['missing_iords'][osort]
    
    def map_preserving_order(self,i):
        i *= -1 #again, we multiply the input particle iords by -1 to turn them into ascending order
        offset_array = self.do_map(i)
        return offset_array


class IordToOffset(abc.ABC):
    @abc.abstractmethod
    def map_ignoring_order(self, i: np.ndarray | int) -> np.ndarray | int:
        """Given an array of iord values, return the corresponding fpos values.

        Warning: The returned values are not guaranteed to be in the same order as the input iord array."""
        pass


class IordToOffsetDense(IordToOffset):
    def __init__(self, iord_array, max_iord=None):
        if max_iord is None:
            max_iord = int(iord_array.max())
        self._iord_to_offset = np.empty(max_iord + 1, dtype=np.int64)
        self._iord_to_offset.fill(-1)
        self._iord_to_offset[iord_array] = np.arange(len(iord_array), dtype=np.int64)

    def map_ignoring_order(self, i):
        return self._iord_to_offset[i]


class IordToOffsetSparse(IordToOffset):
    """Class for efficiently mapping from iords to offsets in the iord array, even if iord values are large.

    WARNING: if a query is made with iords that are not themselves in ascending order, a sort takes place
    ahead of the query and therefore the set returned is correct but the ordering is not preserved."""
    def __init__(self, iord_array):
        self._iord = iord_array
        self._iord_argsort = np.argsort(iord_array)

    def map_ignoring_order(self, iord_values: np.ndarray | int) -> np.ndarray | int:
        if not hasattr(iord_values, "__len__"):
            iord_values = np.array([iord_values])
            singleton = True
        else:
            iord_values = np.asarray(iord_values)
            singleton = False

            if is_sorted(iord_values) != 1:
                iord_values = np.sort(iord_values)

        result = binary_search(np.asarray(iord_values), self._iord, self._iord_argsort)

        if singleton:
            return result[0]
        else:
            return result


class IordOffsetModifier(IordToOffset):
    """A wrapper around an IordToOffset which adds a constant offset to the result of the underlying mapping.

    Useful if the iord values e.g. are only available for a single family; then the fpos_offset will correspond
    to the first index of that family in the pynbody snapshot.
    """

    def __init__(self, iord_to_offset: IordToOffset, fpos_offset: int):
        self._underlying = iord_to_offset
        self._fpos_offset = fpos_offset

    def map_ignoring_order(self, i: np.ndarray | int) -> np.ndarray | int:
        result = self._underlying.map_ignoring_order(i)
        result += self._fpos_offset
        return result

def make_iord_to_offset_mapper_from_file(iord_info) -> IordToOffsetFromFile:
    """
    Given iords that are preserved in order within the output file, produce a mapping of iord to index
    within a pynbody SimSnap.

    Input is a list of dictionary objects for each block of particles such that
    iord_info = [block1, block2,...]
    and blockN = {'ind_start': starting indices of block in SimSnap,
            'ind_end':ending indice of block,
            'first_iord':first iord,
            'last_iord':last iord,
            'missing_iords':a list of missing iords for the block(e.g. corresponding to deleted particles)}
    """
    if iord_info[0]['first_iord']<iord_info[0]['last_iord']:
        return IordToOffsetFromFileAscending(iord_info)
    else:
        return IordToOffsetFromFileDescending(iord_info)


def make_iord_to_offset_mapper(iord: np.ndarray) -> IordToOffset:
    """Given an array of unique integers, iord, make an object which maps from an iord value to offset in the array.

    i.e. given an iord array and a subset of values my_iord_values,

     make_iord_to_offset_mapper(iord).map_ignoring_order(my_iord_values)

    returns the indexes of my_iord_values in the iord array.
    """

    min_iord = int(iord.min())
    max_iord = int(iord.max())

    if (min_iord >= 0) and (max_iord < 2 * len(iord)):
        # maximum iord is not very big, just do a direct in-memory mapping for speed
        return IordToOffsetDense(iord, max_iord)
    else:
        # maximum iord is large, so we'll use util.binary_search to save memory at the cost of speed
        return IordToOffsetSparse(iord)
