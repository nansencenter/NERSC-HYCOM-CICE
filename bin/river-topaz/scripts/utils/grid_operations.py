#!/usr/bin/env python
# vim: ts=4 expandtab sw=4 sts=4
import logging
from typing import List, Tuple

import numpy as np
from scipy import spatial
from xarray import DataArray, Dataset, broadcast, open_dataset, zeros_like

log = logging.getLogger()





def transform_coordinates(coords):
    """Transform coordinates from geodetic to cartesian

    Keyword arguments:
    coords - a set of lan/lon coordinates (e.g. a tuple or
             an array of tuples)
    """
    # WGS 84 reference coordinate system parameters
    A = 6378.137  # major axis [km]
    E2 = 6.69437999014e-3  # eccentricity squared

    coords = np.asarray(coords).astype("f4")

    # is coords a tuple? Convert it to an one-element array of tuples
    if coords.ndim == 1:
        coords = np.array([coords])

    # convert to radiants
    lat_rad = np.radians(coords[:, 0])
    lon_rad = np.radians(coords[:, 1])

    # convert to cartesian coordinates
    r_n = A / (np.sqrt(1 - E2 * (np.sin(lat_rad) ** 2)))
    x = r_n * np.cos(lat_rad) * np.cos(lon_rad)
    y = r_n * np.cos(lat_rad) * np.sin(lon_rad)
    z = r_n * (1 - E2) * np.sin(lat_rad)

    return np.column_stack((x, y, z))




def refine_grid(
    lon_in: DataArray,
    lat_in: DataArray,
    lon_out: DataArray,
    lat_out: DataArray,
    input_shape: Tuple,
    output_shape: Tuple,
):
    """
    Refine grid with kdtree
    Computed by identifying nearest neighbour (in input grid)
    of all the output grid points.
    """

    # stack coordinates
    try:
        coords_in = np.column_stack((lat_in.values.ravel(), lon_in.values.ravel()))
        coords_out = np.column_stack((lat_out.values.ravel(), lon_out.values.ravel()))
    except:
        coords_in = np.column_stack((lat_in.ravel(), lon_in.ravel()))
        coords_out = np.column_stack((lat_out.ravel(), lon_out.ravel()))
    # log.debug(f"{coords_in=}")
    # log.debug(f"{coords_out=}")
    # log.debug(f"{coords_in.shape=}")
    # log.debug(f"{coords_out.shape=}")

    coords_in = np.where(np.isnan(coords_in), 0, coords_in)

    # construct KD-tree
    ground_pixel_tree = spatial.cKDTree(transform_coordinates(coords_in))
    # ground_pixel_tree = spatial.cKDTree(coords_in)

    # get distances and indices
    log.debug(f"{input_shape=}")
    log.debug(f"{output_shape=}")
    distances, indices = get_closests_points_kdtree(
        ground_pixel_tree, coords_out, input_shape
    )
    # reshape if user specified to
    distances, indices = reshape_matrices_(distances, indices, output_shape)
    return distances, indices


def reshape_matrices_(distances, indices, output_shape):
    distances = distances.reshape(output_shape)
    log.debug(len(indices))
    indices = list(indices)
    for dim in range(len(indices)):
        log.debug(f"before reshape {indices[dim].shape=}")
        indices[dim] = indices[dim].reshape(output_shape)
        log.debug(f"after reshape {indices[dim].shape=}")
    indices = tuple(indices)
    log.debug(f"{indices=}")
    return distances, indices



def get_closests_points_kdtree(tree, lat_lon_targets, input_shape):
    distances, indices = tree.query(transform_coordinates(lat_lon_targets))
    log.debug(f"{indices=}")
    indices = np.unravel_index(indices, input_shape)
    log.debug(f"{indices=}")
    return distances, indices

