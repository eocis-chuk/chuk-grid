#!/usr/bin/env python
# -*- coding: utf-8 -*-

#     create prototype EOCIS high resolution lat-lon grids for the British Isles
#
#     Copyright (C) 2023  EOCIS and National Centre for Earth Observation (NCEO)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

import xarray as xr
import numpy as np
import pyproj
import datetime


def check_bnds(bnds,j=3,i=6):
    """Perform a sanity check on a 2d bounds matrix according to follwing psuedo-code

    For 0 < j < n and 0 < i < m,
      If cells (j,i) and (j,i+1) are contiguous, then
        bnd(j,i,1)=bnd(j,i+1,0)
        bnd(j,i,2)=bnd(j,i+1,3)
      If cells (j,i) and (j+1,i) are contiguous, then
        bnd(j,i,3)=bnd(j+1,i,0) and bnd(j,i,2)=bnd(j+1,i,1)

    (from http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries)"""
    if bnds[j, i, 1] != bnds[j, i + 1, 0]:
        print(bnds[j, i, 1],bnds[j, i + 1, 0])
        raise Exception("check_bounds A")
    if bnds[j, i, 2] != bnds[j, i + 1, 3]:
        raise Exception("check_bounds B")
    if bnds[j, i, 3] != bnds[j + 1, i, 0]:
        raise Exception("check_bounds C")
    if bnds[j, i, 2] != bnds[j + 1, i, 1]:
        raise Exception("check_bounds D")


def create_grid(min_n, min_e, max_n, max_e, spacing_m, version, include_bounds):
    """
    Create an xarray Dataset describing a grid on EPSG:27700

    :param min_n: minimum BNG northing value (m) of area extent
    :param min_e: minimum BNG easting value (m) of area extent
    :param max_n: maximum BNG northing value (m) of area extent
    :param max_e: maximum BNG easting value (m) of area extent
    :param spacing_m: grid spacing (m)
    :param version: version string to add as an attribute
    :param include_bounds: whether to output lat_bnds, lon_bnds, x_bnds, y_bnds

    :return: xarray dataset consisting of the following arrays

    latitude    (nj,ni) array of WGS84 latitude values calculated via pyproj
    longitude   (nj,ni) array of WGS84 longitude values calculated via pyproj
    northings   (nj) array contiaining BNG northings
    eastings     (ni) array containing BNG eastings
    """

    # Northings from north to south (decreasing)
    n_val = np.linspace(max_n, min_n, int((max_n-min_n)/spacing_m)+1, dtype=np.int32)
    # Eastings from west to east (increasing)
    e_val = np.linspace(min_e, max_e, int((max_e - min_e) / spacing_m) + 1, dtype=np.int32)

    shape = len(n_val), len(e_val)

    northings = np.broadcast_to(n_val[None].T, shape)
    eastings = np.broadcast_to(e_val, shape)

    transformer = pyproj.Transformer.from_crs(27700, 4326)
    lats, lons = transformer.transform(eastings, northings)

    if include_bounds:
        n_bnds = np.zeros((len(n_val), 2), dtype=np.int32)
        n_bnds[:, 0] = n_val + spacing_m / 2
        n_bnds[:, 1] = n_val - spacing_m / 2

        e_bnds = np.zeros((len(e_val), 2), dtype=np.int32)
        e_bnds[:, 0] = e_val - spacing_m / 2
        e_bnds[:, 1] = e_val + spacing_m / 2
        n_northings = np.broadcast_to(n_bnds[:, 0][None].T, shape)
        s_northings = np.broadcast_to(n_bnds[:, 1][None].T, shape)
        w_eastings = np.broadcast_to(e_bnds[:, 0], shape)
        e_eastings = np.broadcast_to(e_bnds[:, 1], shape)

        # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
        latlon_bnds_shape = len(n_val), len(e_val), 4 # last dimension orders corners as NW, NE, SE, SW

        lat_bnds = np.zeros(latlon_bnds_shape)
        lon_bnds = np.zeros(latlon_bnds_shape)

        for (idx,e,n) in [(0,w_eastings, n_northings),(1,e_eastings, n_northings),(2,e_eastings,s_northings),(3,w_eastings,s_northings)]:
            print(f"Creating bounds for index {idx}")
            lats1, lons1 = transformer.transform(e, n)
            lat_bnds[..., idx] = lats1
            lon_bnds[..., idx] = lons1

        # sanity check some bounds according to the CF doc
        for bnds in [lat_bnds, lon_bnds]:
            check_bnds(bnds)

    print("Lat/Lon Bounding Box: Lat: %.4f to %.4f, Lon: %.4f to %.4f" % (np.min(lats),np.max(lats),np.min(lons),np.max(lons)))

    spec = {
        'lat': (['y', 'x'], lats),
        'lon': (['y', 'x'], lons),
        'crsOSGB': ([], np.int32(0), pyproj.Proj(27700).crs.to_cf())
    }

    if include_bounds:
        spec.update({
            'x_bnds': (['x', 'bnds'], e_bnds),
            'y_bnds': (['y', 'bnds'], n_bnds),
            'lat_bnds': (['y', 'x', 'latlon_bnds'], lat_bnds),
            'lon_bnds': (['y', 'x', 'latlon_bnds'], lon_bnds)
        })

    ds = xr.Dataset(spec,
        coords={
            'x': e_val,
            'y': n_val,
        })

    ds.lat.attrs.update(standard_name='latitude', units='degrees_north',bounds='lat_bnds')
    ds.lon.attrs.update(standard_name='longitude', units='degrees_east',bounds='lon_bnds')
    ds.x.attrs.update(long_name='Easting', standard_name='projection_x_coordinate', units='m',bounds="x_bnds")
    ds.y.attrs.update(long_name='Northing', standard_name='projection_y_coordinate', units='m',bounds="y_bnds")
    if include_bounds:
        ds.x_bnds.attrs.update(long_name='Easting boundaries', units='m')
        ds.y_bnds.attrs.update(long_name='Northing boundaries', units='m')
        ds.lon_bnds.attrs.update(long_name='Longitude cell boundaries', units='degrees_east')
        ds.lat_bnds.attrs.update(long_name='Latitude cell boundaries', units='degrees_north')

    ds.attrs.update(
        title = f"prototype EOCIS CHUK grid at {spacing_m}m resolution",
        institution = "EOCIS CHUK",
        date_created = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        version = version,
        Convention = "CF-1.8"
    )

    return ds


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()

    # default extent to somewhat (but not completely) encompass bounding box from https://eocis.org/wp-content/uploads/2023/05/EOCIS-CHUK-Format-AO-For-Issue.pdf
    # Min.Lat 47.5°N Min.Long 12.6°W(-12.6°)
    # Max.Lat 61.3°N Max.Long 3.4°E
    parser.add_argument("output_path", help="path to write the output netcdf4 file to")
    parser.add_argument("--resolution", type=int, help="grid resolution in metres", default=1000)
    parser.add_argument("--min-northing", type=int, help="minimum northing (metres)", default=-267000)
    parser.add_argument("--max-northing", type=int, help="maximum northing (metres)", default=1250000)
    parser.add_argument("--min-easting", type=int, help="minimum easting (metres)", default=-332000)
    parser.add_argument("--max-easting", type=int, help="maximum easting (metres)", default=765000)
    parser.add_argument("--precision", help="set output precision to single or double", default="single")
    parser.add_argument("--version", help="set the version of the grid as an attribute in the output file", default="0.4 (provisional)")
    parser.add_argument("--include-bounds", action="store_true", help="include lat/lon bnds in the generated grid")

    args = parser.parse_args()

    ds = create_grid(min_n=args.min_northing, min_e=args.min_easting,
                     max_n=args.max_northing, max_e=args.max_easting,
                     spacing_m=args.resolution, version=args.version,
                     include_bounds=args.include_bounds)

    precision = "float32" if args.precision == "single" else "float64"

    encoding = {
        "lat": {"dtype": precision, "zlib": True, "complevel": 5},
        "lon": {"dtype": precision, "zlib": True, "complevel": 5}
    }

    if args.include_bounds:
        encoding.update({
            "lat_bnds": {"dtype": precision, "zlib": True, "complevel": 5, "_FillValue": None},
            "lon_bnds": {"dtype": precision, "zlib": True, "complevel": 5, "_FillValue": None},
            "y_bnds":  {"_FillValue": None},
            "x_bnds": {"_FillValue": None}
        })

    ds.to_netcdf(args.output_path, encoding=encoding)
