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


def create_grid(min_n, min_e, max_n, max_e, spacing_m):
    """
    Create an xarray Dataset describing a grid on EPSG:27700

    :param min_n: minimum BNG northing value (m) of area extent
    :param min_e: minimum BNG easting value (m) of area extent
    :param max_n: maximum BNG northing value (m) of area extent
    :param max_e: maximum BNG easting value (m) of area extent
    :param spacing_m: grid spacing (m)

    :return: xarray dataset consisting of the following arrays

    latitude    (nj,ni) array of WGS84 latitude values calculated via pyproj
    longitude   (nj,ni) array of WGS84 longitude values calculated via pyproj
    northings   (nj) array contiaining BNG northings
    easting     (ni) array containing BNG eastings
    """

    # Northings from north to south
    n_val = np.linspace(max_n, min_n, int((max_n-min_n)/spacing_m)+1)
    e_val = np.linspace(min_e, max_e, int((max_e-min_e)/spacing_m)+1)

    shape = len(n_val), len(e_val)
    northings = np.broadcast_to(n_val[None].T, shape)
    eastings = np.broadcast_to(e_val, shape)

    transformer = pyproj.Transformer.from_crs(27700, 4326)
    lats, lons = transformer.transform(eastings, northings)

    print("Lat/Lon Bounding Box: Lat: %.4f to %.4f, Lon: %.4f to %.4f" % (np.min(lats),np.max(lats),np.min(lons),np.max(lons)))

    ds = xr.Dataset({
            'lat': (['y', 'x'], lats),
            'lon': (['y', 'x'], lons),
            'crsOSGB': ([], np.int32(0), pyproj.Proj(27700).crs.to_cf())
        },
        coords={
            'x': e_val,
            'y': n_val,
        })

    ds.lat.attrs.update(standard_name='latitude', units='degrees_north')
    ds.lon.attrs.update(standard_name='longitude', units='degrees_east')
    ds.x.attrs.update(long_name='Easting', standard_name='projection_x_coordinate', units='m')
    ds.y.attrs.update(long_name='Northing', standard_name='projection_y_coordinate', units='m')

    ds.attrs.update(
        title = f"prototype EOCIS CHUK grid at {spacing_m}m resolution",
        institution = "EOCIS CHUK",
        date_created = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))

    return ds


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()

    # default extent taken from tas_hadukgrid_uk_1km_ann-30y_198101-201012.nc
    parser.add_argument("output_path", help="path to write the output netcdf4 file to")
    parser.add_argument("--resolution", type=int, help="grid resolution in metres", default=1000)
    parser.add_argument("--min-northing", type=int, help="minimum northing (metres)", default=-199500)
    parser.add_argument("--max-northing", type=int, help="maximum northing (metres)", default=1249500)
    parser.add_argument("--min-easting", type=int, help="minimum easting (metres)", default=-199500)
    parser.add_argument("--max-easting", type=int, help="maximum easting (metres)", default=699500)
    parser.add_argument("--precision", help="set output precision to single or double", default="double")

    args = parser.parse_args()

    ds = create_grid(min_n=args.min_northing, min_e=args.min_easting,
                     max_n=args.max_northing, max_e=args.max_easting,
                     spacing_m=args.resolution)

    precision = "float32" if args.precision == "single" else "float64"

    ds.to_netcdf(args.output_path, encoding={
        "lat": {"dtype": precision, "zlib": True, "complevel": 5},
        "lon": {"dtype": precision, "zlib": True, "complevel": 5}
    })
