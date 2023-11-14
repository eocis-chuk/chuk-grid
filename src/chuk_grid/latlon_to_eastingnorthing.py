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

import numpy as np
import pyproj

def calculate_eastings_northings(min_lat, min_lon, max_lat, max_lon):
    """
    Get the BNG eastings/northings bounding box that encloses the given lat/lon box

    :param min_lat: minimum lat of area extent
    :param min_lon: minimum lon of area extent
    :param max_lat: maximum lat of area extent
    :param max_lon: maximum lon of area extent
    """

    lat_val = np.linspace(max_lat, min_lat, 100)
    lon_val = np.linspace(min_lon, max_lon, 100)

    shape = len(lat_val), len(lon_val)

    lats = np.broadcast_to(lat_val[None].T, shape)
    lons = np.broadcast_to(lon_val, shape)

    transformer = pyproj.Transformer.from_crs(4326,27700)
    eastings, northings = transformer.transform(lats,lons)

    min_n = 1000*np.floor(np.min(northings)/1000)
    max_n = 1000*np.ceil(np.max(northings)/1000)

    min_e = 1000 * np.floor(np.min(eastings) / 1000)
    max_e = 1000 * np.ceil(np.max(eastings) / 1000)

    return min_e, max_e, min_n, max_n


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument("--min-lat", type=float, help="minimum lat")
    parser.add_argument("--max-lat", type=float, help="maximum lat")
    parser.add_argument("--min-lon", type=float, help="minimum lon")
    parser.add_argument("--max-lon", type=float, help="maximum lon")

    args = parser.parse_args()

    min_e, max_e, min_n, max_n = calculate_eastings_northings(min_lat=args.min_lat, min_lon=args.min_lon,
                     max_lat=args.max_lat, max_lon=args.max_lon)

    print("E/N Bounding Box: N: %.0f to %.0f, E: %.0f to %.0f" % (min_n,max_n,min_e,max_e))

