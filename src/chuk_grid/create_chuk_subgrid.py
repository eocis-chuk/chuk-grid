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

from .latlon_to_eastingnorthing import calculate_eastings_northings
from .create_chuk_grid import create_grid

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="path to an input file in netcdf4 format")
    parser.add_argument("output_path",
                        help="path to output a file with a chuk grid subset that encloses the area of the input file")
    parser.add_argument("--resolution", type=int, help="grid resolution in metres", default=100)
    parser.add_argument("--version", help="set the version of the grid as an attribute in the output file",
                        default="0.4 (provisional)")
    parser.add_argument("--min-northing", type=int, help="minimum northing (metres)", default=-267000)
    parser.add_argument("--max-northing", type=int, help="maximum northing (metres)", default=1250000)
    parser.add_argument("--min-easting", type=int, help="minimum easting (metres)", default=-332000)
    parser.add_argument("--max-easting", type=int, help="maximum easting (metres)", default=765000)

    args = parser.parse_args()

    input_ds = xr.open_dataset(args.input_path)

    lat_min = float(input_ds.lat.min())
    lat_max = float(input_ds.lat.max())
    lon_min = float(input_ds.lon.min())
    lon_max = float(input_ds.lon.max())


    min_e, max_e, min_n, max_n = calculate_eastings_northings(lat_min, lon_min, lat_max, lon_max)

    # does the bounding box overlap with the min/max
    if max_e < args.min_easting or min_e > args.max_easting or max_n < args.min_northing or min_n > args.max_northing:
        print("No overlap with CHUK bounding box, no output file produced")
        return

    output_ds = create_grid(min_n, min_e, max_n, max_e, spacing_m=args.resolution,
                            version=args.version, include_bounds=False)

    encoding = {
        "lat": {"dtype": "float32", "zlib": True, "complevel": 5},
        "lon": {"dtype": "float32", "zlib": True, "complevel": 5}
    }


    output_ds.to_netcdf(args.output_path, encoding=encoding)

if __name__ == '__main__':
    main()