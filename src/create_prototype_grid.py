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


class BNG(object):

    @staticmethod
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
        n_list = []
        n = min_n
        idx = 0
        while n <= max_n:
            n_list.append(n)
            idx += 1
            n = min_n + spacing_m * idx

        e_list = []
        e = min_e
        idx = 0
        while e <= max_e:
            e_list.append(e)
            idx += 1
            e = min_e + spacing_m * idx

        n_val = np.flip(np.array(n_list, dtype=np.float64))  # arrange from north to south
        e_val = np.array(e_list, dtype=np.float64)

        n_da = xr.DataArray(data=n_val, dims=("n",)).expand_dims({"e": len(e_list)}).transpose("n", "e")
        e_da = xr.DataArray(data=e_val, dims=("e",)).expand_dims({"n": len(n_list)})

        northings = n_da.values.flatten()
        eastings = e_da.values.flatten()

        transformer = pyproj.Transformer.from_crs(27700, 4326)
        r = transformer.transform(eastings, northings)
        lats = r[0]
        lons = r[1]

        lats_da = xr.DataArray(data=np.reshape(lats, (len(n_list), len(e_list))), dims=("nj", "ni"))
        lons_da = xr.DataArray(data=np.reshape(lons, (len(n_list), len(e_list))), dims=("nj", "ni"))
        n_da = xr.DataArray(data=n_val, dims=("nj",))
        e_da = xr.DataArray(data=e_val, dims=("ni",))

        ds = xr.Dataset()
        ds["latitude"] = lats_da
        ds["longitude"] = lons_da
        ds["northings"] = n_da
        ds["eastings"] = e_da
        ds.attrs["comment"] = f"prototype EOCIS CHUK grid at {spacing_m}m resolution"
        ds.attrs["institution"] = "EOCIS CHUK"
        ds.attrs["creation_date"] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")  # 2019-09-02T15:19:26
        return ds


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()

    # default extent taken from tas_hadukgrid_uk_1km_ann-30y_198101-201012.nc
    parser.add_argument("output_path", help="path to write the output netcdf4 file to")
    parser.add_argument("--resolution", type=int, help="grid resolution in metres", default=1000)
    parser.add_argument("--min-northing", help="minimum northing (metres)", default=-199500)
    parser.add_argument("--max-northing", help="minimum northing (metres)", default=1249500)
    parser.add_argument("--min-easting", help="minimum easting (metres)", default=-199500)
    parser.add_argument("--max-easting", help="maximum easting (metres)", default=699500)

    args = parser.parse_args()

    ds = BNG.create_grid(min_n=args.min_northing, min_e=args.min_easting,
                         max_n=args.max_northing, max_e=args.max_easting,
                         spacing_m=args.resolution)
    ds.to_netcdf(args.output_path)
