import xarray as xr
import numpy as np
import argparse
import uuid
import pyproj
import numpy as np

# http://cfconventions.org/cf-conventions/

parser = argparse.ArgumentParser()
parser.add_argument("input_path",help="path to input CHUK file")
parser.add_argument("output_path",help="path to output CHUK file")

uuid = str(uuid.uuid4())

args = parser.parse_args()

print(f"Creating metadata example: {args.input_path} => {args.output_path}")

ds = xr.open_dataset(args.input_path)

ds['crsWGS84'] = xr.DataArray(np.int32(0), dims=[], attrs=pyproj.Proj(4326).crs.to_cf())

ds.attrs["summary"] = "A summary of this dataset"
ds.attrs["license"] = "Creative Commons Licence by attribution (https://creativecommons.org/licenses/by/4.0/)"
ds.attrs["uuid"] = uuid
ds.attrs["spatial_resolution"] = "100 m"
ds.attrs["history"] = "describe the history of this product"
ds.attrs["title"] = "title for this dataset"
ds.attrs["comment"] = "a useful comment about this dataset"
ds.attrs["Conventions"] = "CF-1.8"

ds.attrs["creator_url"] = "the creator's web URL"
ds.attrs["creator_name"] = "the creator's name"
ds.attrs["creator_email"] = "general email address for creator NOT a named individual"
ds.attrs["creator_processing_institution"]  = "the creator's institution"

ds.attrs["creation_date"] = "2022-07-16T03:57:11Z"
ds.attrs["publisher_url"] = "https://eocis.org"
ds.attrs["publisher_name"] = "EOCIS"
ds.attrs["publisher_email"] = "EOCIS@reading.ac.uk"
ds.attrs["acknowledgement"] = "Funded by UK EOCIS. Use of these data should acknowledge EOCIS"

# create a dummy continuous fraction variable

fraction_array = np.zeros(ds["lat"].shape,dtype=np.float32)
ds["land_fraction"] = xr.DataArray(fraction_array,dims=("y","x"))
ds["land_fraction"].attrs["standard_name"] = "land_area_fraction"
ds["land_fraction"].attrs["long_name"] = "fraction of cell covered by land"
ds["land_fraction"].attrs["source"] = "more details about the source of this variable"
ds["land_fraction"].attrs["comment"] = "if applicable, helpful comments"
ds["land_fraction"].attrs["units"] = "1"
ds["land_fraction"].attrs["grid_mapping"] = "crsOSGB: x y crsWGS84: lat lon"
ds["land_fraction"].attrs["coordinates"] = "lat lon"

# create a dummy mask type field
mask_array = np.zeros(ds["lat"].shape,dtype=np.int16)
ds["land_class"] = xr.DataArray(mask_array,dims=("y","x"))
ds["land_class"].attrs["standard_name"] = "area_type"
ds["land_class"].attrs["long_name"] = "predominantly land or sea"
ds["land_class"].attrs["flag_masks"] = [0x01, 0x02]
# http://cfconventions.org/Data/area-type-table/current/build/area-type-table.html
ds["land_class"].attrs["flag_meanings"] = "sea land"
ds["land_class"].attrs["source"] = "more details about the source of this variable"
ds["land_class"].attrs["comment"] = "if applicable, helpful comments"
ds["land_class"].attrs["grid_mapping"] = "crsOSGB: x y crsWGS84: lat lon"
ds["land_class"].attrs["coordinates"] = "lat lon"

encoding={
    "land_fraction":{
        "chunksizes":[500,500],
        "dtype":"float32",
        "zlib":True,"complevel":5
    },
    "land_class": {
        "chunksizes":[500,500],
        "dtype":"byte",
        "zlib":True,"complevel":5
    },
    "lat": {
        "chunksizes":[500,500],
        "dtype":"float32",
        "zlib":True,"complevel":5
    },
    "lon": {
        "chunksizes":[500,500],
        "dtype":"float32",
        "zlib":True,"complevel":5
    }
}

ds.to_netcdf(args.output_path,encoding=encoding)