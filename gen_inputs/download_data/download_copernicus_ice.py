from copernicusmarine import subset
from datetime import datetime, timedelta
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# -------- CONFIG --------
# Model Output
dataset_id = "cmems_mod_arc_phy_anfc_nextsim_hm"
variables = [
    "si_ridge_ratio", "siage", "sialb", "siconc", "siconc_my",
    "siconc_young", "sisnthick", "sithick", "vxsi", "vysi"
]
lat_min, lat_max = 68, 73
lon_min, lon_max = -156, -137
output_directory = "/pscratch/sd/b/bundzis/External_data/Copernicus/Sea_ice"
#os.makedirs(output_directory, exist_ok=True)

# -------- TIME CHUNKS --------
def generate_month_ranges(start, end):
    current = start
    while current <= end:
        next_month = (current.replace(day=28) + timedelta(days=4)).replace(day=1)
        yield current, min(next_month, end)
        current = next_month

# -------- DOWNLOAD FUNCTION --------
def download_month(start_date, end_date):
    out_file = f"CMEMS_{start_date.strftime('%Y_%m')}.nc"
    out_path = os.path.join(output_directory, out_file)

    if os.path.exists(out_path):
        return f"✅ Already exists: {out_file}"

    attempt = 0
    while attempt < 3:
        try:
            subset(
                dataset_id=dataset_id,
                variables=variables,
                minimum_longitude=lon_min,
                maximum_longitude=lon_max,
                minimum_latitude=lat_min,
                maximum_latitude=lat_max,
                start_datetime=start_date.strftime('%Y-%m-%dT00:00:00'),
                end_datetime=end_date.strftime('%Y-%m-%dT00:00:00'), # see if this prevents duplicate times if set to 23 hours, otherwise set back to 0
                output_filename=out_file,
                output_directory=output_directory
            )
            return f"✅ Downloaded: {out_file}"
        except Exception as e:
            attempt += 1
            wait_time = 60 * attempt
            print(f"⚠️ Attempt {attempt} failed for {out_file}: {e}")
            print(f"⏳ Retrying in {wait_time} seconds...")
            time.sleep(wait_time)
    return f"❌ Failed after 3 attempts: {out_file}"

# -------- MAIN EXECUTION --------
if __name__ == "__main__":
    start = datetime(2019, 8, 1)
    end   = datetime(2024, 12, 31)

    month_ranges = list(generate_month_ranges(start, end))

    max_threads = 2  # Be respectful to the CMEMS servers
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = [executor.submit(download_month, start, end) for start, end in month_ranges]

        for future in as_completed(futures):
            print(future.result())

