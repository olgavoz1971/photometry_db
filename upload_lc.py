import json
from time import time
import csv
import itertools
from os.path import basename
import psycopg2
import psycopg2.extensions
import logging
from os import getenv
from dotenv import load_dotenv

load_dotenv()   # load variables from .env file
logging.basicConfig(level=logging.INFO)

filename_csv = '/home/olgavoz/CSV/Obs_UPJS_Photometry.csv'

obs_type_dict = {
    1: 7,
    2: 8,
    3: 9,
    4: 99,
    216: 99,
    276: 2,
    604: 5,
    605: 3,
    1135: 4
}

band_dict = {
    1: 'U',
    2: 'B',
    3: 'V',
    4: 'R',
    5: 'I',
    6: 'u_sdss',
    7: 'g_sdss',
    8: 'r_sdss',
    9: 'i_sdss',
    10: 'z_sdss',
    99: 'NO_FILTER'

}


def timeit(func):
    """
        Decorator for measuring the running time of a function
    """
    def measure_time(*args, **kwargs):
        start_time = time()
        result = func(*args, **kwargs)
        logging.info(f'Processing time of {func.__qualname__} is {(time() - start_time): .2f} seconds')
        return result

    return measure_time


class TerminateScript(Exception):
    pass


def clean_fits_filename(filename):
    return filename.strip().replace('.fz', '')


def upload_lc(cursor: psycopg2.extensions.cursor, row: list):
    (phot_id, jd, tframe, mag, mag_err, detrend_mag, detrend_mag_err, type_detrend, num_calib, type_calib, detrend,
     diff_mag, diff_mag_err, aver_comp, num_comp, differential, comp_stars, object_id_arch, obs_type_id, observer,
     exp_time, binning, temperature, quality, preview_link, image_link, privileges) = row

    logging.info(f'{phot_id=}')

    cursor.execute('SELECT id FROM lightcurves WHERE id_arch = %s', (phot_id,))
    row = cursor.fetchone()

    if row is not None:
        logging.info(f'Row with id_phot_arch = {phot_id} is already present in the DB; skipping it.')
        return

    # ----------------- observation_id
    fits_file_name_short = clean_fits_filename(basename(image_link))

    cursor.execute('SELECT id FROM observations WHERE filename = %s', (fits_file_name_short,))
    row = cursor.fetchone()  # Use fetchone() to retrieve a single row

    observation_id = row[0] if row is not None else 'NULL'

    dateobs_str = f'J{jd}-12'  # !!! Postgres interprets Julian dates in its own way

    # ----------------- data_source
    cursor.execute('SELECT id FROM datasources WHERE progid = %s', (observer,))
    row = cursor.fetchone()
    if row is None:
        raise TerminateScript(f'Lightcurve old id {phot_id}: Unknown data source with ProgID: {observer}')

    datasource_id = row[0]  # Extract the datasource_id from the result

    # ----------------- Object
    cursor.execute('SELECT id FROM objects WHERE id_arch = %s', (object_id_arch,))
    row = cursor.fetchone()  # Use fetchone() to retrieve a single row
    if row is None:
        raise TerminateScript(f'Lightcurve old id {phot_id}: Incorrect id_object_arch {object_id_arch}')
    object_id = row[0]

    # ----------------- Photometric band
    # Retrieve the observation type ID from the dictionary
    obs_type_id = obs_type_dict.get(int(obs_type_id), None)
    band = band_dict.get(obs_type_id, None)
    if band is None:
        raise TerminateScript(f'Unknown band.id "observation" (old) = {obs_type_id}')

    # TODO: Temporary solution until photosys concept is developed
    cursor.execute('SELECT id FROM photosys WHERE band = %s', (band,))
    row = cursor.fetchone()
    if row is None:
        raise TerminateScript(f'Lightcurve old id {phot_id}: Band {band} not found in photosys table')

    photosys_id = row[0]  # Extract the photosys_id from the result

    # Insert a new lightcurve entry into the DB
    cursor.execute(
        '''
        INSERT INTO lightcurves (
            object_id, dateobs, magnitude, mag_err, mag_diff, mag_diff_err,
            datasource_id, photosys_id, observation_id, id_arch, quality
        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ''',
        (object_id, dateobs_str, mag, mag_err, diff_mag, diff_mag_err,
         datasource_id, photosys_id, observation_id, phot_id, quality)
    )

    # Retrieve the last inserted lightcurve ID
    cursor.execute('SELECT id FROM lightcurves ORDER BY id DESC LIMIT 1')
    lightcurve_id = cursor.fetchone()[0]

    # ---------------- Upload related Metadata --------------------------

    comp_stars_list = eval(comp_stars)  # Convert the string representation of comparison stars into a dictionary
    params_dict = {
        'detrend_mag': float(detrend_mag),
        'detrend_mag_err': float(detrend_mag_err),
        'type_detrend': type_detrend,
        'num_calib': int(num_calib),
        'type_calib': type_calib,
        'detrend': detrend,
        'aver_comp': float(aver_comp),
        'num_comp': int(num_comp),
        'differential': int(differential),
        'exp_time': float(exp_time),
        'binning': binning,
        'temperature': float(temperature),
        'comp_stars': comp_stars_list
    }
    js = json.dumps(params_dict)

    logging.info(f'{lightcurve_id}')

    cursor.execute(
        '''
        INSERT INTO metadata_lc_json (lightcurve_id, data)
        VALUES (%s, %s)
        ''',
        (lightcurve_id, js)
    )


@timeit
def do_your_stuff():
    # Connect to the database using credentials from environment variables
    conn = psycopg2.connect(
        host=getenv("DB_HOST"),
        dbname=getenv("DB_NAME"),
        user=getenv("DB_WRITER_USER"),
        password=getenv("DB_WRITER_PASS")
    )
    cursor = conn.cursor()

    with open(filename_csv) as f:
        csv_reader = csv.reader(f)
        # row_header = next(csv_reader)  # header
        # print(row_header)
        # We have 1,466,178 rows in Obs_UPJS_Photometry.csv. The first row is the header.
        # In this example, we attempt to upload only one row, the last one:
        n_skip = 1466177  # skip first n_skip rows
        n_max = n_skip + 1  # the last number, if None, then iteration continues until the iterator is exhausted
        csv_reader_sliced = itertools.islice(csv_reader, n_skip, n_max)  # slice to iterate over the next n_max-n rows
        for row in csv_reader_sliced:
            upload_lc(cursor, row)
            # Commit the transaction after each insert or update
            # This ensures data is saved incrementally, which is useful if the process is interrupted
            conn.commit()
    cursor.close()
    conn.close()


if __name__ == '__main__':
    do_your_stuff()
