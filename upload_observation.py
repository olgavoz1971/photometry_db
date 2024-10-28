import glob
import os
import sys
import warnings

import astropy
# noinspection PyUnresolvedReferences
from astropy.units import deg as u_deg, hourangle, arcsecond

import psycopg2
import psycopg2.extensions
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from numpy import radians
import json
import logging
from os import getenv
from dotenv import load_dotenv

load_dotenv()   # load variables from .env file
logging.basicConfig(level=logging.DEBUG)


imagetype_dict = {'Light Frame': 'light', 'Dark Frame': 'dark'}
frame_dict = {'FK5': 'fk5', 'ICRS': 'icrs'}
band_dict = {'U': 'u', 'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I', 'u': 'u_sdss', 'g': 'g_sdss', 'r': 'r_sdss',
             'i': 'i_sdss', 'z': 'z_sdss'}
stop_list = ['COMMENT', 'HISTORY', 'SIMPLE', 'DATASUM', 'CHECKSUM', 'XTENSION', 'PCOUNT', 'GCOUNT']


class TerminateScript(Exception):
    pass


def clean_fits_filename(filename):
    return filename.strip().replace('.fz', '')


def upload_metadata_json(header: astropy.io.fits.hdu, observation_id: int, cursor):
    params_di = {}
    for key in list(header.keys()):
        if key.upper() in stop_list:
            continue
        if key.startswith('TR') and (key.count('_') == 1):
            continue
        if key.startswith('A_') or key.startswith('B_') or key.startswith('AP_') or key.startswith('BP_'):
            continue
        if key.lstrip() == '':
            continue
        params_di[key] = header[key]
    js = json.dumps(params_di)
    cursor.execute('UPDATE metadata_obs_json SET data = %s WHERE observation_id = %s', (js, observation_id))
    return


def upload_fits(cursor: psycopg2.extensions.cursor, filename_full):
    with fits.open(filename_full) as hdul:
        header = hdul[-1].header

    dateobs_str = header['DATE-OBS'] + 'UTC'
    telescope = header['TELESCOP']
    instrument = header['INSTRUME']
    nx = header['NAXIS1']
    ny = header['NAXIS2']
    logging.info(f'{dateobs_str=}')

    # Calculate pixel scales
    try:
        scale_x = 206265 / header['FOCALLEN'] * header['XPIXSZ'] / 1000 / 3600  # deg/px
        scale_y = 206265 / header['FOCALLEN'] * header['YPIXSZ'] / 1000 / 3600
    except Exception as e:
        logging.error(f'Missing header key?: {str(e)}')
        raise TerminateScript('scale_x, scale_y problem')
    deg_x = nx * scale_x
    deg_y = ny * scale_y
    mode = imagetype_dict[header['IMAGETYP']]
    if mode != 'dark':
        band_orig = header['FILTER']
        # band = f'\'{band_dict.get(band_orig, None)}\''
        band = band_dict.get(band_orig, None)
        if band is None:
            raise TerminateScript(f'Unknown photometric band {band_orig} in {filename_full}')

        wcs = None
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=FITSFixedWarning)
                # Check if there are WCS-related keywords in the header
                if any(key in header for key in ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2']):
                    wcs = WCS(header, relax=True)
                    radec = SkyCoord.from_pixel(nx / 2 - 0.5, ny / 2 - 0.5, wcs=wcs, origin=0)
                    coordequsrc_str = 'wcs'
                else:
                    raise ValueError("WCS information is missing from the FITS header.")

        except Exception as e:
            logging.warning(f'wcs failure: {str(e)}')
            try:
                logging.info('Try OBJCTRA/OBJCTDEC...')
                ra_str = header['OBJCTRA']
                dec_str = header['OBJCTDEC']
                frame = header.get('RADECSYS', 'ICRS')
                radec = SkyCoord(ra_str + ' ' + dec_str, unit=(hourangle, u_deg), frame=frame_dict.get(frame))
                coordequsrc_str = 'objectequ'
            except Exception as e1:
                logging.warning(f'RA DEC failure:{str(e1)}')
                raise TerminateScript('RA DEC failure')

        # Check scale consistency
        if coordequsrc_str == 'wcs':
            scale_tol = 0.05
            scale_x_wcs, scale_y_wcs = proj_plane_pixel_scales(wcs)
            if (2 * abs(scale_x_wcs - scale_x) / (scale_x_wcs + scale_x) > scale_tol) or \
                    (2 * abs(scale_y_wcs - scale_y) / (scale_y_wcs + scale_y) > scale_tol):
                logging.info(
                    f'scale_x_wcs = {scale_x_wcs} scale_x = {scale_x} scale_y_wcs = {scale_y_wcs} scale_y = {scale_y}')
                with open('observations.log', 'a') as f:
                    f.write(f'{filename_full}: The scales differ too much\n'
                            f'scale_x_wcs = {scale_x_wcs} scale_x = {scale_x} '
                            f'scale_y_wcs = {scale_y_wcs} scale_y = {scale_y}\n')
                # raise TerminateScript('The scales differ too much')

        # coordequ_str = f'\'({radec.icrs.ra.deg}d, {radec.icrs.dec.deg}d)\''
        coordequ_str = f'({radec.icrs.ra.deg}d, {radec.icrs.dec.deg}d)'
        # coordequsrc_str = f'\'{coordequsrc_str}\''
        coordequsrc_str = f'{coordequsrc_str}'
    else:
        band = 'NULL'
        coordequ_str = 'NULL'
        coordequsrc_str = 'NULL'

    img_type = 'image'
    accumtime = header['EXPTIME']

    # --------------
    # Fetch telescope ID
    cursor.execute('SELECT id FROM telescopes WHERE name = %s', (telescope,))
    row = cursor.fetchone()
    if row is None:
        raise TerminateScript(f'Unknown telescope {telescope}')
    telescope_id = row[0]

    # Fetch instrument ID
    cursor.execute('SELECT id FROM instruments WHERE name = %s', (instrument,))
    row = cursor.fetchone()
    if row is None:
        raise TerminateScript(f'Unknown instrument {instrument}')
    instrument_id = row[0]

    observation_id = None
    path_to_fits, fits_name = os.path.split(filename_full)
    filename_short = clean_fits_filename(fits_name)

    # Check if the observation already exists
    cursor.execute('SELECT id FROM observations WHERE filename = %s', (filename_short,))
    row = cursor.fetchone()
    if row:
        observation_id = row[0]
        logging.warning(f'DB already contains file {filename_short}. Updating...\n')

    # Insert new observation if it doesn't exist
    if observation_id is None:  # new observation
        cursor.execute(
            '''
            INSERT INTO observations (dateobs, accumtime, telescope_id, instrument_id, fov, mode, band, coordequ,
                                      coordequsrc, type, filename, path_to_fits) 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            ''',
            (dateobs_str, accumtime, telescope_id, instrument_id,
             f'({radians(deg_x)},{radians(deg_y)})', mode, band, coordequ_str,
             coordequsrc_str, img_type, filename_short, path_to_fits)
        )

        # Retrieve the ID of the newly inserted observation
        cursor.execute('SELECT id FROM observations ORDER BY id DESC LIMIT 1')
        row = cursor.fetchone()
        observation_id = row[0]

    else:
        query = '''
            UPDATE observations 
            SET dateobs = %s, accumtime = %s, telescope_id = %s, 
                instrument_id = %s, fov = %s, mode = %s, 
                band = %s, coordequ = %s, coordequsrc = %s, 
                type = %s, filename = %s, path_to_fits = %s 
            WHERE id = %s
            '''
        cursor.execute(query, (
            dateobs_str, accumtime, telescope_id, instrument_id,
            f'({radians(deg_x)},{radians(deg_y)})', mode,
            band, coordequ_str, coordequsrc_str, img_type,
            filename_short, path_to_fits, observation_id
        ))

    upload_metadata_json(header, observation_id, cursor)


def do_your_stuff(path_to_dir):
    # Connect to the database using credentials from environment variables
    conn = psycopg2.connect(
        host=getenv("DB_HOST"),
        dbname=getenv("DB_NAME"),
        user=getenv("DB_WRITER_USER"),
        password=getenv("DB_WRITER_PASS")
    )
    cursor = conn.cursor()
    # Create a sorted list of all files found in the specified directory and its subdirectories
    # by walking through the directory structure
    list_of_files = sorted(filter(os.path.isfile, glob.glob(path_to_dir + '/**/*', recursive=True)))

    for filename_full in list_of_files:
        ext = os.path.splitext(filename_full)[-1]
        if (ext != '.fz') and (ext != '.fit') and (ext != '.fits'):
            continue
        logging.info(filename_full)
        upload_fits(cursor, filename_full)
        # Commit the transaction after each insert or update
        # This ensures data is saved incrementally, which is useful if the process is interrupted
        conn.commit()
    cursor.close()
    conn.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('Usage: ./upload_observation.py /home/skvo/data/upjs/Alica')
    mypath_ = sys.argv[1]
    do_your_stuff(mypath_)
