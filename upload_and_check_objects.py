import csv
import itertools
import warnings
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
# noinspection PyUnresolvedReferences
from astropy.units import deg as u_deg, hourangle, arcsecond
import psycopg2
from typing import Optional
import logging
from os import getenv
from dotenv import load_dotenv

load_dotenv()   # load variables from .env file
logging.basicConfig(filename='objects.log', level=logging.INFO)

Simbad.TIMEOUT = 120
Vizier.TIMEOUT = 120


class TerminateScript(Exception):
    pass


def recover_gaia_name(name):
    return None if name is None else f'Gaia DR3 {name}'


def recover_ucac4_name(name):
    return None if name is None else f'UCAC4 {name}'


def get_vsx_name(names):
    for name in filter(lambda _name: _name is not None, names):
        try:
            table_list = Vizier.query_object(name, catalog='VSX')
        except Exception as e:
            logging.warning(f' !!!!!!!!!!!! Vizier connection error {e}\n Retry...')
            table_list = Vizier.query_object(name, catalog='VSX')
        if len(table_list) < 1:
            continue
        elif len(table_list) > 1:
            raise TerminateScript(f'{name}: vsx has returned more than 1 object')
        return table_list[0]['Name'].value[0]
    return None


def _get_simbad_object(name) -> Optional[tuple[str, SkyCoord]]:
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        try:
            so = Simbad.query_object(name)
        except Exception as e:
            logging.warning(f' !!!!!!!!!!!! Simbad connection error {e}\n Retry...')
            so = Simbad.query_object(name)
        if so is None:
            return None
        return so['MAIN_ID'][0], SkyCoord(ra=so['RA'][0], dec=so['DEC'][0], unit=[hourangle, u_deg])


def get_simbad_object(names: list[str]) -> Optional[tuple[str, SkyCoord]]:
    for name in filter(lambda _name: _name is not None, names):
        simbad_obj = _get_simbad_object(name)
        if simbad_obj is not None:
            return simbad_obj
    return None


def parse_row(row):
    """
    Parse a row from the CSV file and retrieve necessary information.
    Args:
        row (list): A list of values from the CSV row.
    Returns:
        tuple: A tuple containing the processed values ready for database insertion.
    """
    old_id, gaia_name, ucac4_name, apass_name, vsx_name, ra_icrs, de_icrs, object_class, object_description = row

    # Normalize the values, setting any placeholders ('--') to None
    gaia_name, ucac4_name, apass_name, vsx_name, object_class, object_description = \
        [None if x == '--' else x.strip()
         for x in [gaia_name, ucac4_name, apass_name, vsx_name, object_class, object_description]]

    simbad_name = None

    # Retrieve the Simbad object based on the available names
    simbad_obj = get_simbad_object([recover_gaia_name(gaia_name),  recover_ucac4_name(ucac4_name), simbad_name,
                                    apass_name, vsx_name])

    radec = SkyCoord(ra_icrs + ' ' + de_icrs, unit=(u_deg, u_deg), frame='icrs')

    if simbad_obj is not None:
        # Check coords separations (between simbad.radec and radec)
        simbad_name = simbad_obj[0]
        sep_tol = 1 * arcsecond
        radec_simbad = simbad_obj[1]
        sep = radec.separation(radec_simbad)
        logging.info(f'{simbad_name}: separation = {sep.to(arcsecond):.3f}')

        if sep > sep_tol:
            with open('objects.log', 'a') as f:
                f.write(f'Simbad {simbad_name} gaia {gaia_name}: The coordinates discrepancy is considerable. '
                        f'Separation = {sep.to(arcsecond)}\n'
                        f'radec_icrs = {ra_icrs} {de_icrs}\n')

    if vsx_name is None:
        vsx_name = get_vsx_name([simbad_name, recover_gaia_name(gaia_name), recover_ucac4_name(ucac4_name)])
    if vsx_name is not None:
        logging.info(f'----------------------{simbad_name}: vsx_name = {vsx_name}')

    coordequ_str = f'\'({radec.icrs.ra.deg}d, {radec.icrs.dec.deg}d)\''

    simbad_name, gaia_name, apass_name, vsx_name, object_class, object_description = \
        ['NULL' if x is None else f'\'{x}\''
         for x in [simbad_name, gaia_name, apass_name, vsx_name, object_class, object_description]]

    return (old_id, gaia_name, simbad_name, ucac4_name, apass_name, vsx_name, coordequ_str, object_class,
            object_description)


def do_your_stuff(objects_filename_csv: str):
    # Connect to the database using credentials from environment variables
    conn = psycopg2.connect(
        host=getenv("DB_HOST"),
        dbname=getenv("DB_NAME"),
        user=getenv("DB_WRITER_USER"),
        password=getenv("DB_WRITER_PASS")
    )
    cursor = conn.cursor()

    # Open the CSV file containing object data
    with open(objects_filename_csv) as f:
        csv_reader = csv.reader(f)
        # Define the range of rows to process
        n_skip = 10
        n_max = 11
        csv_reader_sliced = itertools.islice(csv_reader, n_skip, n_max)  # iterator of the next n_max-n_skip rows
        for row in csv_reader_sliced:
            old_id, gaia_name, simbad_name, ucac4_name, apass_name, vsx_name, coordequ_str, object_class, \
                object_description = parse_row(row)

            # Check if an object with the same `gaia_name` already exists in the database
            query = 'SELECT * FROM objects WHERE gaia_name = %s'
            cursor.execute(query, (gaia_name,))
            row = cursor.fetchone()
            if len(row) > 0:
                object_id = row[0]
                logging.info(f'object with gaia_name={gaia_name} is already present in the DB with id={object_id}')

                # Update existing object record with new values for `vsx_name` and `simbad_name`
                query = '''
                    UPDATE objects 
                    SET vsx_name = %s, simbad_name = %s 
                    WHERE id = %s
                '''
                cursor.execute(query, (vsx_name, simbad_name, object_id))

            else:
                # Insert a new object record if it does not already exist
                query = '''
                    INSERT INTO objects (gaia_name, simbad_name, ucac4_name, apass_name, vsx_name, 
                                         coordequ, class, description, id_arch) 
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                '''
                cursor.execute(query, (gaia_name, simbad_name, ucac4_name, apass_name,
                                       vsx_name, coordequ_str, object_class, object_description, old_id))
            # Commit the transaction after each insert or update
            # This ensures data is saved incrementally, which is useful if the process is interrupted
            conn.commit()
    cursor.close()
    conn.close()


if __name__ == '__main__':
    do_your_stuff('CSV/Object.csv')
