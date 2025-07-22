#!/usr/bin/env python

import sqlite3
import pandas as pd
import sys


def find_field_season_radec(ra, dec, con):

    """
    Query the sky_extents table to find all fields/seasons containing the given
    (ra, dec).

    Parameters:
    -----------
    ra : float
        Right Ascension in degrees (0 ≤ ra < 360)
    dec : float
        Declination in degrees (-90 ≤ dec ≤ +90)
    con : db connection handle

    Returns:
    --------
    matches : list of dict
        List of matching rows where the given position falls within the extent
    """

    if ra < 0:
        exit("ERROR: Please provide RA>0 and RA<360")

    if ra > 180:
        ra180 = 360 - ra
    else:
        ra180 = ra

    QUERY_FIELD_SEASON_RADEC = """
    select NAME, TYPE from SKY_EXTENTS
    WHERE (CROSSRA0='N' AND ({RA} BETWEEN RACMIN and RACMAX) AND ({DEC} BETWEEN DECCMIN and DECCMAX)) OR
          (CROSSRA0='Y' AND ({RA180} BETWEEN RACMIN-360 and RACMAX) AND ({DEC} BETWEEN DECCMIN and DECCMAX))
    """

    query = QUERY_FIELD_SEASON_RADEC.format(RA=ra, DEC=dec, RA180=ra180)

    cursor = con.cursor()
    cursor.execute(query)
    rows = cursor.fetchall()
    if not rows:
        return pd.DataFrame(columns=["name", "type"])
    return pd.DataFrame(rows, columns=["name", "type"])


# Open db connection
db_path = "fields.db"
con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
position_csv = sys.argv[1]
df = pd.read_csv(position_csv)
for index, row in df.iterrows():
    ra = row["RA"]
    dec = row["DEC"]
    print(f"# Row {index}: RA = {ra}, Dec = {dec}")
    r = find_field_season_radec(ra, dec, con)
    print(r)
    print("# --- ")
exit()

#con = sqlite3.connect(db_path)

results = find_field_season_radec(32.69245, -51.017, con)
print(results)
results = find_field_season_radec(15, -25, con)
print(results)
exit()
for r in results:
    print(f"{r['type']} '{r['name']}' contains the point.")
