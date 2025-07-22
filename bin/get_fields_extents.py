#!/usr/bin/env python

from spt3g_cutter.fields import get_all_fields, get_field_extent, get_all_seasons
import sqlite3


conn = sqlite3.connect("fields.db")
cursor = conn.cursor()

# Updated schema with crosses_ra0
cursor.execute("""
    CREATE TABLE IF NOT EXISTS sky_extents (
        name TEXT PRIMARY KEY,
        type TEXT,
        racmin REAL,
        racmax REAL,
        deccmin REAL,
        deccmax REAL,
        crossra0 TEXT
    );
""")

for field in get_all_fields():
    ra, dec = get_field_extent(field)
    racmin, racmax = ra
    deccmin, deccmax = dec
    if racmin > racmax:
        crossra0 = 'Y'
    else:
        crossra0 = 'N'
    cursor.execute("""
        INSERT OR REPLACE INTO sky_extents (
            name, type, racmin, racmax, deccmin, deccmax, crossra0
        ) VALUES (?, ?, ?, ?, ?, ?, ?);
    """, (field, "field", racmin, racmax, deccmin, deccmax, crossra0))

for field in get_all_seasons():
    ra, dec = get_field_extent(field)
    racmin, racmax = ra
    deccmin, deccmax = dec
    if racmin > racmax:
        crossra0 = 'Y'
    else:
        crossra0 = 'N'
    cursor.execute("""
        INSERT OR REPLACE INTO sky_extents (
            name, type, racmin, racmax, deccmin, deccmax, crossra0
        ) VALUES (?, ?, ?, ?, ?, ?, ?);
    """, (field, "season", racmin, racmax, deccmin, deccmax, crossra0))

conn.commit()
conn.close()
