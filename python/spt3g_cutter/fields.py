# Function lifted from:
# https://github.com/SouthPoleTelescope/spt3g_software/blob/8a24bc5c73ff3654bfa70e460ff25f0911c2ee77/sources/python/fields.py

"""
fields.py

Functions for handling observing field definitions and corresponding source
lists.
"""
import numpy as np

core_G3Units_deg = 0.017453292519943295
core_G3Units_rad = 1


# from spt3g import core
__all__ = [
    "get_field_name",
    "get_field_season",
    "get_all_fields",
    "get_all_seasons",
    "get_fields_in_season",
    "is_season",
    "is_subfield",
    "get_field_extent",
    "get_boresight_field_extent",
    "get_field_center",
    "sample_field",
]


def get_field_name(field):
    """
    Return the standard name for the given observing field.

    Arguments
    ---------
    field : str
        Field name, may be an alias

    Returns
    -------
    field : str
        Standard field name
    """
    field = field.lower()
    if field in ["spt3g", "winter"]:
        field = "spt3g-winter"
    elif field in ["summera", "summer-a", "summer", "spt3g-summer"]:
        field = "spt3g-summera"
    elif field in ["summerb", "summer-b"]:
        field = "spt3g-summerb"
    elif field in ["summerc", "summer-c"]:
        field = "spt3g-summerc"
    elif field in ["sgrastar"]:
        field = "spt3g-sgrastar"
    elif field in ["ra17h46dec-29"]:
        field = "j1745-290"
    elif field in ["galaxy"]:
        field = "spt3g-galaxy"
    elif field in ["100d", "ra23h30dec-55"]:
        field = "sptpol-100d"
    elif field in ["500d", "ra0hdec-57.5"]:
        field = "sptpol-500d"
    elif field in ["kids", "ra0p75hdec-31"]:
        field = "sptpol-kids"
    elif field in ["ecs", "sptpol-summer"]:
        field = "sptpol-ecs"
    elif field in ["ecs-back", "sptpol-summer-back"]:
        field = "sptpol-ecs-back"
    elif field in ["widea", "wide-a"]:
        field = "spt3g-widea"
    elif field in ["wideb", "wide-b"]:
        field = "spt3g-wideb"
    elif field in ["widec", "wide-c"]:
        field = "spt3g-widec"
    elif field in ["wided", "wide-d"]:
        field = "spt3g-wided"
    elif field in ["widee", "wide-e"]:
        field = "spt3g-widee"
    elif field in ["widef", "wide-f"]:
        field = "spt3g-widef"
    elif field in ["wideg", "wide-g"]:
        field = "spt3g-wideg"
    elif field in ["wideh", "wide-h"]:
        field = "spt3g-wideh"
    elif field in ["widei", "wide-i"]:
        field = "spt3g-widei"
    elif field in ["a1"]:
        field = "ra4h03dec-21"
    elif field in ["a2"]:
        field = "ra4h03dec-23"
    elif field in ["a3"]:
        field = "ra4h03dec-25"
    elif field in ["a4"]:
        field = "ra4h03dec-27"
    elif field in ["b1"]:
        field = "ra21h53dec-21"
    elif field in ["b2"]:
        field = "ra21h53dec-23"
    elif field in ["b3"]:
        field = "ra21h53dec-25"
    elif field in ["b4"]:
        field = "ra21h53dec-27"
    elif field in ["c1"]:
        field = "ra12hdec-21"
    elif field in ["c2"]:
        field = "ra12hdec-23"
    elif field in ["c3"]:
        field = "ra12hdec-25"
    elif field in ["c4"]:
        field = "ra12hdec-27"
    elif field in ["d1"]:
        field = "ra21h10dec-29.75"
    elif field in ["d2"]:
        field = "ra21h10dec-33.25"
    elif field in ["d3"]:
        field = "ra21h10dec-36.75"
    elif field in ["d4"]:
        field = "ra21h10dec-40.25"
    elif field in ["e1"]:
        field = "ra18h40dec-45.5"
    elif field in ["e2"]:
        field = "ra18h40dec-52.5"
    elif field in ["e3"]:
        field = "ra18h40dec-59.5"
    elif field in ["e4"]:
        field = "ra18h40dec-66.5"
    elif field in ["f1"]:
        field = "ra12h30dec-45"
    elif field in ["f2"]:
        field = "ra12h30dec-51"
    elif field in ["g1"]:
        field = "ra6h40dec-66.5"
    elif field in ["g2"]:
        field = "ra6h40dec-72.5"
    elif field in ["g3"]:
        field = "ra6h40dec-77.5"
    elif field in ["h1"]:
        field = "ra23hdec-72.5"
    elif field in ["h2"]:
        field = "ra23hdec-77.5"
    elif field in ["i1"]:
        field = "ra14h20dec-72.5"
    elif field in ["i2"]:
        field = "ra14h20dec-77.5"
    elif field in ["edfs"]:
        field = "ra4h04dec-48"

    return field


# dictionary of seasons associated with each field
field_seasons = {
    # # sptsz fields
    # "ra21hdec-42.5": "sptsz",
    # "ra21hdec-50": "sptsz",
    # "ra21hdec-60": "sptsz",
    # "ra22h30dec-55": "sptsz",
    # "ra23hdec-45": "sptsz",
    # "ra23hdec-62.5": "sptsz",
    # "ra23h30dec-55": "sptsz",
    # "ra0h50dec-50": "sptsz",
    # "ra1hdec-42.5": "sptsz",
    # "ra1hdec-60": "sptsz",
    # "ra2h30dec-50": "sptsz",
    # "ra3h30dec-42.5": "sptsz",
    # "ra3h30dec-60": "sptsz",
    # "ra4h10dec-50": "sptsz",
    # "ra5h30dec-45": "sptsz",
    # "ra5h30dec-55": "sptsz",
    # "ra6hdec-62.5": "sptsz",
    # "ra6h30dec-45": "sptsz",
    # "ra6h30dec-55": "sptsz",
    # # sptpol fields
    # "sptpol-100d": "sptpol-100d",
    # "sptpol-500d": "sptpol-500d",
    # "sptpol-kids": "sptpol-kids",
    # "ra23hdec-25": "sptpol-ecs",
    # "ra23hdec-35": "sptpol-ecs",
    # "ra1hdec-25": "sptpol-ecs",
    # "ra1hdec-35": "sptpol-ecs",
    # "ra3hdec-25": "sptpol-ecs",
    # "ra3hdec-35": "sptpol-ecs",
    # "ra5hdec-25": "sptpol-ecs",
    # "ra5hdec-35": "sptpol-ecs",
    # "ra11hdec-25": "sptpol-ecs-back",
    # "ra13hdec-25": "sptpol-ecs-back",
    # spt3g fields
    "ra0hdec-44.75": "spt3g-winter",
    "ra0hdec-52.25": "spt3g-winter",
    "ra0hdec-59.75": "spt3g-winter",
    "ra0hdec-67.25": "spt3g-winter",
    "ra5hdec-24.5": "spt3g-summera",
    "ra5hdec-31.5": "spt3g-summera",
    "ra5hdec-38.5": "spt3g-summera",
    "ra5hdec-45.5": "spt3g-summera",
    "ra5hdec-52.5": "spt3g-summera",
    "ra5hdec-59.5": "spt3g-summera",
    "ra5hdec-29.75": "spt3g-summera",
    "ra5hdec-33.25": "spt3g-summera",
    "ra5hdec-36.75": "spt3g-summera",
    "ra5hdec-40.25": "spt3g-summera",
    "ra1h40dec-29.75": "spt3g-summerb",
    "ra1h40dec-33.25": "spt3g-summerb",
    "ra1h40dec-36.75": "spt3g-summerb",
    "ra1h40dec-40.25": "spt3g-summerb",
    "ra12h30dec-29.75": "spt3g-summerc",
    "ra12h30dec-33.25": "spt3g-summerc",
    "ra12h30dec-36.75": "spt3g-summerc",
    "ra12h30dec-40.25": "spt3g-summerc",
    "j1745-290": "spt3g-sgrastar",
    "ra17h45dec-29": "spt3g-galaxy",
    "ra17h37dec-32": "spt3g-galaxy",
    "ra17h27dec-35": "spt3g-galaxy",
    "ra4h03dec-21": "spt3g-widea",
    "ra4h03dec-23": "spt3g-widea",
    "ra4h03dec-25": "spt3g-widea",
    "ra4h03dec-27": "spt3g-widea",
    "ra21h53dec-21": "spt3g-wideb",
    "ra21h53dec-23": "spt3g-wideb",
    "ra21h53dec-25": "spt3g-wideb",
    "ra21h53dec-27": "spt3g-wideb",
    "ra12hdec-21": "spt3g-widec",
    "ra12hdec-23": "spt3g-widec",
    "ra12hdec-25": "spt3g-widec",
    "ra12hdec-27": "spt3g-widec",
    "ra21h10dec-29.75": "spt3g-wided",
    "ra21h10dec-33.25": "spt3g-wided",
    "ra21h10dec-36.75": "spt3g-wided",
    "ra21h10dec-40.25": "spt3g-wided",
    "ra18h40dec-45.5": "spt3g-widee",
    "ra18h40dec-52.5": "spt3g-widee",
    "ra18h40dec-59.5": "spt3g-widee",
    "ra18h40dec-66.5": "spt3g-widee",
    "ra12h30dec-45": "spt3g-widef",
    "ra12h30dec-51": "spt3g-widef",
    "ra6h40dec-66.5": "spt3g-wideg",
    "ra6h40dec-72.5": "spt3g-wideg",
    "ra6h40dec-77.5": "spt3g-wideg",
    "ra23hdec-72.5": "spt3g-wideh",
    "ra23hdec-77.5": "spt3g-wideh",
    "ra14h20dec-72.5": "spt3g-widei",
    "ra14h20dec-77.5": "spt3g-widei",
    "ra4h04dec-48": "spt3g-edfs",
}

# dictionary of fields associated with each season
season_fields = {}
for k, v in field_seasons.items():
    d = season_fields.setdefault(v, [])
    d += [k]


def get_all_fields(season_prefix=None):
    """
    Return a list of all available observing fields, optionally applying a
    season filter.

    Arguments
    ---------
    season_prefix : str
        Filter field list to include only those whose corresponding season
        starts with this string.

    Returns
    -------
    fields : list of str
        List of matching fields.
    """
    fields = list(field_seasons)
    if season_prefix is not None:
        fields = [
            f for f in fields if get_field_season(f).startswith(season_prefix)
        ]
    return fields


def get_all_seasons(season_prefix=None, remove_prefix=False):
    """
    Return a list of all available observing seasons, optionally applying a
    season filter.

    Arguments
    ---------
    season_prefix : str
        Filter season list to include only those that start with this string.
    remove_prefix : bool
        If True, remove the prefix from season names before returning.

    Returns
    -------
    seasons : list of str
        List of matching seasons.
    """
    seasons = list(season_fields)
    if season_prefix is not None:
        seasons = [s for s in seasons if s.startswith(season_prefix)]
    if remove_prefix:
        seasons = [s.replace(season_prefix, "") for s in seasons]
    return seasons


def get_field_season(field):
    """
    Return the season in which a subfield is observed.

    Arguments
    ---------
    field : str
        Field name, may be an alias or an observing season itself.

    Returns
    -------
    season : str
        Observing season
    """
    field = get_field_name(field)

    if field in field_seasons.values():
        return field

    return field_seasons[field]


def get_fields_in_season(season):
    """
    Return a list of all fields observed as part of the given season.

    Arguments
    ---------
    season : str
        Season name

    Returns
    -------
    fields : list of str
        Fields in the season
    """
    season = get_field_name(season)
    return season_fields[season]


def is_season(field):
    """
    Return True if the field name corresponds to a season, otherwise False.
    """
    return get_field_name(field) in season_fields


def is_subfield(field):
    """
    Return True if the field name corresponds to a subfield in a season,
    otherwise False.
    """
    return get_field_name(field) in field_seasons


def get_field_extent(
    field,
    ra_pad=3 * core_G3Units_deg,
    dec_pad=2 * core_G3Units_deg,
    sky_pad=True,
):
    """
    Get the extent of the given field.

    RA angles are always given between 0 and 360 degrees, with the left edge of
    the field given first. Dec angles are given between -90 and 90 degrees,
    with the lower edge of the field given first.

    For example, the SPT-3G winter field "ra0hdec-44.75" has ra range (310, 50)
    degrees, and dec range (-47.5, -42) degrees.

    The default padding parameters ensure that the extent of the SPT-3G fields
    covers the area that any bolometer on the focal plane may touch.  Setting
    the padding parameters to 0 instead returns the part of the field that has
    uniformly weighted coverage.  Padding in ra is defined in terms of degrees
    on the sky if ``sky_pad`` is True.

    The same padding considerations are not taken into account for the SPT-SZ
    and SPTpol fields.

    Arguments
    ---------
    field : str
        Name of the field
    ra_pad : float
        Padding to apply to extent in ra, in G3Units.  Values > 0 extend the
        field outward, and values < 0 shrink the field inward.  Set to 0 to
        select the uniformly weighted coverage region of the field.  If
        ``sky_pad`` is True, the padding is determined in terms of true degrees
        on the sky.
    dec_pad : float
        Padding to apply to extent in dec, in G3Units. Values > 0 extend the
        field outward, and values < 0 shrink the field inward.  Set to 0 to
        select the uniformly weighted coverage region of the field.
    sky_pad : bool
        If True, ``ra_pad`` is in terms of true degrees on the sky.  Otherwise,
        it is applied directly to the ra extent without correction.

    Returns
    -------
    ra_range : 2-tuple of floats
        Tuple of (min, max), the field extent in ra, in G3Units.
    dec_range : 2-tuple of floats, in G3Units
        Tuple of (min, max), the field extent in dec, in G3Units.
    """
    field = get_field_name(field)

    extents = {
        # sptsz fields
        "ra21hdec-42.5": ((300, 330), (-45, -40)),
        "ra21hdec-50": ((300, 330), (-55, -45)),
        "ra21hdec-60": ((300, 330), (-65, -55)),
        "ra22h30dec-55": ((330, 345), (-60, -50)),
        "ra23hdec-45": ((330, 360), (-50, -40)),
        "ra23hdec-62.5": ((330, 360), (-65, -60)),
        "ra23h30dec-55": ((345, 360), (-60, -50)),
        "ra0h50dec-50": ((0, 25), (-55, -45)),
        "ra1hdec-42.5": ((0, 30), (-45, -40)),
        "ra1hdec-60": ((0, 30), (-65, -55)),
        "ra2h30dec-50": ((25, 50), (-55, -45)),
        "ra3h30dec-42.5": ((30, 75), (-45, -40)),
        "ra3h30dec-60": ((30, 75), (-65, -55)),
        "ra4h10dec-50": ((50, 75), (-55, -45)),
        "ra5h30dec-45": ((75, 90), (-50, -40)),
        "ra5h30dec-55": ((75, 90), (-60, -50)),
        "ra6hdec-62.5": ((75, 105), (-65, -60)),
        "ra6h30dec-45": ((90, 105), (-50, -40)),
        "ra6h30dec-55": ((90, 105), (-60, -50)),
        "sptsz": ((300, 105), (-65, -40)),
        # sptpol fields
        "sptpol-100d": ((345, 360), (-60, -50)),
        "sptpol-500d": ((330, 30), (-65, -50)),
        "sptpol-kids": ((330, 52.5), (-36, -26)),
        "ra23hdec-25": ((330, 360), (-30, -20)),
        "ra23hdec-35": ((330, 360), (-40, -30)),
        "ra1hdec-25": ((0, 30), (-30, -20)),
        "ra1hdec-35": ((0, 30), (-40, -30)),
        "ra3hdec-25": ((30, 60), (-30, -20)),
        "ra3hdec-35": ((30, 60), (-40, -30)),
        "ra5hdec-25": ((60, 90), (-30, -20)),
        "ra5hdec-35": ((60, 90), (-40, -30)),
        "sptpol-ecs": ((330, 90), (-40, -20)),
        "ra11hdec-25": ((150, 180), (-30, -20)),
        "ra13hdec-25": ((180, 210), (-30, -20)),
        "sptpol-ecs-back": ((150, 210), (-30, -20)),
        # spt3g fields
        "ra0hdec-44.75": ((310, 50), (-47.5, -42)),
        "ra0hdec-52.25": ((310, 50), (-55, -49.5)),
        "ra0hdec-59.75": ((310, 50), (-62.5, -57)),
        "ra0hdec-67.25": ((310, 50), (-70, -64.5)),
        "spt3g-winter": ((310, 50), (-70, -42)),
        "ra5hdec-24.5": ((50, 100), (-27, -22)),
        "ra5hdec-31.5": ((50, 100), (-34, -29)),
        "ra5hdec-38.5": ((50, 100), (-41, -36)),
        "ra5hdec-45.5": ((50, 100), (-48, -43)),
        "ra5hdec-52.5": ((50, 100), (-55, -50)),
        "ra5hdec-59.5": ((50, 100), (-62, -57)),
        "ra5hdec-29.75": ((50, 100), (-30.5, -29)),
        "ra5hdec-33.25": ((50, 100), (-34, -32.5)),
        "ra5hdec-36.75": ((50, 100), (-37.5, -36)),
        "ra5hdec-40.25": ((50, 100), (-41, -39.5)),
        "spt3g-summera": ((50, 100), (-62, -29)),
        "ra1h40dec-29.75": ((0, 50), (-30.5, -29)),
        "ra1h40dec-33.25": ((0, 50), (-34, -32.5)),
        "ra1h40dec-36.75": ((0, 50), (-37.5, -36)),
        "ra1h40dec-40.25": ((0, 50), (-41, -39.5)),
        "spt3g-summerb": ((0, 50), (-41, -29)),
        "ra12h30dec-29.75": ((150, 225), (-30.5, -29)),
        "ra12h30dec-33.25": ((150, 225), (-34, -32.5)),
        "ra12h30dec-36.75": ((150, 225), (-37.5, -36)),
        "ra12h30dec-40.25": ((150, 225), (-41, -39.5)),
        "spt3g-summerc": ((150, 225), (-41, -29)),
        "j1745-290": ((266.1067, 266.7267), (-29.0028, -29.0128)),
        "spt3g-sgrastar": ((266.1067, 266.7267), (-29.0028, -29.0128)),
        "ra17h45dec-29": ((259.25, 273.25), (-30.5, -27.5)),
        "ra17h37dec-32": ((257.25, 271.25), (-33.5, -30.5)),
        "ra17h27dec-35": ((254.75, 268.75), (-36.5, -33.5)),
        "spt3g-galaxy": ((254.75, 273.25), (-36.5, -27.5)),
        "ra4h03dec-21": ((14.5, 107), (-21.5, -20.5)),
        "ra4h03dec-23": ((14.5, 107), (-23.5, -22.5)),
        "ra4h03dec-25": ((14.5, 107), (-25.5, -24.5)),
        "ra4h03dec-27": ((14.5, 107), (-27.5, -26.5)),
        "spt3g-widea": ((14.5, 107), (-27.5, -20.5)),
        "ra21h53dec-21": ((-78, 14.5), (-21.5, -20.5)),
        "ra21h53dec-23": ((-78, 14.5), (-23.5, -22.5)),
        "ra21h53dec-25": ((-78, 14.5), (-25.5, -24.5)),
        "ra21h53dec-27": ((-78, 14.5), (-27.5, -26.5)),
        "spt3g-wideb": ((-78, 14.5), (-27.5, -20.5)),
        "ra12hdec-21": ((120, 240), (-21.5, -20.5)),
        "ra12hdec-23": ((120, 240), (-23.5, -22.5)),
        "ra12hdec-25": ((120, 240), (-25.5, -24.5)),
        "ra12hdec-27": ((120, 240), (-27.5, -26.5)),
        "spt3g-widec": ((120, 240), (-27.5, -20.5)),
        "ra21h10dec-29.75": ((-85, 0), (-30.5, -29)),
        "ra21h10dec-33.25": ((-85, 0), (-34, -32.5)),
        "ra21h10dec-36.75": ((-85, 0), (-37.5, -36)),
        "ra21h10dec-40.25": ((-85, 0), (-41, -39.5)),
        "spt3g-wided": ((-85, 0), (-41, -29)),
        "ra18h40dec-45.5": ((-110, -50), (-48, -43)),
        "ra18h40dec-52.5": ((-110, -50), (-55, -50)),
        "ra18h40dec-59.5": ((-110, -50), (-62, -57)),
        "ra18h40dec-66.5": ((-110, -50), (-69, -64)),
        "spt3g-widee": ((-110, -50), (-69, -43)),
        "ra12h30dec-45": ((150, 225), (-47, -43)),
        "ra12h30dec-51": ((150, 225), (-53, -49)),
        "spt3g-widef": ((150, 225), (-53, -43)),
        "ra6h40dec-66.5": ((50, 150), (-69, -64)),
        "ra6h40dec-75": ((50, 150), (-79, -71)),
        "ra6h40dec-72.5": ((50, 150), (-74, -71)),
        "ra6h40dec-77.5": ((50, 150), (-79, -76)),
        "spt3g-wideg": ((50, 150), (-79, -64)),
        "ra23hdec-72.5": ((-80, 50), (-74, -71)),
        "ra23hdec-77.5": ((-80, 50), (-79, -76)),
        "spt3g-wideh": ((-80, 50), (-79, -71)),
        "ra14h20dec-72.5": ((150, 280), (-74, -71)),
        "ra14h20dec-77.5": ((150, 280), (-79, -76)),
        "spt3g-widei": ((150, 280), (-79, -71)),
        "ra4h04dec-48": ((55.5, 67), (-51.5, -45.5)),
        "spt3g-edfs": ((55.5, 67), (-51.5, -45.5)),
    }

    # Extra custom field for yearly -- added by FM
    # extents["yearly"] = extents["spt3g-winter"]

    ra, dec = extents[field]

    # padding in true degrees ra
    deg = core_G3Units_deg
    if sky_pad:
        ra_pad /= float(np.cos(np.max(np.abs(dec)) * deg / core_G3Units_rad))

    # apply padding
    ra = (
        (ra[0] * deg - ra_pad) % (360 * deg),
        (ra[1] * deg + ra_pad) % (360 * deg),
    )
    dec = (dec[0] * deg - dec_pad, dec[1] * deg + dec_pad)

    # IMPORTANT: Added for cutter -- we convert back to degrees
    ra = (ra[0]/deg, ra[1]/deg)
    dec = (dec[0]/deg, dec[1]/deg)

    return ra, dec


def get_boresight_field_extent(field):
    """
    Return the extent over which the telescope boresight moves in order to
    cover the field.  This includes padding in both RA and dec relative to the
    uniform coverage region returned by `get_field_extent`.

    Arguments
    ---------
    field : str
        Name of the field

    Returns
    -------
    ra_range : 2-tuple of floats
        Tuple of (min, max), the boresight extent in ra, in G3Units.
    dec_range : 2-tuple of floats, in G3Units
        Tuple of (min, max), the boresight extent in dec, in G3Units.
    """
    season = get_field_season(field).replace("spt3g-", "")
    dec_pad = (0.5 if season in ["widea", "wideb", "widec"] else 1.0) * core_G3Units_deg

    return get_field_extent(
        field, ra_pad=1.5 * core_G3Units_deg, dec_pad=dec_pad, sky_pad=False
    )


def get_field_center(
    field,
    ra_pad=3 * core_G3Units_deg,
    dec_pad=2 * core_G3Units_deg,
    sky_pad=True,
    extent=None,
    return_dims=False,
):
    """
    Return the center ra and dec for a given observing field.

    Arguments
    ---------
    field : str
        Name of the field
    ra_pad : float
        Padding to apply to extent in ra, in G3Units.  Values > 0 extend the
        field outward, and values < 0 shrink the field inward.  Set to 0 to
        select the uniformly weighted coverage region of the field.  If
        ``sky_pad`` is True, the padding is determined in terms of true degrees
        on the sky.
    dec_pad : float
        Padding to apply to extent in dec, in G3Units. Values > 0 extend the
        field outward, and values < 0 shrink the field inward.  Set to 0 to
        select the uniformly weighted coverage region of the field.
    sky_pad : bool
        If True, ``ra_pad`` is in terms of true degrees on the sky.  Otherwise,
        it is applied directly to the ra extent without correction.
    extent : list of 2-tuples
        Tuples of RA and dec ranges, as returned by ``get_field_extent``.  If
        supplied, then the ``field``, ``ra_pad`` and ``dec_pad`` arguments are
        ignored.
    return_dims : bool
        If True, also return width and height of the field given its extent.

    Returns
    -------
    ra0, dec0 : scalar
        Sky coordinates of the center of the observing field
    width, height : scalar
        Width and height of the observing field, if ``return_dims`` is True.
    """
    if extent is None:
        ra, dec = get_field_extent(field)
    else:
        ra, dec = extent

    ra0 = np.mean(ra)
    dec0 = np.mean(dec)

    if ra[0] > ra[1]:
        ra0 = (ra0 + 180 * core_G3Units_deg) % (360 * core_G3Units_deg)

    out = (ra0, dec0)

    if return_dims:
        width = np.diff(ra)[0]
        if width < 0:
            width += 360 * core_G3Units_deg
        height = np.diff(dec)[0]
        out += (width, height)

    return tuple(map(float, out))


def sample_field(field, n=1, **pad_kwargs):
    """
    Sample n points uniformly in the specified field.

    Arguments
    ---------
    field : str
        Name of the field or subfield to sample.
    n : int
        Number of points to generate.
    pad_kwargs : keyword arguments (ra_pad, dec_pad, sky_pad)
        Padding arguments to pass to ``get_field_extent``

    Returns
    -------
    random_ra : scalar or numpy array
        RA of uniformly sampled point(s)
    random_dec : scalar or numpy array
        Dec of uniformly sampled point(s)
    """

    (lo_ra, hi_ra), (lo_dec, hi_dec) = get_field_extent(field, **pad_kwargs)
    if lo_ra > hi_ra:
        lo_ra -= 360*core_G3Units_deg
    if lo_dec > hi_dec:
        print("dec bounds must be in increasing order")
    rand_ra = np.random.uniform(lo_ra, hi_ra, size=n)
    rand_dec = np.arcsin(np.random.uniform(np.sin(lo_dec), np.sin(hi_dec), size=n))
    if n == 1:
        return (rand_ra[0], rand_dec[0])
    return (rand_ra, rand_dec)
