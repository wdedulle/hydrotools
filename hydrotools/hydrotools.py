"""Main module."""

from math import (
    atan2,
    cos,
    sin,
    degrees,
    modf,
    hypot,
    isnan,
    radians,
    atan,
    tan,
    sqrt,
    asin,
)
import re
import json


# Determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs. This function
# returns True or False.  The algorithm is called
# "Ray Casting Method".
def point_in_poly(x, y, poly):
    if (x, y) in poly:
        return True
    for I in range(len(poly)):
        p1 = None
        p2 = None
        if I == 0:
            p1 = poly[0]
            p2 = poly[1]
        else:
            p1 = poly[I - 1]
            p2 = poly[I]
        if (
            p1[1] == p2[1]
            and p1[1] == y
            and x > min(p1[0], p2[0])
            and x < max(p1[0], p2[0])
        ):
            return True

        n = len(poly)
        inside = False

        p1x, p1y = poly[0]
        for i in range(n + 1):
            p2x, p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y

    return inside


# Calculate line direction
def bearing(x1, y1, x2, y2):
    angle = atan2(
        cos(y1) * sin(y2) - sin(y1) * cos(y2) * cos(x2 - x1), sin(x2 - x1) * cos(y2)
    )
    return (degrees(angle) + 360) % 360


# decimals to degrees
def dd2dms(lat, lon):
    latf, latn = modf(lat)
    lonf, lonn = modf(lon)
    latd = int(latn)
    latm = int(latf * 60)
    lats = (lat - latd - latm / 60) * 3600.00
    lond = int(lonn)
    lonm = int(lonf * 60)
    lons = (lon - lond - lonm / 60) * 3600.00
    compass = {"lat": ("N", "S"), "lon": ("E", "W")}

    lat_compass = compass["lat"][0 if latd >= 0 else 1]
    lon_compass = compass["lon"][0 if lond >= 0 else 1]

    return f"{abs(latd)}°{abs(latm):.2f}'{abs(lats):.2f}\" {lat_compass},{abs(lond)}°{abs(lonm):.2f}'{abs(lons):.2f}\" {lon_compass}"


# degrees to decimals
def dms2dd(lat, lon):
    latd, latm, lats, lat_compass = re.split(r"[^\d\.A-Z]+", lat)
    lond, lonm, lons, lon_compass = re.split(r"[^\d\.A-Z]+", lon)
    latdd = float(latd) + float(latm) / 60 + float(lats) / 3600
    londd = float(lond) + float(lonm) / 60 + float(lons) / 3600
    if lat_compass == "S":
        latdd *= -1
    if lon_compass == "W":
        londd *= -1

    return (latdd, londd)


def calculate_flat_perimeter(geojson_str):
    """
    Berekent de omtrek van een GeoJSON-polygoon op een plat vlak
    met behulp van math.hypot.
    """
    # Laad de GeoJSON string naar een Python dictionary
    data = json.loads(geojson_str)

    # Controleer of het om een Feature of direct om een Geometry gaat
    if data.get("type") == "Feature":
        geometry = data.get("geometry", {})
    else:
        geometry = data

    if geometry.get("type") != "Polygon":
        raise ValueError("De meegegeven GeoJSON is geen geldige 'Polygon'.")

    # Haal de buitenste ring op (index 0 van de coördinaten)
    exterior_ring = geometry["coordinates"][0]

    perimeter = 0.0

    for i in range(len(exterior_ring) - 1):
        x1, y1 = exterior_ring[i]
        x2, y2 = exterior_ring[i + 1]

        # math.hypot(dx, dy) berekent de wortel van (dx² + dy²)
        perimeter += hypot(x2 - x1, y2 - y1)

    return perimeter


def distance_pythagoras(x1, y1, x2, y2):
    x_dist = radians(x1 - x2)
    y_dist = radians(y1 - y2)

    dist_sqr = x_dist**2 + y_dist**2
    dist_rad = sqrt(dist_sqr)
    dist = dist_rad * 6371
    return dist


def distance_halversine(x1, y1, x2, y2):
    x_dist = radians(x1 - x2)
    y_dist = radians(y1 - y2)
    y1_rad = radians(y1)
    y2_rad = radians(y2)

    a = sin(y_dist / 2) ** 2 + cos(y1_rad) * cos(y2_rad) * sin(x_dist / 2) ** 2
    c = 2 * asin(sqrt(a))
    dist = c * 6371
    return dist


def distance_vincenty(x1, y1, x2, y2):
    distance = None
    a = 6378137
    f = 1 / 298.257223563
    b = abs((f * a) - a)

    L = radians(x2 - x1)
    U1 = atan((1 - f) * tan(radians(y1)))
    U2 = atan((1 - f) * tan(radians(y2)))
    sinU1 = sin(U1)
    cosU1 = cos(U1)
    sinU2 = sin(U2)
    cosU2 = cos(U2)
    lam = L
    for i in range(100):
        sinLam = sin(lam)
        cosLam = cos(lam)
        sinSigma = sqrt(
            (cosU2 * sinLam) ** 2 + (cosU1 * sinU2 - sinU1 * cosU2 * cosLam) ** 2
        )

        if sinSigma == 0:
            distance = 0
            break
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLam
        sigma = atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLam / sinSigma
        cosSqAlpha = 1 - sinAlpha**2
        cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha

        if isnan(cos2SigmaM):
            cos2SigmaM = 0
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        LP = lam
        lam = L + (1 - C) * f * sinAlpha * (
            sigma
            + C
            * sinSigma
            * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM))
        )
        if not abs(lam - LP) > 1e-12:
            break

    uSq = cosSqAlpha * (a**2 - b**2) / b**2
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))

    deltaSigma = (
        B
        * sinSigma
        * (
            cos2SigmaM
            + B
            / 4
            * (
                cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)
                - B
                / 6
                * cos2SigmaM
                * (-3 + 4 * sinSigma * sinSigma)
                * (-3 + 4 * cos2SigmaM * cos2SigmaM)
            )
        )
    )
    s = b * A * (sigma - deltaSigma)
    distance = s

    return distance / 1000
