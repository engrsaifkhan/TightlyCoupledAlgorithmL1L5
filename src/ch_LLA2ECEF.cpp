#include "ch_LLA2ECEF.h"
#include <math.h>

void ch_LLA2ECEF(double lat, double lon, double h, double ECEF_XYZ[3])
{
    const double R0 = 6378137.0;          // WGS84 equatorial radius (m)
    const double e  = 0.0818191908425;    // WGS84 eccentricity

    const double sinLat = sin(lat);
    const double cosLat = cos(lat);
    const double sinLon = sin(lon);
    const double cosLon = cos(lon);

    const double RE = R0 / sqrt(1.0 - (e * e * sinLat * sinLat));

    ECEF_XYZ[0] = (RE + h) * cosLat * cosLon;
    ECEF_XYZ[1] = (RE + h) * cosLat * sinLon;
    ECEF_XYZ[2] = ((1.0 - e * e) * RE + h) * sinLat;
}