/*
    luna.c  -  by Stephen F. Booth - 2023-10-20

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates the Moon's appearance for a given date and time,
    or the computer's current date and time if none is given.
*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "astro_demo_common.h"

static const char *PhaseAngleName(double eclipticPhaseAngle)
{
    assert(eclipticPhaseAngle >= 0);
    assert(eclipticPhaseAngle <= 360);

    int phase = (int)floor(eclipticPhaseAngle / (360 / 8));
    switch (phase)
    {
        case 0:
        case 8:
            return "New";
        case 1:
            return "Waxing Crescent";
        case 2:
            return "First Quarter";
        case 3:
            return "Waxing Gibbous";
        case 4:
            return "Full";
        case 5:
            return "Waning Gibbous";
        case 6:
            return "Third Quarter";
        case 7:
            return "Waning Crescent";
        default:
            return "INVALID ECLIPTIC PHASE ANGLE";
    }
}

int main(int argc, const char *argv[])
{
    int error;
    astro_observer_t observer;
    astro_time_t time;
    astro_angle_result_t moon_phase;
    astro_illum_t moon_illumination;
    astro_equatorial_t moon_equator_of_date;
    astro_horizon_t moon_horizontal_coordinates;
    astro_libration_t moon_libration;
    astro_func_result_t moon_hour_angle;
    double H, φ, δ, α;
    double q;
    double moon_parallactic_angle;
    astro_equatorial_t sun_equator_of_date;
    double δ0, α0;
    double χ;
    double moon_position_angle;
    double moon_zenith_angle;

    error = ParseArgs(argc, argv, &observer, &time);
    if (error)
        return error;

    puts("+==============+");
    puts("|     Luna     |");
    puts("+==============+");

    puts("");

    printf("Observer position = %s %.2lf %s %.2lf degrees\n", observer.latitude >= 0 ? "N" : "S", fabs(observer.latitude), observer.longitude >= 0 ? "E" : "W", fabs(observer.longitude));

    printf("UTC date = ");
    PrintTime(time);
    printf("\n");

    puts("");

    /*
        Calculate the Moon's ecliptic phase angle,
        which ranges from 0 to 360 degrees.

          0 = new moon,
         90 = first quarter,
        180 = full moon,
        270 = third quarter.
    */
    moon_phase = Astronomy_MoonPhase(time);
    if (moon_phase.status != ASTRO_SUCCESS)
    {
        printf("Astronomy_MoonPhase error %d\n", moon_phase.status);
        return 1;
    }

    printf("Ecliptic phase angle = %0.3lf degrees\n", moon_phase.angle);
    printf("Appearance = %s\n", PhaseAngleName(moon_phase.angle));

    /*
        Calculate the percentage of the Moon's disc that is illuminated
        from the Earth's point of view.
    */
    moon_illumination = Astronomy_Illumination(BODY_MOON, time);
    if (moon_illumination.status != ASTRO_SUCCESS)
    {
        printf("Astronomy_Illumination error %d\n", moon_illumination.status);
        return 1;
    }

    printf("Illuminated fraction = %0.2lf%%\n", 100.0 * moon_illumination.phase_fraction);
    printf("Magnitude = %0.2lf\n", moon_illumination.mag);

    puts("");

    /*
        Calculate the Moon's horizontal coordinates for the azimuth and altitude.
    */
    moon_equator_of_date = Astronomy_Equator(BODY_MOON, &time, observer, EQUATOR_OF_DATE, ABERRATION);
    if (moon_equator_of_date.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "ERROR: Astronomy_Equator returned status %d trying to get coordinates of date.\n", moon_equator_of_date.status);
        return 1;
    }

    printf("Right Ascension = %.2lf hours\n", moon_equator_of_date.ra);
    printf("Declination = %.2lf degrees\n", moon_equator_of_date.dec);

    moon_horizontal_coordinates = Astronomy_Horizon(&time, observer, moon_equator_of_date.ra, moon_equator_of_date.dec, REFRACTION_NORMAL);

    printf("Azimuth = %.2lf degrees\n", moon_horizontal_coordinates.azimuth);
    printf("Altitude = %.2lf degrees\n", moon_horizontal_coordinates.altitude);

    puts("");

    moon_libration = Astronomy_Libration(time);

    printf("Distance = %.2lf kilometers\n", moon_libration.dist_km);
    printf("Apparent Diameter = %.2lf degrees\n", moon_libration.diam_deg);

    puts("");

    /*
        Calculate the Moon's parallactic angle (q) using
        Meeus equation 14.1 (p. 98).
    */
    moon_hour_angle = Astronomy_HourAngle(BODY_MOON, &time, observer);
    if (moon_hour_angle.status != ASTRO_SUCCESS)
    {
        printf("Astronomy_HourAngle error %d\n", moon_hour_angle.status);
        return 1;
    }

    H = moon_hour_angle.value * HOUR2RAD;
    φ = observer.latitude * DEG2RAD;
    δ = moon_equator_of_date.dec * DEG2RAD;
    α = moon_equator_of_date.ra * HOUR2RAD;

    q = atan2(sin(H), tan(φ) * cos(δ) - sin(δ) * cos(H));

    moon_parallactic_angle = q * RAD2DEG;

    /*
        Calculate the position angle of the Moon's bright limb (χ) using
        Meeus equation 48.5 (p. 346).

        The angle χ is the position angle of the midpoint of the illuminated limb
        of the Moon reckoned eastward from the North Point of the disk (not from the axis
        of rotation of the lunar globe). The position angles of the cusps are χ ± 90 °.
    */
    sun_equator_of_date = Astronomy_Equator(BODY_SUN, &time, observer, EQUATOR_OF_DATE, ABERRATION);
    if (sun_equator_of_date.status != ASTRO_SUCCESS)
    {
        printf("Astronomy_Equator error %d\n", sun_equator_of_date.status);
        return 1;
    }

    δ0 = sun_equator_of_date.dec * DEG2RAD;
    α0 = sun_equator_of_date.ra * HOUR2RAD;

    χ = atan2(cos(δ0) * sin(α0 - α), sin(δ0) * cos(δ) - cos(δ0) * sin(δ) * cos(α0 - α));

    moon_position_angle = χ * RAD2DEG;
    while (moon_position_angle < 0)
        moon_position_angle += 360;

    /*
        The angle χ is not measured from the direction of the observer's zenith.
        The zenith angle of the bright limb is χ - q.
    */
    moon_zenith_angle = moon_position_angle - moon_parallactic_angle;

    printf("Parallactic angle = %0.3lf degrees\n", moon_parallactic_angle);
    printf("Bright limb position angle = %0.3lf degrees\n", moon_position_angle);
    printf("Zenith angle = %0.3lf degrees\n", moon_zenith_angle);

    return 0;
}
