/*
    luna.c  -  by Stephen F. Booth - 2023-10-20

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program calculates the Moon's appearance for a given date and time,
    or the computer's current date and time if none is given.
*/

#include <stdio.h>
#include <math.h>
#include "astro_demo_common.h"

#if 0

static astro_angle_result_t AngleError(astro_status_t status)
{
    astro_angle_result_t result;
    result.status = status;
    result.angle = NAN;
    return result;
}

/**
 * @brief Calculates the parallactic angle `q` for a body.
 *
 * @param body
 *     The celestial body to be observed. Not allowed to be `BODY_EARTH`.
 * @param time
 *     The date and time at which the observation takes place.
 * @param observer
 *     A location on or near the surface of the Earth.
 * @param aberration
 *     Selects whether or not to correct for aberration.
 * @param refraction
 *     Selects whether to correct for atmospheric refraction, and if so, which model to use.
 *     The recommended value for most uses is `REFRACTION_NORMAL`.
 *     See function remarks for more details.
 *
 * @return
 *     The body's parallactic angle.
 */
static astro_angle_result_t Astronomy_ParallacticAngle(astro_body_t body, astro_time_t *time, astro_observer_t observer, astro_aberration_t aberration, astro_refraction_t refraction)
{
    astro_equatorial_t body_equator_of_date;
//    astro_equatorial_t body_equator_j2000;
    astro_func_result_t body_hour_angle;
    double H, φ, δ;
    double q;
    astro_angle_result_t result;

    body_equator_of_date = Astronomy_Equator(body, time, observer, EQUATOR_OF_DATE, aberration);
    if (body_equator_of_date.status != ASTRO_SUCCESS)
        return AngleError(body_equator_of_date.status);

//    body_equator_j2000 = Astronomy_Equator(body, &time, observer, EQUATOR_J2000, aberration);
//    if (body_equator_j2000.status != ASTRO_SUCCESS)
//        return AngleError(body_equator_j2000.status);

    /*
     Calculate the parallactic angle (q) using
     Meeus equation 14.1 (p. 98).
     */
    body_hour_angle = Astronomy_HourAngle(body, time, observer);
    if (body_hour_angle.status != ASTRO_SUCCESS)
        return AngleError(body_hour_angle.status);

    H = body_hour_angle.value * HOUR2RAD;
    φ = observer.latitude * DEG2RAD;
    δ = body_equator_of_date.dec * DEG2RAD;
//    δ = body_equator_j2000.dec * DEG2RAD;

    q = atan2(sin(H), tan(φ) * cos(δ) - sin(δ) * cos(H));

    result.angle = q * RAD2DEG;
    result.status = ASTRO_SUCCESS;

    return result;
}

/**
 * @brief Calculates the position angle `χ` of a body's bright limb.
 *
 * The angle `χ` is the position angle of the midpoint of the illuminated limb
 * of the body reckoned eastward from the North Point of the disk (not from the axis
 * of rotation of the globe). The position angles of the cusps are χ ± 90 °.
 *
 * @param body
 *     The celestial body to be observed. Not allowed to be `BODY_EARTH`.
 * @param time
 *     The date and time at which the observation takes place.
 * @param observer
 *     A location on or near the surface of the Earth.
 * @param aberration
 *     Selects whether or not to correct for aberration.
 * @param refraction
 *     Selects whether to correct for atmospheric refraction, and if so, which model to use.
 *     The recommended value for most uses is `REFRACTION_NORMAL`.
 *     See function remarks for more details.
 *
 * @return
 *     The body's position angle.
 */
static astro_angle_result_t Astronomy_PositionAngle(astro_body_t body, astro_time_t *time, astro_observer_t observer, astro_aberration_t aberration, astro_refraction_t refraction)
{
    astro_equatorial_t body_equator_of_date;
    double α;
    astro_equatorial_t sun_equator_of_date;
    astro_equatorial_t sun_equator_j2000;
    double δ0, α0;
    double χ;
    double body_position_angle;
    astro_angle_result_t result;

    /*
     Calculate the position angle of the body's bright limb (χ) using
     Meeus equation 48.5 (p. 346).
     */
    body_equator_of_date = Astronomy_Equator(body, time, observer, EQUATOR_OF_DATE, aberration);
    if (body_equator_of_date.status != ASTRO_SUCCESS)
        return AngleError(body_equator_of_date.status);

    α = body_equator_of_date.ra * HOUR2RAD;

    sun_equator_of_date = Astronomy_Equator(BODY_SUN, time, observer, EQUATOR_OF_DATE, aberration);
    if (sun_equator_of_date.status != ASTRO_SUCCESS)
        return AngleError(sun_equator_of_date.status);

    sun_equator_j2000 = Astronomy_Equator(BODY_SUN, time, observer, EQUATOR_J2000, aberration);
    if (sun_equator_j2000.status != ASTRO_SUCCESS)
        return AngleError(sun_equator_j2000.status);

    δ0 = sun_equator_of_date.dec * DEG2RAD;
    α0 = sun_equator_of_date.ra * HOUR2RAD;
//    δ0 = sun_equator_j2000.dec * DEG2RAD;
//    α0 = sun_equator_j2000.ra * HOUR2RAD;

    χ = atan2(cos(δ0) * sin(α0 - α), sin(δ0) * cos(α) - cos(δ0) * sin(α) * cos(α0 - α));

    body_position_angle = χ * RAD2DEG;
    while (body_position_angle < 0)
        body_position_angle += 360;

    result.angle = body_position_angle;
    result.status = ASTRO_SUCCESS;

    return result;
}

#endif

static const char *PhaseAngleName(double eclipticPhaseAngle)
{
    long int phase = lround(eclipticPhaseAngle);
    if (phase >= 0 && phase <= 5)
        return "New";
    else if (phase >= 6 && phase <= 84)
        return "Waxing Crescent";
    else if (phase >= 85 && phase <= 95)
        return "First Quarter";
    else if (phase >= 96 && phase <= 174)
        return "Waxing Gibbous";
    else if (phase >= 175 && phase <= 185)
        return "Full";
    else if (phase >= 186 && phase <= 264)
        return "Waning Gibbous";
    else if (phase >= 265 && phase <= 275)
        return "Third Quarter";
    else if (phase >= 276 && phase <= 354)
        return "Waning Crescent";
    else if (phase >= 355 && phase <= 360)
        return "New";
    else
        return "INVALID ECLIPTIC PHASE ANGLE";
}

int main(int argc, const char *argv[])
{
    int error;
    astro_observer_t observer;
    astro_time_t time;
    astro_angle_result_t moon_phase;
    astro_illum_t moon_illumination;
    astro_equatorial_t moon_equator_of_date;
    astro_equatorial_t moon_equator_j2000;
    astro_horizon_t moon_horizontal_coordinates;
    astro_func_result_t moon_hour_angle;
    double H, φ, δ, α;
    double q;
    double moon_parallactic_angle;
//    astro_equatorial_t sun_equator_of_date;
    astro_equatorial_t sun_equator_j2000;
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

    puts("");

    /*
       Calculate the Moon's J2000 coordinates for the right ascension and declination.
     */
	moon_equator_j2000 = Astronomy_Equator(BODY_MOON, &time, observer, EQUATOR_J2000, ABERRATION);
    if (moon_equator_j2000.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "ERROR: Astronomy_Equator returned status %d trying to get J2000 coordinates.\n", moon_equator_j2000.status);
        return 1;
    }

    printf("Right Ascension = %.2lf hours\n", moon_equator_j2000.ra);
    printf("Declination = %.2lf degrees\n", moon_equator_j2000.dec);

    /*
       Calculate the Moon's horizontal coordinates for the azimuth and altitude.
     */
	moon_equator_of_date = Astronomy_Equator(BODY_MOON, &time, observer, EQUATOR_OF_DATE, ABERRATION);
    if (moon_equator_of_date.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "ERROR: Astronomy_Equator returned status %d trying to get coordinates of date.\n", moon_equator_of_date.status);
        return 1;
    }

    moon_horizontal_coordinates = Astronomy_Horizon(&time, observer, moon_equator_of_date.ra, moon_equator_of_date.dec, REFRACTION_NORMAL);

    printf("Azimuth = %.2lf degrees\n", moon_horizontal_coordinates.azimuth);
    printf("Altitude = %.2lf degrees\n", moon_horizontal_coordinates.altitude);

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
//    δ = moon_equator_of_date.dec * DEG2RAD;
//    α = moon_equator_of_date.ra * HOUR2RAD;
    δ = moon_equator_j2000.dec * DEG2RAD;
    α = moon_equator_j2000.ra * HOUR2RAD;

    q = atan2(sin(H), tan(φ) * cos(δ) - sin(δ) * cos(H));

    moon_parallactic_angle = q * RAD2DEG;

    /*
       Calculate the position angle of the Moon's bright limb (χ) using
       Meeus equation 48.5 (p. 346).

       The angle χ is the position angle of the midpoint of the illuminated limb
       of the Moon reckoned eastward from the North Point of the disk (not from the axis
       of rotation of the lunar globe). The position angles of the cusps are χ ± 90 °.
    */
//    sun_equator_of_date = Astronomy_Equator(BODY_SUN, &time, observer, EQUATOR_OF_DATE, ABERRATION);
//    if (sun_equator_of_date.status != ASTRO_SUCCESS)
//    {
//        printf("Astronomy_Equator error %d\n", sun_equator_of_date.status);
//        return 1;
//    }

    sun_equator_j2000 = Astronomy_Equator(BODY_SUN, &time, observer, EQUATOR_J2000, ABERRATION);
    if (sun_equator_j2000.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "ERROR: Astronomy_Equator returned status %d trying to get J2000 coordinates.\n", sun_equator_j2000.status);
        return 1;
    }

//    δ0 = sun_equator_of_date.dec * DEG2RAD;
//    α0 = sun_equator_of_date.ra * HOUR2RAD;
    δ0 = sun_equator_j2000.dec * DEG2RAD;
    α0 = sun_equator_j2000.ra * HOUR2RAD;

    χ = atan2(cos(δ0) * sin(α0 - α), sin(δ0) * cos(α) - cos(δ0) * sin(α) * cos(α0 - α));

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
