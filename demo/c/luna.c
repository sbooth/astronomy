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

// MARK: -

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
 * The parallactic angle `q` is the angle between a body's zenith (the uppermost point of the
 * body in the sky as seen by an observer) and the Northern celestial pole.
 * `q` is the angle between a body's vertical and its hour circle.
 *
 * @param body
 *     The celestial body to be observed. Not allowed to be `BODY_EARTH`.
 * @param time
 *     The date and time at which the observation takes place.
 * @param observer
 *     A location on or near the surface of the Earth.
 * @param aberration
 *     Selects whether or not to correct for aberration.
 *
 * @return
 *     The body's parallactic angle.
 */
astro_angle_result_t Astronomy_ParallacticAngle(astro_body_t body, astro_time_t *time, astro_observer_t observer, astro_aberration_t aberration)
{
    astro_equatorial_t body_equator_of_date;
    astro_func_result_t body_hour_angle;
    double H, φ, δ;
    double q;
    astro_angle_result_t result;

    if (time == NULL)
        return AngleError(ASTRO_INVALID_PARAMETER);

    body_equator_of_date = Astronomy_Equator(body, time, observer, EQUATOR_OF_DATE, aberration);
    if (body_equator_of_date.status != ASTRO_SUCCESS)
        return AngleError(body_equator_of_date.status);

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

    q = atan2(sin(H), tan(φ) * cos(δ) - sin(δ) * cos(H));

    result.angle = q * RAD2DEG;
    result.status = ASTRO_SUCCESS;

    return result;
}

/**
 * @brief Calculates the position angle `χ` of a body's bright limb.
 *
 * The angle `χ` is the position angle of the midpoint of the illuminated limb
 * of the body reckoned eastward from the North point of the disk (not from the axis
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
 *
 * @return
 *     The body's bright limb position angle.
 */
astro_angle_result_t Astronomy_PositionAngle(astro_body_t body, astro_time_t *time, astro_observer_t observer, astro_aberration_t aberration)
{
    astro_equatorial_t body_equator_of_date;
    double δ, α;
    astro_equatorial_t sun_equator_of_date;
    double δ0, α0;
    double χ;
    double body_position_angle;
    astro_angle_result_t result;

    if (time == NULL)
        return AngleError(ASTRO_INVALID_PARAMETER);

    /*
     Calculate the position angle of the body's bright limb (χ) using
     Meeus equation 48.5 (p. 346).
     */
    body_equator_of_date = Astronomy_Equator(body, time, observer, EQUATOR_OF_DATE, aberration);
    if (body_equator_of_date.status != ASTRO_SUCCESS)
        return AngleError(body_equator_of_date.status);

    δ = body_equator_of_date.dec * DEG2RAD;
    α = body_equator_of_date.ra * HOUR2RAD;

    sun_equator_of_date = Astronomy_Equator(BODY_SUN, time, observer, EQUATOR_OF_DATE, aberration);
    if (sun_equator_of_date.status != ASTRO_SUCCESS)
        return AngleError(sun_equator_of_date.status);

    δ0 = sun_equator_of_date.dec * DEG2RAD;
    α0 = sun_equator_of_date.ra * HOUR2RAD;

    χ = atan2(cos(δ0) * sin(α0 - α), sin(δ0) * cos(δ) - cos(δ0) * sin(δ) * cos(α0 - α));

    body_position_angle = χ * RAD2DEG;
    while (body_position_angle < 0)
        body_position_angle += 360;

    result.angle = body_position_angle;
    result.status = ASTRO_SUCCESS;

    return result;
}

// MARK: -

static const char *PhaseAngleName(double eclipticPhaseAngle)
{
    assert(eclipticPhaseAngle >= 0);
    assert(eclipticPhaseAngle <= 360);

    int phase = (int)floor(eclipticPhaseAngle / (360 / 8));
    switch (phase)
    {
        case 0:
        case 8:
            return "🌑 New";
        case 1:
            return "🌒 Waxing Crescent";
        case 2:
            return "🌓 First Quarter";
        case 3:
            return "🌔 Waxing Gibbous";
        case 4:
            return "🌕 Full";
        case 5:
            return "🌖 Waning Gibbous";
        case 6:
            return "🌗 Third Quarter";
        case 7:
            return "🌘 Waning Crescent";
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
    astro_angle_result_t moon_parallactic_angle;
    astro_angle_result_t moon_position_angle;
    double moon_zenith_angle;

    error = ParseArgs(argc, argv, &observer, &time);
    if (error)
        return error;

    puts("           ┏━━━━━━━━━━━━━━━━━━━┓");
    puts("           ┃     The Moon      ┃");
    puts("           ┗━━━━━━━━━┳━━━━━━━━━┛");

    printf("%-*s ┃ %s %.2lf %s %.2lf degrees\n", 20, "Observer position", observer.latitude >= 0 ? "N" : "S", fabs(observer.latitude), observer.longitude >= 0 ? "E" : "W", fabs(observer.longitude));

    printf("%-*s ┃ ", 20, "UTC date");
    PrintTime(time);
    printf("\n");

    puts("                ━━━━━╋━━━━━");

    /*
        Calculate the Moon's ecliptic phase angle,
        which ranges from 0 to 360 degrees.
    */
    moon_phase = Astronomy_MoonPhase(time);
    if (moon_phase.status != ASTRO_SUCCESS)
    {
        printf("Astronomy_MoonPhase error %d\n", moon_phase.status);
        return 1;
    }

    printf("%-*s ┃ %0.2lf degrees\n", 20, "Ecliptic phase angle", moon_phase.angle);
    printf("%-*s ┃ %s\n", 20, "Appearance", PhaseAngleName(moon_phase.angle));

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

    printf("%-*s ┃ %0.2lf%%\n", 20, "Illuminated fraction", 100.0 * moon_illumination.phase_fraction);
    printf("%-*s ┃ %0.2lf\n", 20, "Magnitude", moon_illumination.mag);

    puts("                ━━━━━╋━━━━━");

    /*
        Calculate the Moon's horizontal coordinates for the azimuth and altitude.
    */
    moon_equator_of_date = Astronomy_Equator(BODY_MOON, &time, observer, EQUATOR_OF_DATE, ABERRATION);
    if (moon_equator_of_date.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "ERROR: Astronomy_Equator returned status %d trying to get coordinates of date.\n", moon_equator_of_date.status);
        return 1;
    }

    printf("%-*s ┃ %.2lf hours\n", 20, "Right ascension", moon_equator_of_date.ra);
    printf("%-*s ┃ %.2lf degrees\n", 20, "Declination", moon_equator_of_date.dec);

    moon_horizontal_coordinates = Astronomy_Horizon(&time, observer, moon_equator_of_date.ra, moon_equator_of_date.dec, REFRACTION_NORMAL);

    printf("%-*s ┃ %.2lf degrees\n", 20, "Azimuth", moon_horizontal_coordinates.azimuth);
    printf("%-*s ┃ %.2lf degrees\n", 20, "Altitude", moon_horizontal_coordinates.altitude);

    puts("                ━━━━━╋━━━━━");

    moon_libration = Astronomy_Libration(time);

    printf("%-*s ┃ %ld kilometers\n", 20, "Distance", lround(moon_libration.dist_km));
    printf("%-*s ┃ %.2lf degrees\n", 20, "Apparent diameter", moon_libration.diam_deg);

   puts("                ━━━━━╋━━━━━");

    /*
        Calculate the Moon's parallactic angle (q)
    */
    moon_parallactic_angle = Astronomy_ParallacticAngle(BODY_MOON, &time, observer, ABERRATION);
    if (moon_parallactic_angle.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "ERROR: Astronomy_ParallacticAngle returned status %d trying to get Moon's parallactic angle.\n", moon_parallactic_angle.status);
        return 1;
    }

    /*
        Calculate the position angle of the Moon's bright limb (χ)
    */
	moon_position_angle = Astronomy_PositionAngle(BODY_MOON, &time, observer, ABERRATION);
    if (moon_position_angle.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "ERROR: Astronomy_PositionAngle returned status %d trying to get Moon's position angle.\n", moon_position_angle.status);
        return 1;
    }

    /*
        The angle χ is not measured from the direction of the observer's zenith.
        The zenith angle of the bright limb is χ - q.
    */
    moon_zenith_angle = moon_position_angle.angle - moon_parallactic_angle.angle;

    printf("%-*s ┃ %0.2lf degrees\n", 20, "Parallactic angle", moon_parallactic_angle.angle);
    printf("%-*s ┃ %0.2lf degrees\n", 20, "Position angle", moon_position_angle.angle);
    printf("%-*s ┃ %0.2lf degrees\n", 20, "Zenith angle", moon_zenith_angle);

    puts("                ━━━━━┻━━━━━");

    return 0;
}
