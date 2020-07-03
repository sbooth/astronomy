/*
    top2013.c  -  Don Cross  -  2020-06-30
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "codegen.h"
#include "top2013.h"

static const double dpi = 6.283185307179586476925287;

static int FormatTermLine(int lnum, char *line, size_t size, const top_term_t *term);


void TopInitModel(top_model_t *model)
{
    memset(model, 0, sizeof(top_model_t));
}


void TopFreeModel(top_model_t *model)
{
    int f, s;
    for (f=0; f < TOP_NCOORDS; ++f)
        for (s=0; s < model->formula[f].nseries_total; ++s)
            free(model->formula[f].series[s].terms);

    TopInitModel(model);
}


static int RoundingAdjustment(char original, char regen, int *diff)
{
    switch (original - regen)
    {
    case 0:
        *diff = 0;
        return 0;

    case -1:
    case +9:
        *diff = -1;
        return 0;

    case +1:
    case -9:
        *diff = +1;
        return 0;

    default:
        fprintf(stderr, "RoundingAdjustment: original=%c, regen=%c\n", original, regen);
        *diff = 0;
        return 1;
    }
}


int TopLoadModel(top_model_t *model, const char *filename, int planet)
{
    int error = 1;
    FILE *infile = NULL;
    int nterms_remaining, lnum, nscanned, check_planet, check_var, tpower;
    int count = 0;      /* total number of terms we found for this planet */
    top_formula_t *formula = NULL;
    top_series_t *series = NULL;
    top_term_t *term = NULL;
    char line[100];
    char gline[100];

    TopInitModel(model);
    model->planet = planet;

    infile = fopen(filename, "rt");
    if (infile == NULL)
        FAIL("TopLoadModel: cannot open file: %s\n", filename);

    lnum = 0;
    nterms_remaining = 0;
    check_planet = 0;
    check_var = 0;
    while (fgets(line, sizeof(line), infile))
    {
        ++lnum;
        if (nterms_remaining == 0)
        {
            nscanned = sscanf(line,
                " TOP2013ELL PLANET %d VARIABLE %d T**%d %d term(s)",
                &check_planet, &check_var, &tpower, &nterms_remaining);

            if (nscanned != 4)
                FAIL("TopLoadModel(%s line %d): invalid data format.\n", filename, lnum);

            if (check_var < 1 || check_var > TOP_NCOORDS)
                FAIL("TopLoadModel(%s line %d): invalid variable number %d\n", filename, lnum, check_var);

            --check_var;    /* convert one-based to zero-based indexing */

            if (tpower < 0 || tpower >= TOP_MAX_SERIES)
                FAIL("TopLoadModel(%s line %d): invalid power of t: %d\n", filename, lnum, tpower);

            if (check_planet == planet)
            {
                if (series != NULL)
                {
                    /* make sure the previous series is complete */
                    if (series->nterms_calc != series->nterms_total)
                        FAIL("TopLoadModel(%s line %d): previous series has %d terms; expected %d\n", filename, lnum, series->nterms_calc, series->nterms_total);
                }

                /* Allocate room for the new terms. */
                formula = &model->formula[check_var];
                formula->nseries_calc = formula->nseries_total = tpower + 1;
                series = &formula->series[tpower];
                series->nterms_total = nterms_remaining;
                series->nterms_calc = 0;
                series->terms = calloc(sizeof(top_term_t), nterms_remaining);
                if (series->terms == NULL)
                    FAIL("TopLoadModel(%s line %d): out of memory\n", filename, lnum);
            }
        }
        else
        {
            --nterms_remaining;

            if (check_planet == planet)
            {
                if (series == NULL)
                    FAIL("TopLoadModel(%s line %d): series == NULL\n", filename, lnum);

                if (series->nterms_calc >= series->nterms_total)
                    FAIL("TopLoadModel(%s line %d): too many terms\n", filename, lnum);

                term = &series->terms[series->nterms_calc++];

                if (strlen(line) < 61)
                    FAIL("TopLoadModel(%s line %d) line is too short.\n", filename, lnum);

                /* patch in 'e' to make numbers in scientific notation. */
                if (line[31] != ' ' || line[57] != ' ')
                    FAIL("TopLoadModel(%s line %d): expected spaces between mantissas and exponents.\n", filename, lnum);

                line[31] = line[57] = 'e';
                nscanned = sscanf(line, "%lf %lf %lf %lf", &term->k, &term->c, &term->s, &term->p);
                if (nscanned == 3)
                    term->p = 0.0;
                else if (nscanned != 4)
                    FAIL("TopLoadModel(%s line %d): invalid term data format.\n", filename, lnum);

                /*
                    Super ugly hack: the weird exponential notation used in TOP2013.dat
                    makes it really difficult to regenerate exactly in all cases.
                    So I generate each line back as text, and remember what I will have to do later
                    to reconstruct the input exactly.
                */
                CHECK(FormatTermLine(lnum, gline, sizeof(gline), term));
                CHECK(RoundingAdjustment(line[30], gline[30], &term->rc));
                CHECK(RoundingAdjustment(line[56], gline[56], &term->rs));

                /* Sanity check the rounding hack. */
                CHECK(FormatTermLine(lnum, gline, sizeof(gline), term));
                line[31] = line[57] = ' ';
                if (strcmp(line, gline))
                {
                    fprintf(stderr, "INPUT:%s", line);
                    fprintf(stderr, "REGEN:%s", gline);
                    FAIL("TopLoadModel(%s line %d): unable to reconstruct identical term line.\n", filename, lnum);
                }

                ++count;
            }
        }
    }

    if (nterms_remaining != 0)
        FAIL("TopLoadModel(%s): missing %d terms at the end.\n", filename, nterms_remaining);

    if (count == 0)
        FAIL("TopLoadModel(%s): could not find any terms for planet %d.\n", filename, planet);

    error = 0;
fail:
    if (infile != NULL) fclose(infile);
    if (error) TopFreeModel(model);
    return error;
}


static int AppendTrigCoeff(char *line, int lnum, double x, int rounding_adjust)
{
    int error = 1;
    int m, length, exponent;
    char polarity;
    char buffer[40];

    if (rounding_adjust < -1 || rounding_adjust > +1)
        FAIL("AppendTrigCoeff(%d): invalid rounding_adjust = %d\n", lnum, rounding_adjust);

    /*
        The data format for TOP2013 has a weird form of scientific notation.
        It always has a 0 to the left of the decimal point, which wastes a digit of precision.
        "-5.2026032025158849e+00"       <== C scientific notation
        "   -0.5202603202515885 +01"    <== output
    */

    snprintf(buffer, sizeof(buffer), "%23.16le", x);
    length = (int)strlen(buffer);
    if (length != 23)
        FAIL("AppendTrigCoeff(%d): output string '%s' has incorrect length %d.\n", lnum, buffer, length);

    if (buffer[19] != 'e')
        FAIL("AppendTrigCoeff(%d): expected 'e' at index 19 in string '%s'\n", lnum, buffer);

    buffer[19] = '\0';      /* truncate the string at the 'e' */

    if (1 != sscanf(buffer+20, "%d", &exponent))
        FAIL("AppendTrigCoeff(%d): cannot scan exponent from '%s'\n", lnum, buffer+20);

    ++exponent;
    if (x == 0.0)
        exponent = 0;
    if (exponent >= 0)
    {
        polarity = '+';
    }
    else
    {
        polarity = '-';
        exponent = -exponent;
    }

    /* Copy digits and shift decimal point */
    buffer[22] = '\0';
    for (m = 17; m >= 0 && buffer[m] != '.'; --m)
        buffer[m+4] = buffer[m];

    if (m != 2 || buffer[2] != '.')
        FAIL("AppendTrigCoeff(%d): decimal point is in the wrong place: '%s'\n", lnum, buffer);

    buffer[6] = buffer[1];
    buffer[5] = '.';
    buffer[4] = '0';
    buffer[3] = buffer[0];
    buffer[2] = ' ';
    buffer[1] = ' ';
    buffer[0] = ' ';

    if (rounding_adjust != 0)
    {
        for (m=21; m >= 0; --m)
        {
            if (buffer[m] != '.')
            {
                if (buffer[m] < '0' || buffer[m] > '9')
                    FAIL("AppendTrigCoeff(%d): rounding failure\n", lnum);
                buffer[m] += (char)rounding_adjust;
                if (buffer[m] < '0')
                    buffer[m] = '9';
                else if (buffer[m] > '9')
                    buffer[m] = '0';
                else
                    break;
            }
        }
    }

    sprintf(buffer+22, " %c%02d", polarity, exponent);
    length = (int)strlen(buffer);
    if (length != 26)
        FAIL("AppendTrigCoeff(%d): generated incorrect length %d in string '%s' for x=%lg\n", lnum, length, buffer, x);

    strcat(line, buffer);
    error = 0;
fail:
    return error;
}


static int FormatTermLine(int lnum, char *line, size_t size, const top_term_t *term)
{
    int error = 1;
    int length;

    snprintf(line, size, "%9.0lf", term->k);
    CHECK(AppendTrigCoeff(line, lnum, term->c, term->rc));
    CHECK(AppendTrigCoeff(line, lnum, term->s, term->rs));
    length = (int)strlen(line);
    if (61 != length)
        FAIL("FormatTermLine(%d): incorrect output line length = %d.\n", lnum, length);
    if (term->k != 0)
        snprintf(line + length, size - (size_t)length, " %11.6lf", term->p);
    strcat(line, "\n");
    error = 0;
fail:
    return error;
}


int TopWriteModel(const top_model_t *model, FILE *outfile)
{
    int error = 1;
    int f, s, t;
    int lnum = 0;
    char line[100];

    for (f=0; f < TOP_NCOORDS; ++f)
    {
        const top_formula_t *formula = &model->formula[f];
        for (s=0; s < formula->nseries_calc; ++s)
        {
            const top_series_t *series = &formula->series[s];

            if (series->nterms_calc == 0)
                continue;

            ++lnum;
            if (0 > fprintf(outfile, " TOP2013ELL    PLANET %d    VARIABLE %d    T**%02d %7d term(s)\n", model->planet, f+1, s, series->nterms_calc))
                FAIL("TopWriteModel(%d): error writing header record to output stream.\n", lnum);

            for (t=0; t < series->nterms_calc; ++t)
            {
                const top_term_t *term = &series->terms[t];
                ++lnum;
                CHECK(FormatTermLine(lnum, line, sizeof(line), term));
                if (0 > fprintf(outfile, "%s", line))
                    FAIL("TopWriteModel(%d): error writing term record to output stream.\n", lnum);
            }
        }
    }
    error = 0;
fail:
    return error;
}


int TopSaveModel(const top_model_t *model, const char *filename)
{
    int error;
    FILE *outfile = fopen(filename, "wt");
    if (outfile == NULL)
        FAIL("TopSaveModel: Cannot open output file: %s\n", filename);
    CHECK(TopWriteModel(model, outfile));
fail:
    if (outfile) fclose(outfile);
    return error;
}


void TopResetModel(top_model_t *model)
{
    int f, s;

    /* Undo any truncation. */

    for (f=0; f < TOP_NCOORDS; ++f)
    {
        top_formula_t *formula = &model->formula[f];
        formula->nseries_calc = formula->nseries_total;
        for (s=0; s < formula->nseries_total; ++s)
        {
            top_series_t *series = &formula->series[s];
            series->nterms_calc = series->nterms_total;
        }
    }
}


static int ContribCompare(const void *aptr, const void *bptr)
{
    const top_contrib_t *a = aptr;
    const top_contrib_t *b = bptr;

    if (a->magnitude < b->magnitude)
        return -1;

    if (a->magnitude > b->magnitude)
        return +1;

    /* tie-breaker: keep elements in reverse order as much as possible */

    if (a->s != b->s)
        return b->s - a->s;

    return b->t - a->t;
}


static int MakeContribList(const top_formula_t *formula, top_contrib_list_t *list, double millennia)
{
    int error;
    int nterms, s, t;
    double tpower;

    /* Count up the number of terms in the elliptical coordinate formula 'f'. */

    nterms = 0;
    for (s=0; s < formula->nseries_calc; ++s)
        nterms += formula->series[s].nterms_calc;

    /* Allocate space for the array. */
    list->array = calloc(nterms, sizeof(top_contrib_t));
    if (list->array == NULL)
        FAIL("MakeContribList: out of memory allocating %d elements.\n", nterms);

    /* Fill up the array. */
    list->nterms = nterms;
    nterms = 0;
    tpower = 1.0;
    millennia = fabs(millennia);
    for (s=0; s < formula->nseries_calc; ++s)
    {
        const top_series_t *series = &formula->series[s];
        for (t=0; t < series->nterms_calc; ++t)
        {
            top_contrib_t *contrib;
            const top_term_t *term = &series->terms[t];
            if (nterms >= list->nterms)
                FAIL("MakeContribList: internal error -- overflow %d items\n", nterms);
            contrib = &list->array[nterms++];
            contrib->s = s;
            contrib->t = t;
            contrib->magnitude = tpower * sqrt(term->c*term->c + term->s*term->s);
        }
        tpower *= millennia;
    }

    if (nterms != list->nterms)
        FAIL("MakeContribList: internal error: expected %d terms, found %d\n", list->nterms, nterms);

    /* Sort the list in ascending order of magnitude. */
    qsort(list->array, list->nterms, sizeof(list->array[0]), ContribCompare);

    error = 0;
fail:
    return error;
}


void TopInitContribMap(top_contrib_map_t *map)
{
    memset(map, 0, sizeof(top_contrib_map_t));
}


void TopFreeContribMap(top_contrib_map_t *map)
{
    int f;

    for (f=0; f < TOP_NCOORDS; ++f)
        free(map->list[f].array);

    TopInitContribMap(map);
}


int TopMakeContribMap(const top_model_t *model, top_contrib_map_t *map, double millennia)
{
    int error, f;

    TopInitContribMap(map);
    for (f=0; f < TOP_NCOORDS; ++f)
        CHECK(MakeContribList(&model->formula[f], &map->list[f], millennia));

fail:
    if (error) TopFreeContribMap(map);
    return error;
}


int TopCalcElliptical(const top_model_t *model, double tt, top_elliptical_t *ellip)
{
    /* Translated from: TOP2013.f */
    /* See: https://github.com/cosinekitty/ephemeris/tree/master/top2013 */
    /* Copied from: ftp://ftp.imcce.fr/pub/ephem/planets/top2013 */
    static const double freq[] =
    {
        0.5296909622785881e+03,
        0.2132990811942489e+03,
        0.7478166163181234e+02,
        0.3813297236217556e+02,
        0.2533566020437000e+02
    };
    int error = 1;
    int i, f, s, t;
    double time[TOP_MAX_SERIES];
    double el[6];
    double arg, dmu, xl;

    /* Check planet index */
    if (model->planet < 5 || model->planet > 9)
        FAIL("TopCalcElliptical: invalid planet = %d\n", model->planet);

    /* Time */
    time[0] = 1.0;
    time[1] = tt / 365250.0;
    for (i=1; i < TOP_MAX_SERIES; ++i)
        time[i] = time[i-1] * time[1];

    dmu = (freq[0] - freq[1]) / 880.0;

    for (f=0; f < TOP_NCOORDS; ++f)
    {
        const top_formula_t *formula = &model->formula[f];
        el[f] = 0.0;
        for (s=0; s < formula->nseries_calc; ++s)
        {
            const top_series_t *series = &formula->series[s];
            for (t=0; t < series->nterms_calc; ++t)
            {
                const top_term_t *term = &series->terms[t];
                if (f==1 && s==1 && term->k==0)
                    continue;
                arg = term->k * dmu * time[1];
                el[f] += time[s] * (term->c*cos(arg) + term->s*sin(arg));
            }
        }
    }

    xl = el[1] + freq[model->planet - 5] * time[1];
    xl = fmod(xl, dpi);
    if (xl < 0.0)
        xl += dpi;
    el[1] = xl;

    /* Convert elliptical elements from array 'el' to friendly struct layout. */
    ellip->a      = el[0];
    ellip->lambda = el[1];
    ellip->k      = el[2];
    ellip->h      = el[3];
    ellip->q      = el[4];
    ellip->p      = el[5];

    error = 0;
fail:
    return error;
}


int TopEcliptic(int planet, const top_elliptical_t *ellip, top_rectangular_t *ecl)
{
    static const double gmp[] =
    {
        0, /* placeholder */
        4.9125474514508118699e-11,
        7.2434524861627027000e-10,
        8.9970116036316091182e-10,
        9.5495351057792580598e-11,
        2.8253458420837780000e-07,
        8.4597151856806587398e-08,
        1.2920249167819693900e-08,
        1.5243589007842762800e-08,
        2.1886997654259696800e-12
    };
    static const double gmsol = 2.9591220836841438269e-04;
    double rgm;
    double xa, xl, xk, xh, xq, xp;
    double xfi, xki, u, ex, ex2, ex3;
    double zr, zi;
    double z1r, z1i;
    double z2r, z2i;
    double z3r, z3i;
    double zteta_r, zteta_i;
    double zto_r, zto_i;
    double gl, gm, e, dl, rsa;
    double xcw, xsw, xm, xr, xms, xmc, xn;
    int error = 1;

    if (planet < 1 || planet > 9)
        FAIL("TopEcliptic: invalid planet = %d\n", planet);

    rgm = sqrt(gmp[planet] + gmsol);
    xa = ellip->a;
    xl = ellip->lambda;
    xk = ellip->k;
    xh = ellip->h;
    xq = ellip->q;
    xp = ellip->p;

    xfi = sqrt(1.0 - xk*xk - xh*xh);
    xki = sqrt(1.0 - xq*xq - xp*xp);
    zr = xk; zi = xh;       /* z = dcmplx(xk,xh) */
    u = 1.0 / (1.0 + xfi);
    ex2 = zr*zr + zi*zi;
    ex = sqrt(ex2);         /* ex = cdabs(z) */
    ex3 = ex * ex2;
    z1r = zr; z1i = -zi;

    gl = fmod(xl, dpi);
    gm = gl - atan2(xh, xk);
    e = gl + (ex - 0.125*ex3)*sin(gm) + 0.5*ex2*sin(2.0*gm) + 0.375*ex3*sin(3.0*gm);

    do
    {
        z2r = 0.0; z2i = e;
        zteta_r = cos(z2i);
        zteta_i = sin(z2i);
        z3r = z1r*zteta_r - z1i*zteta_i;
        z3i = z1r*zteta_i + z1i*zteta_r;
        dl = gl - e + z3i;
        rsa = 1.0 - z3r;
        e += dl/rsa;
    } while (fabs(dl) >= 1.0e-15);

    z1r = z3i * u * zr;
    z1i = z3i * u * zi;
    z2r = +z1i;
    z2i = -z1r;
    zto_r = (-zr + zteta_r + z2r) / rsa;
    zto_i = (-zi + zteta_i + z2i) / rsa;
    xcw = zto_r;
    xsw = zto_i;
    xm = xp*xcw - xq*xsw;
    xr = xa*rsa;

    ecl->x = xr*(xcw - 2.0*xp*xm);
    ecl->y = xr*(xsw + 2.0*xq*xm);
    ecl->z = -2.0*xr*xki*xm;

    xms = xa*(xh + xsw)/xfi;
    xmc = xa*(xk + xcw)/xfi;
    xn = rgm / (xa * sqrt(xa));

    ecl->vx = xn*((2.0*xp*xp - 1.0)*xms + 2.0*xp*xq*xmc);
    ecl->vy = xn*((1.0 - 2.0*xq*xq)*xmc - 2.0*xp*xq*xms);
    ecl->vz = 2.0*xn*xki*(xp*xms + xq*xmc);

    error = 0;
fail:
    return error;
}


int TopEquatorial(const top_rectangular_t *ecl, top_rectangular_t *equ)
{
    double rot[3][3];
    const double pi = dpi / 2.0;
    const double dgrad = pi / 180.0;
    const double sdrad = dgrad / 3600.0;
    const double eps = (23.0 + 26.0/60.0 + 21.41136/3600.0)*dgrad;
    const double phi = -0.05188 * sdrad;
    const double ceps = cos(eps);
    const double seps = sin(eps);
    const double cphi = cos(phi);
    const double sphi = sin(phi);

    rot[0][0] =  cphi;
    rot[0][1] = -sphi*ceps;
    rot[0][2] =  sphi*seps;
    rot[1][0] =  sphi;
    rot[1][1] =  cphi*ceps;
    rot[1][2] = -cphi*seps;
    rot[2][0] =  0.0;
    rot[2][1] =  seps;
    rot[2][2] =  ceps;

    equ->x = (rot[0][0] * ecl->x) + (rot[0][1] * ecl->y) + (rot[0][2] * ecl->z);
    equ->y = (rot[1][0] * ecl->x) + (rot[1][1] * ecl->y) + (rot[1][2] * ecl->z);
    equ->z = (rot[2][0] * ecl->x) + (rot[2][1] * ecl->y) + (rot[2][2] * ecl->z);

    equ->vx = (rot[0][0] * ecl->vx) + (rot[0][1] * ecl->vy) + (rot[0][2] * ecl->vz);
    equ->vy = (rot[1][0] * ecl->vx) + (rot[1][1] * ecl->vy) + (rot[1][2] * ecl->vz);
    equ->vz = (rot[2][0] * ecl->vx) + (rot[2][1] * ecl->vy) + (rot[2][2] * ecl->vz);

    return 0;
}
