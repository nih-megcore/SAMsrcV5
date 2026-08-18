#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <samlib.h>
#include <sam_parse.h>
#include <siglib.h>

#define EPS 1.0e-10

static int failures;

#define CHECK(condition, message) do { \
    if (!(condition)) { \
        fprintf(stderr, "FAIL: %s (line %d)\n", message, __LINE__); \
        failures++; \
    } \
} while (0)

static int close_enough(double actual, double expected, double tolerance)
{
    return fabs(actual - expected) <= tolerance;
}

static double dot3(const double *a, const double *b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void test_coordinates(void)
{
    double cart[3] = { 1.25, -2.5, 4.0 };
    double spherical[3];
    double roundtrip[3];
    double radial[3] = { 1.0, 2.0, 3.0 };
    double xaxis[3];
    double yaxis[3];
    int i;

    CtoS(cart, spherical);
    StoC(spherical, roundtrip);
    for (i = 0; i < 3; i++)
        CHECK(close_enough(roundtrip[i], cart[i], EPS), "Cartesian/spherical round trip");

    OrthoPlane(radial, xaxis, yaxis);
    CHECK(close_enough(dot3(radial, radial), 1.0, EPS), "radial vector normalized");
    CHECK(close_enough(dot3(xaxis, xaxis), 1.0, EPS), "first tangent normalized");
    CHECK(close_enough(dot3(yaxis, yaxis), 1.0, EPS), "second tangent normalized");
    CHECK(close_enough(dot3(radial, xaxis), 0.0, EPS), "radial perpendicular to first tangent");
    CHECK(close_enough(dot3(radial, yaxis), 0.0, EPS), "radial perpendicular to second tangent");
    CHECK(close_enough(dot3(xaxis, yaxis), 0.0, EPS), "tangent vectors perpendicular");
}

static void test_rotations_and_fields(void)
{
    double euler[3] = { 0.2, -0.3, 0.4 };
    double rotation[3][3];
    double row_dot;
    double values[6] = { 1., 2., 3., -4., 5., -6. };
    double copied[6] = { 0. };
    FIELD field;
    gsl_vector *vector = gsl_vector_alloc(6);
    int i, j, k;

    EtoR(euler, rotation);
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            row_dot = 0.;
            for (k = 0; k < 3; k++)
                row_dot += rotation[i][k] * rotation[j][k];
            CHECK(close_enough(row_dot, i == j ? 1.0 : 0.0, EPS), "rotation matrix orthonormal");
        }
    }

    A6toFIELD(values, &field);
    FIELDtoA6(&field, copied);
    for (i = 0; i < 6; i++)
        CHECK(close_enough(copied[i], values[i], EPS), "FIELD/array round trip");

    FIELD2Vector(&field, vector);
    field.p[0] = field.p[1] = field.p[2] = 0.;
    field.v[0] = field.v[1] = field.v[2] = 0.;
    Vector2FIELD(vector, &field);
    FIELDtoA6(&field, copied);
    for (i = 0; i < 6; i++)
        CHECK(close_enough(copied[i], values[i], EPS), "FIELD/GSL vector round trip");
    gsl_vector_free(vector);
}

static void test_signals(void)
{
    double centered[4] = { 1., 2., 3., 4. };
    double trend[5] = { 3., 5., 7., 9., 11. };
    double sequence[12];
    double window[5];
    double increasing[5] = { 1., 2., 3., 4., 5. };
    double decreasing[5] = { 5., 4., 3., 2., 1. };
    double sum;
    int i;

    demean(centered, 4);
    for (i = 0, sum = 0.; i < 4; i++)
        sum += centered[i];
    CHECK(close_enough(sum, 0.0, EPS), "demean removes mean");

    detrend(trend, 5);
    for (i = 0; i < 5; i++)
        CHECK(close_enough(trend[i], 0.0, EPS), "detrend removes linear signal");

    for (i = 0; i < 12; i++)
        sequence[i] = (double)i;
    CHECK(close_enough(power(sequence, 0, 11), 143.0 / 12.0, EPS), "power computes population variance");

    Hanning(window, 5);
    CHECK(close_enough(window[0], 0.0, EPS), "Hanning starts at zero");
    CHECK(close_enough(window[1], 0.5, EPS), "Hanning quarter point");
    CHECK(close_enough(window[2], 1.0, EPS), "Hanning midpoint");
    CHECK(close_enough(window[3], 0.5, EPS), "Hanning symmetry");
    CHECK(close_enough(window[4], 0.0, EPS), "Hanning ends at zero");

    CHECK(close_enough(kendall(increasing, increasing, 5), 1.0, EPS), "Kendall positive rank correlation");
    CHECK(close_enough(kendall(increasing, decreasing, 5), -1.0, EPS), "Kendall negative rank correlation");
}

static void test_linear_algebra(void)
{
    gsl_matrix *matrix = gsl_matrix_calloc(2, 2);
    gsl_matrix *inverse = gsl_matrix_calloc(2, 2);
    gsl_vector *lead = gsl_vector_alloc(2);
    gsl_vector *weights = gsl_vector_alloc(2);
    double source_power;
    double noise_power;

    gsl_matrix_set(matrix, 0, 0, 2.0);
    gsl_matrix_set(matrix, 1, 1, 4.0);
    pinv(matrix, inverse);
    CHECK(close_enough(gsl_matrix_get(inverse, 0, 0), 0.5, EPS), "pseudoinverse first diagonal");
    CHECK(close_enough(gsl_matrix_get(inverse, 1, 1), 0.25, EPS), "pseudoinverse second diagonal");
    CHECK(close_enough(gsl_matrix_get(inverse, 0, 1), 0.0, EPS), "pseudoinverse off diagonal");

    gsl_vector_set_all(lead, 1.0);
    SAMsolve(matrix, inverse, lead, weights, 9.0, &source_power, &noise_power);
    CHECK(close_enough(gsl_vector_get(weights, 0), 2.0 / 3.0, EPS), "SAM first weight");
    CHECK(close_enough(gsl_vector_get(weights, 1), 1.0 / 3.0, EPS), "SAM second weight");
    CHECK(close_enough(source_power, 4.0 / 3.0, EPS), "SAM source power");
    CHECK(close_enough(noise_power, 5.0, EPS), "SAM noise power");

    gsl_vector_free(weights);
    gsl_vector_free(lead);
    gsl_matrix_free(inverse);
    gsl_matrix_free(matrix);
}

static void test_sam_directories(void)
{
    PARMINFO params;
    char path[512];
#ifndef _WIN32
    char temp[] = "/tmp/sam-path-test-XXXXXX";
    char nested[512];
#endif
    char *input;
    char *output;
    char *argv[] = {
        "legacy", "-r", "dataset", "-i_SAMdir", "input-root",
        "-o_SAMdir", "output-root", "-v", NULL
    };
    int argc = 8;

    new_params(&params);
    params.DataSetName = "/data/example.ds";
    GetSAMPath(path, sizeof(path), &params, SAM_INPUT);
    CHECK(strcmp(path, "/data/example.ds/SAM") == 0, "default SAM input root");
    GetSAMPath(path, sizeof(path), &params, SAM_OUTPUT);
    CHECK(strcmp(path, "/data/example.ds/SAM") == 0, "default SAM output root");

    params.InputSAMDirectory = "separate-input";
    params.OutputSAMDirectory = "separate-output";
    GetSAMPath(path, sizeof(path), &params, SAM_INPUT);
    CHECK(strcmp(path, "separate-input") == 0, "configured SAM input root");
    GetSAMPath(path, sizeof(path), &params, SAM_OUTPUT);
    CHECK(strcmp(path, "separate-output") == 0, "configured SAM output root");

#ifndef _WIN32
    CHECK(mkdtemp(temp) != NULL, "create SAM directory test root");
    snprintf(nested, sizeof(nested), "%s/one/two/SAM", temp);
    CHECK(makedirs(nested) == 0, "recursively create SAM output root");
    CHECK(direxists(nested), "recursive SAM output root exists");

    input = output = NULL;
    parse_samdir_args(&argc, argv, &input, &output);
    CHECK(argc == 4, "legacy SAM flags removed before getopt");
    CHECK(strcmp(input, "input-root") == 0, "legacy SAM input flag parsed");
    CHECK(strcmp(output, "output-root") == 0, "legacy SAM output flag parsed");
    CHECK(strcmp(argv[1], "-r") == 0 && strcmp(argv[3], "-v") == 0,
          "unrelated legacy options preserved");

    rmdir(nested);
    snprintf(nested, sizeof(nested), "%s/one/two", temp);
    rmdir(nested);
    snprintf(nested, sizeof(nested), "%s/one", temp);
    rmdir(nested);
    rmdir(temp);
#endif
}

int main(void)
{
    test_coordinates();
    test_rotations_and_fields();
    test_signals();
    test_linear_algebra();
    test_sam_directories();

    if (failures != 0) {
        fprintf(stderr, "%d core unit test(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    printf("core C unit tests passed\n");
    return EXIT_SUCCESS;
}
