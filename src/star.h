#pragma once

//#define DEBUG_TEXT
//#define PROMPT_USER
#define SUBTRACT_DARK
//#define OUTPUT_IMAGES
#define PRINT_TOP_VOTES
//#define ONLY_ONE_STAGE
//#define WEIGH_STAR_ERROR
#define USE_PRIOR

#define M_PI			3.14159265358979323846264338327950288

#define PRIOR_TIMEOUT		10 /* seconds */

#define LOW_VOTES_COEFF		0.6f
#define MAX_VOTES_COEFF		0.75f

#define GL_SAMPLE_SKIP		16
#define LOCAL_WIDTH			64
#define LOCAL_HEIGHT		64
#define LOCAL_SAMPLE_SKIP	4
#define STAR_MIN_OUTSTND	4
#define STAR_MAX_OUTSTND	100

#define MAX_IMAGE_STARS		20
#define DOUBLE_STAR_PIXEL_RADIUS	10
#define IMAGE_STARS_USE_MAX_CUTOFF	15

#define MAX_MAGNITUDE	10.0f
#define MAX_ANG_DIST	0.30f
#define MAX_RADIUS		5.0f
#define MAX_PRIOR_DIST	(1.5f * MAX_ANG_DIST)

#define ERROR_FUNCTION(phi) ( .000002 / sqrt(phi) )

#define IMAGE_WIDTH		1280
#define IMAGE_HEIGHT	960
#define IMAGE_PIXELS	(IMAGE_WIDTH*IMAGE_HEIGHT)

#define INDEX_INVALID		(-1)
#define INDEX_FALSE_STAR	(-2)

#define STAR_USE_DOUBLE

typedef
#ifdef STAR_USE_DOUBLE
	double
#else
	float
#endif
	sfloat;

typedef short index_t;

struct catalog_star {
	sfloat RA, Dec;
	float app_mag;
	sfloat r[3];
	char name[11];
};

struct catalog_pair {
	index_t star1, star2;
	sfloat distance;
};

struct image_star {
	sfloat centroid_x, centroid_y;
	sfloat r_prime[3];
	sfloat r[3];
	sfloat error;
	sfloat radius;
	index_t identity;
};

struct image_pair {
	index_t star1, star2;
	sfloat lower, upper;
};

struct coordinates {
	sfloat RA, DEC;
};

struct coords_discrete {
	char RA_hr, RA_min;
	float RA_sec;
	char DEC_deg, DEC_min;
	float DEC_sec;
};

struct prior_s {
	coordinates coords;
	sfloat phi_max;
	unsigned int tickCount;
};

#ifndef DEBUG_TEXT
#define debug_printf( ... ) ((void)0)
#else
int debug_printf(char *format, ...);
/* print a message, if it is considered significant enough.
      Adapted from [K&R2], p. 174 */
#endif

// Define math functions and macros

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define SQUARE(x) ((x)*(x))

double dot_product(double a[3], double b[3]);
void cross_product(double a[3], double b[3], double out[3]);
void matrix_mult(double M[3][3], double x[3], double out[3]);
double magnitude(double x[3]);
double normalize(double x[3], double out[3]);

// Rotates x about r by theta (right handed). Assumes r is a unit vector.
void rotate(double x[3], double r[3], double theta, double out[3]);

void GetSphericalFromImage(double x, double y, float FocalLength, double rOut[3]);
void GetCoordsFromSpherical(double r[3], coordinates* Out);

void GetDiscreteCoords(coordinates* coords, coords_discrete* out);
