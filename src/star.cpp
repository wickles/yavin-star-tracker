/* star.cpp
Alexander Wickes and Gil Tabak
September 20, 2010
*/

#include "star.h"
#include <cmath>
#include <cstdio>
#include <cstdarg>

//distortion function.  This one doesn't make a very big difference.
static inline sfloat distortion_function(sfloat phi)
{
	return ( 0.9987f*phi - 0.0084f*phi*phi + 0.4382f*phi*phi*phi - 3.5699f*phi*phi*phi*phi + 11.634f*phi*phi*phi*phi*phi );
}

#ifdef DEBUG_TEXT
int debug_printf(char* format, ...) {
#ifndef DEBUG_TEXT
	/* Empty body, so a good compiler will optimise calls
	   to debug_printf away */
#else
        va_list args;
        va_start(args, format);
        int ret = vfprintf(stderr, format, args);
        va_end(args);
		return ret;
#endif /* DEBUG_TEXT */
}
#endif

double dot_product(double a[3], double b[3])						//returns a . b
{
	return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

void cross_product(double a[3], double b[3], double out[3])			//out = a x b
{
	double tmp[3];
	tmp[0] = a[1]*b[2] - a[2]*b[1];
	tmp[1] = a[2]*b[0] - a[0]*b[2];
	tmp[2] = a[0]*b[1] - a[1]*b[0];

	int i;
	for ( i = 0; i < 3; i++ )
		out[i] = tmp[i];
}

void matrix_mult(double M[3][3], double x[3], double out[3])		//out = M(x)
{
	int i;
	for ( i = 0; i < 3; i++ )
		out[i] = dot_product( M[i], x );
}

double magnitude(double x[3])										//returns the magnitude of x.
{
	return sqrt( dot_product(x,x) );
}

double normalize(double x[3], double out[3])						//out = x/|x| if x != 0. 
{
	double mag = magnitude(x);
	int i;
	if ( abs(mag) < 0.00000001 )									//if the vector is very small
	{
		for ( i = 0; i < 3; i++ )									//check to see it it is == 0
			out[i] = 0.0;
		return 0.0;													//if so, return 0.
	}

	for ( i = 0; i < 3; i++ )										//otherwise, normalize x.
		out[i] = x[i] / mag;

	return mag;
}

void rotate(double x[3], double r[3], double theta, double out[3])	//out = rotation of x around r by theta, assuming |r| = 1
{
	double cross[3];				
	cross_product( r, x, cross );

	double dot = dot_product(r,x);
	double c = cos(theta);
	double s = sin(theta);

	int i;
	for ( i = 0; i < 3; i++ )
		out[i] = x[i]*c + cross[i]*s + r[i]*dot*(1.0 - c);			//use the Rodrigues rotational forumla.
}

//input x and y, cartesian coordinates from the center of an image, as well as the focal length, and make rOut = (x,y,z) on a unit sphere.
void GetSphericalFromImage(double x, double y, float FocalLength, double rOut[3])	
{
	if ( rOut == NULL )
		return;

	double r = sqrt( x*x + y*y );									//distance from the center of the image

	double phi = atan( r / FocalLength );							//angular distance away from the center of the image (point of tangency with the sphere)
	phi = distortion_function( phi );								//take into account radial distortion.
	double theta = atan2( y, x ) + (double)M_PI;					//angle from the x-axis of the planar image to where the star appears on the image.
	double delta = asin( sin(theta) * sin(phi) );					//law of sines yields altitude
	double lambda = acos( cos(phi) / cos(delta) );					//law of cosines yields longitude.
	if ( x < 0 )													//sign ambiguity due to the cosine is resolved.
		lambda = -lambda;

	rOut[0] = cos(delta) * cos(lambda);								//output as cartesian coordinates.
	rOut[1] = cos(delta) * sin(lambda);
	rOut[2] = sin(delta);
}

void GetCoordsFromSpherical(double r[3], coordinates* Out)			//gets RA and DEC in Out from the unit vector r[3].
{
	if ( Out == NULL || r == NULL )
		return;

	Out->RA = atan2( r[1], r[0] );									//get RA from r[0] and r[1]
	if ( Out->RA < 0 )												//make sure RA is in [0,2pi).
		Out->RA += 2 * M_PI;
	Out->DEC = atan( r[2] / sqrt( SQUARE(r[0]) + SQUARE(r[1]) ) );	//find the declination from r.
}


void GetDiscreteCoords(coordinates* coords, coords_discrete* out)
{
	if (coords == NULL || out == NULL)
		return;

	out->RA_hr = int( coords->RA/(2*M_PI)*24 );
	out->RA_min = int( coords->RA/(2*M_PI)*24*60 ) % 60;
	out->RA_sec = fmod( coords->RA/(2*M_PI)*24*60*60, 60 );
	out->DEC_deg = int( coords->DEC*(180/M_PI) );
	out->DEC_min = int( abs(coords->DEC)*(180/M_PI)*60 ) % 60;
	out->DEC_sec = fmod( abs(coords->DEC)*(180/M_PI)*60*60, 60);
}