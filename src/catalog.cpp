#include "catalog.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include "timer.h"

using namespace std;

// determines the julian date from computer clock
static double getJulianDate(SYSTEMTIME* systime)
{
	int A,B,C,E,F;																	//used for computing the Julian Date
	double JDN;

	A = (int) (systime->wYear/100);
	B = (int) (A/4);
	C = 2 - A + B;
	E = (int) (365.25 * (systime->wYear + 4716));
	F = (int) (30.6001 * (systime->wMonth + 1));
	JDN = (double)(C + E + F + systime->wDay - 1524.5);

	return JDN	+ (double)systime->wHour/24 + (double)systime->wMinute/(24*60)
				+ (double)systime->wSecond/(24*60*60) + (double)systime->wMilliseconds/(24*60*60*1000);
}

// arguments in RADIANS (or RADIANS / YEARS). Takes RA,DEC,proper motion, and uses time to find new RA and DEC. Puts that in tempRA and tempDEC
static void updateCoordinates(double* RA, double* Dec, double pmRA, double pmDEC, double t, double zeta, double z, double theta)
{
	//printf("RA = %f Dec = %f\n", RA,Dec);

	*RA += pmRA * t * 100; 															//update RA and DEC due to proper motion.
	*Dec += pmDEC * t * 100;

	double A = cos(*Dec) * sin(*RA + zeta);
	double B = cos(theta) * cos(*Dec) * cos(*RA + zeta) - sin(theta) * sin(*Dec);
	double C = sin(theta) * cos(*Dec) * cos(*RA + zeta) + cos(theta) * sin(*Dec);	//A, B, and C are used for precession.

	//printf("A = %f B = %f C = %f\n",A,B,C);

	*RA = z + atan2(A,B);															//update RA and DEC based on precession.
	if (*RA < 0)
		*RA += 2 * M_PI;
	*Dec = asin(C);
}

void ReadCatalog( vector<catalog_star>& Catalog, const char* filename, SYSTEMTIME* systime )
{
	Timer timer;
	timer.StartTimer();

	//vector<catalog_star> Catalog;
	Catalog.clear();
	Catalog.reserve(10000);

	//the following variables are used to compute precession:
	double julianDate = getJulianDate(systime);
	debug_printf("Julian Date: %f\n", julianDate);
	double t = (julianDate - 2451545.0) / 36525.0; //time since J2000
	double zeta = ( 2306.2181*t + 0.30188*t*t + 0.017998*t*t*t ) / (60*60) * (M_PI/180);
	double z = ( 2306.2181*t + 1.09468*t*t + 0.018203*t*t*t ) / (60*60) * (M_PI/180);
	double theta = ( 2004.3109*t - 0.42665*t*t - 0.041833*t*t*t ) / (60*60) * (M_PI/180);

	//printf("t = %10f zeta = %f z = %f theta = %f\n",t,zeta,z,theta);

	FILE* file = fopen(filename, "r");
	if (file == NULL)
	{
		printf("Error: Could not open catalog.\n");
		//return Catalog;
		return;
	}
	
	char buffer[256];
	while ( fgets(buffer, 256, file) != NULL )
	{
		if ( strlen(buffer) >= 107 )
		{
			// read RA, Dec, apparent magnitude, proper motion
			int RA_hour, RA_min, DEC_deg, DEC_min, DEC_sec;
			float RA_sec, app_mag, pmRA, pmDec;
			char DEC_sign;
			char name[11];
			int ret = 0;
			ret += sscanf(	buffer+4,	"%10c", name );
			name[10] = '\0';
			ret += sscanf(	buffer+75,	"%2d%2d%4f%c%2d%2d%2d", &RA_hour, &RA_min, &RA_sec, &DEC_sign, &DEC_deg, &DEC_min, &DEC_sec );
			ret += sscanf(	buffer+102,	"%4f", &app_mag );
			ret += sscanf(	buffer+148,	"%6f%6f", &pmRA, &pmDec );
/*
			int ret = sscanf(	buffer, "%*75c%2d%2d%4f%c%2d%2d%2d%*13c%4f%*41c%6f%6f",
								&RA_hour, &RA_min, &RA_sec, &DEC_sign, &DEC_deg, &DEC_min, &DEC_sec, &app_mag, &pmRA, &pmDec );
*/

			//printf("(Original) %02d:%02d:%02.1f, %c%02d:%02d:%02d, %.2f, %+.3f, %+.3f\n", RA_hour, RA_min, RA_sec, DEC_sign, DEC_deg, DEC_min, DEC_sec, app_mag, pmRA, pmDec );

			// checks for incomplete data and ignores magnitudes higher than max
			if (ret == 11 && app_mag <= MAX_MAGNITUDE )
			{
				double RA = 2 * M_PI * ( (sfloat)RA_hour/24 + (sfloat)RA_min/(24*60) + (sfloat)RA_sec/(24*60*60) );
				double DEC = M_PI * ( (sfloat)DEC_deg/180 + (sfloat)DEC_min/(180*60) + (sfloat)DEC_sec/(180*60*60) );
				if (DEC_sign == '-')
					DEC = -DEC;

				// units are arcsec/yr, convert arcsec -> deg -> rad
				double pmRA_d = pmRA / (60*60) * (M_PI/180);
				double pmDec_d = pmDec / (60*60) * (M_PI/180);

				//printf( "(Converted) RA: %f, Dec: %f, AppMag: %.2f, pmRA: %.15f, pmDec: %.15f\n", star.RA, star.Dec, star.app_mag, pmRA, pmDec );
				updateCoordinates( &RA, &DEC, pmRA_d, pmDec_d, t, zeta, z, theta );
				//printf( "(Updated) RA: %f, Dec: %f, AppMag: %.2f\n", star.RA, star.Dec, star.app_mag );

				//give each star its properties
				catalog_star star;
				star.RA = (sfloat)RA;
				star.Dec = (sfloat)DEC;
				star.app_mag = app_mag;
				star.r[0] = cos(star.Dec)*cos(star.RA);
				star.r[1] = cos(star.Dec)*sin(star.RA);
				star.r[2] = sin(star.Dec);
				strcpy(star.name, name);

				//put the star in the catalog.
				Catalog.push_back(star);

			}
		}
	}

	fclose(file);

	timer.StopTimer();

#if 0
	vector<catalog_star>::iterator it;
	for ( it = Catalog.begin(); it != Catalog.end(); it++ )
		printf(	"%04d | (RA, Dec) = (%f, %+f) | AppMag: %.2f | r = (%+f,%+f,%+f)\n",
				int(it - Catalog.begin()), it->RA, it->Dec, it->app_mag, it->r[0], it->r[1], it->r[2] );
#endif

	debug_printf("Catalog successfully read, size: %d | Time elapsed: %d ms\n", Catalog.size(), timer.GetTime());
#ifdef PROMPT_USER
	printf("Press ENTER to continue.\n");
	getchar();
#endif

	//return Catalog;
}

static bool compare_catalog(catalog_pair first, catalog_pair second)
{
	return (first.distance < second.distance);
}

// Puts the sorted catalog star pairs into CatalogPairs, excluding double stars and pairs farther than a certain distance away
void GetSortedCatalogPairs( vector<catalog_pair>& CatalogPairs, vector<catalog_star>& Catalog, float FocalLength )
{
	Timer timer;
	timer.StartTimer();

	const sfloat dot_lower = cos( MAX_ANG_DIST );											//bound on large angles can't see both stars in the same image.
	const sfloat dot_upper = cos( atan((sfloat)DOUBLE_STAR_PIXEL_RADIUS/FocalLength) );		//bound on small angles to avoid double stars

	debug_printf("Creating catalog pairs in Range = [%f, %f]\n", dot_lower, dot_upper );

	//vector<catalog_pair> CatalogPairs;
	CatalogPairs.clear();
	CatalogPairs.reserve(1000000);

	index_t i, j;
	for ( i = 0; i < Catalog.size() - 1; i++ )												//first star
	{
		for ( j = i + 1; j < Catalog.size(); j++ )											//second star
		{
			sfloat dot = dot_product(Catalog[i].r, Catalog[j].r);
			if ( dot_lower <= dot && dot <= dot_upper)										//the conditions are met
			{
				catalog_pair pair = { i, j, dot };											//make a pair
				CatalogPairs.push_back(pair);												//add the pair to the vector
			}

			//printf("(%04d,%04d) Dot: %+.17f | AngDist: %+.17f\n", i, j, dot, acos(dot));
		}
	}

	sort( CatalogPairs.begin(), CatalogPairs.end(), compare_catalog );						//sort the vector by distances.

	timer.StopTimer();

#if 0
	vector<catalog_pair>::iterator it;
	for ( it = CatalogPairs.begin(); it != CatalogPairs.end(); it++ )
		printf("%05d (%04d, %04d) AngDist: %.17f\n", int(it - CatalogPairs.begin()), it->star1, it->star2, it->ang_dist);
#endif

	debug_printf("Successfully created and sorted CatalogPairs, size: %d | Time elapsed: %d ms\n", CatalogPairs.size(), timer.GetTime());
#ifdef PROMPT_USER
	printf("Press ENTER to continue.\n");
	getchar();
#endif

	//return CatalogPairs;
}

Catalog::Catalog()
{
	initialized = false;
}

void Catalog::Initialize( const char* filename, SYSTEMTIME* systime, float FocalLength )
{
	ReadCatalog( Stars, filename, systime );
	GetSortedCatalogPairs( SortedPairs, Stars, FocalLength );
	initialized = true;
}
