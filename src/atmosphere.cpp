/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
*/

/* Atmospheric_Correction.cpp
Gil Tabak 2/17/2011
This file uses the system time (to be replaced with GPS time) and the lattitude and longitude 
read by the main program from the Settings.txt file (also to be replaced by GPS readings) to calculate
the azimuth and elevation angles of each of the stars.  An atmospheric correction is applied to find
the new azimuth and elevation angles, which are then transformed back to (corrected) RA and DEC.
 */

//**IMPORTANT NOTE**/
//RA and DEC are in RADIANs; LONG and LAT are in DEGREES (RA and DEC are computed; lat and long are inputs)
//LONG is taken to be positive for east, negative for west

// TEST: RA anything, DEC +- 90 deg. --> AZI 0, ELE = LAT

#include "atmosphere.h"
#include <cmath>

using namespace std;

//more information on http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-sider-time/?searchterm=GAST
double getLST(double LONG, double JDN) //compute the Hour Angle in Radians
{
	double D, GMST, LST;
	D = JDN - 2451545.0;													//Julian Date since J2000 (in days)
	GMST = fmod(18.697374558 + 24.06570982441908 * D, 24);					//compute the Greenwtch Mean Sidereal Time (in hours)
	LST = LONG * 12 / M_PI + GMST;											//compute the mean Local Sidereal Time (in hours)
	if (LST < 0)															//make sure LST in the right range (0-24 hours)
		LST += 24;
	if (LST > 24)
		LST -= 24;
	return LST * M_PI / 12;													//return LST in ***RADIANS***
}

double getHA1(double LST, double RA)
{
	return LST - RA;														//compute the hour angle **In RADIANS**
}

//the transformation can be found at http://www.astrobio.nau.edu/~koerner/ast301/lecture/lecture1/lec1.html
//or derived using spherical trigonometry.
double getEle(double HA, double LAT, double DEC)
{
	return asin(cos(HA) * cos(DEC) * cos(LAT) + sin(DEC) * sin(LAT) );
}

double getAzi(double HA, double LAT, double DEC)
{
	double AZI = atan2( -sin(HA), tan(DEC) * cos(LAT) - sin(LAT) * cos(HA) );
	if (AZI < 0)
		AZI = AZI + 2*M_PI;
	return AZI;
}

//based on the information found at http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
//this model uses two different approximations for different cases, when the elevation
//is greater than 5 degrees or less than 5 degrees.
//This model roughly corrects for pressure and temperature.  A better model would
//correct for humidity and wavelength as well.

//note: height is in METERS, and temp is in DEGREES CELSIUS
double correctEle(double ELE, double height, double temp)
{
	double correction;
	double th = tan(ELE);
	double h = height;
	if (ELE = M_PI/2)
		return ELE;
	else if (ELE > 0.0873)
		correction = ( (58.1/th) - (0.07/(th*th*th)) + (0.000086/(th*th*th*th*th)) ) / 3600.0;
	else 
		correction = (1735.0 - 518.2*h + 103.4*h*h - 12.79*h*h*h + 0.711*h*h*h*h) / 3600.0;

	correction = correction*pow(1-0.0000225577*height,5.25588)*283/(273+temp); // correct for temperature and pressure

	return ELE + correction * (M_PI / 180); 
}

//The backwards transformations can be written essentially the same way.
//You can use this to make the code more elegant, but much more confusing.

double getHA2(double AZI, double LAT, double ELE)
{
	double HA = atan2( -sin(AZI), tan(ELE) * cos(LAT) - sin(LAT) * cos(AZI) );
	if (HA < 0)
		HA += 2*M_PI;
	return HA;
}

double getDEC(double AZI, double LAT, double ELE)
{
	return asin( cos(AZI) * cos(ELE) * cos(LAT) + sin(ELE) * sin(LAT) );
}

double getRA(double LST, double HA)
{
	return LST - HA;
}
