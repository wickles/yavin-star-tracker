#include "star.h"

#pragma once

//more information on http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-sider-time/?searchterm=GAST
//compute the Hour Angle in Radians
double getLST(double RA, double LONG, double relJDN);
double getHA1(double LST, double RA);
double getEle(double HA, double LAT, double DEC);
double getAzi(double HA, double LAT, double DEC);

double correctEle(double ELE, double height, double temp);

double getHA2(double AZI, double LAT, double ELE);
double getRA(double LST, double HA);
double getDEC(double AZI, double LAT, double ELE);