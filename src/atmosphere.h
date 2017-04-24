/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
*/

#include "star.h"

#pragma once

//more information on http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-sider-time/?searchterm=GAST
//compute the Hour Angle in Radians
double getLST(double LONG, double JDN);
double getHA1(double LST, double RA);
double getEle(double HA, double LAT, double DEC);
double getAzi(double HA, double LAT, double DEC);

double correctEle(double ELE, double height, double temp);

double getHA2(double AZI, double LAT, double ELE);
double getRA(double LST, double HA);
double getDEC(double AZI, double LAT, double ELE);