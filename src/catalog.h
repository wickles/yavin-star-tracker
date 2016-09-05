/* catalog.cpp
Alexander Wickes and Gil Tabak
September 20, 2010
*/

#include "star.h"
#include <vector>
#include <Windows.h>

#pragma once

/*
This structure is meant to contain a star catalog and an associated set of star pairs,
which are loaded when Initialize is called. The catalog is loaded from the entries in
catalog.txt, excluding stars with magnitude higher than the value specified in star.h.
Pairs are created for all unique pair of catalog stars 
*/
struct Catalog {
	bool initialized;
	std::vector<catalog_star> Stars;
	std::vector<catalog_pair> SortedPairs;

	Catalog();
	void Initialize( const char* filename, SYSTEMTIME* systime, float FocalLength );
};
