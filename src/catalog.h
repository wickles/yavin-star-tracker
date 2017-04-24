/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
*/

#include "star.h"
#include <vector>

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
