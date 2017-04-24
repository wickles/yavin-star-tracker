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

extern sfloat eigen_vector[4];

// Returns RMS error
double GetAttitude(	std::vector<image_star>& ImageStars, std::vector<catalog_star>& CatalogStars, 
					coordinates* BoresightOut, coordinates* CoordOut, sfloat RMat[3][3], double RQuat[4], double &err);