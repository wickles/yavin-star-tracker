#include "star.h"
#include <vector>

#pragma once

extern sfloat eigen_vector[4];

// Returns RMS error
double GetAttitude(std::vector<image_star>& ImageStars, std::vector<catalog_star>& CatalogStars, coordinates* CoordOut, sfloat RMat[3][3], double RQuat[4]);
