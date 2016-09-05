#include "star.h"
#include <vector>

#pragma once

extern sfloat eigen_vector[4];

// Returns RMS error
double GetAttitude(std::vector<image_star>& ImageStars, std::vector<catalog_star>& CatalogStars, coordinates* CoordOut, sfloat RotationOut[3][3]);
