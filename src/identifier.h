#include "star.h"
#include <vector>
#include "catalog.h"

#pragma once

void IdentifyImageStars(std::vector<image_star>& ImageStars, Catalog& theCatalog, prior_s* Prior);