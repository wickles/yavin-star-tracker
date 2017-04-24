#pragma once

#ifndef IDENTIFIER_H
#define IDENTIFIER_H

#include "star.h"
#include <vector>
#include "catalog.h"

void IdentifyImageStars(std::vector<image_star>& ImageStars, Catalog& theCatalog, prior_s* Prior);

#endif IDENTIFIER_H