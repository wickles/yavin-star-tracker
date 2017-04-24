/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
*/

#pragma once

#ifndef IDENTIFIER_H
#define IDENTIFIER_H

#include "star.h"
#include <vector>
#include "catalog.h"

void IdentifyImageStars(std::vector<image_star>& ImageStars, Catalog& theCatalog, prior_s* Prior);

#endif IDENTIFIER_H