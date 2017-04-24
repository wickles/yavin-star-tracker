#include "star.h"
#include <vector>

#pragma once

struct image_s {
	int width;
	int height;
	int bitmode;
	float FocalLength;
	void* data;
	void* out;
};

struct detector_s {
	int gl_sample_skip;		// distance to next pixel in global sample
	int local_width;		// box width for local mean and standard deviation
	int local_height;		// box height for above
	int local_sample_skip;	// distance to next pixel in local sample
	int star_min_outstnd;
	sfloat mean_sky;
};

// ImageStars -- vector to be cleared and filled with image stars
// Detector -- pointer to struct containing settings for detector
// Image -- pointer to struct containing image data & attributes
// returns number of stars detected in image, not number of stars in data structure (usually gets truncated)
size_t DetectStars(std::vector<image_star>& ImageStars, detector_s* Detector, image_s* Image, double out_radius, double in_radius);
