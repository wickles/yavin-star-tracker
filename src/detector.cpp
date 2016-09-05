/* detector.cpp
Alexander Wickes and Gil Tabak
September 20, 2010
*/

#include "detector.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "timer.h"

using namespace std;

static inline short get_pixel(image_s* Image, int x, int y)
{
	if (Image->bitmode == 8)
		return (short)((unsigned char*)Image->data)[y*Image->width+x];
	return ((short*)Image->data)[y*Image->width+x];
}

static inline void set_pixel(image_s* Image, int x, int y, short val)
{
	if ( Image->out == NULL )
		return;
	if (Image->bitmode == 8)
		((unsigned char*)Image->out)[y*Image->width+x] = (unsigned char)val;
	((short*)Image->out)[y*Image->width+x] = val;
}

static inline bool is_marked(image_s* Image, vector<bool>& marked, int x, int y)
{
	return marked[y*Image->width+x];
}

static inline void mark(image_s* Image, vector<bool>& marked, int x, int y)
{
	marked[y*Image->width+x] = true;
}

sfloat inline error_function(sfloat phi)
{
	return ERROR_FUNCTION(phi);
}

static inline sfloat mean_over_area(image_s* Image, int x, int y, int radius)
{
	int counter = 0;
	int sum = 0;
	int x1, y1;
	for ( x1 = x - radius; x1 <= x + radius; x1++ )
	{
		if ( x1 < 0 || x1 >= Image->width )
			continue;
		for ( y1 = y - radius; y1 <= y + radius; y1++ )
		{
			if ( y1 < 0 || y1 >= Image->height )
				continue;
			sum += get_pixel(Image, x1, y1);
			counter += 1;
		}
	}

	return ( counter > 0 ? (sfloat)sum/counter : 0 );
}

static inline bool isOutstanding( image_s* Image, int x, int y, int low_val )
{
	return ( get_pixel(Image, x, y) >= low_val );
}

//detector_recurse( Image, marked, low_val, sample_min, x, y, &total_count, &outstnd_count, &x_numer, &y_numer, &sum );

static void detector_recurse(	image_s* Image, vector<bool>& marked, int low_val, short sample_min,
								int x, int y, int* total_count, int* outstnd_count,
								long long* x_numer, long long* y_numer, long long* sum )
{
	if ( x <= 0 || x >= Image->width-1 || y <= 0 || y >= Image->height-1 )
		*outstnd_count = -1;

	if ( *outstnd_count < 0 || is_marked(Image, marked, x, y) || *outstnd_count >= STAR_MAX_OUTSTND )
		return;

	short pixel = get_pixel(Image, x, y);
	bool outstanding = isOutstanding( Image, x, y, low_val );

	if ( outstanding || *outstnd_count >= 1 )
	{
		mark(Image, marked, x, y);
		int val = (int)pixel - sample_min;

		*x_numer += val * x;
		*y_numer += val * y;
		*sum += pixel;
		*total_count += 1;

		if ( outstanding )
		{
			*outstnd_count += 1;
			set_pixel(Image, x, y, 0xfff);

			int i,j;
			for ( i = -1; i <= 1; i++ )
			{
				for ( j = -1; j <= 1; j++ )
				{
					if (i == 0 && j == 0)
						continue;
					detector_recurse(	Image, marked, low_val, sample_min, x+i, y+j,
										total_count, outstnd_count, x_numer, y_numer, sum );
				}
			}
		}
	}
}

static inline void CentroidAt(	image_s* Image, vector<bool>& marked, vector<image_star>& ImageStars,
								int low_val, short sample_min, int x, int y )
{
	int total_count = 0;
	int outstnd_count = 0;
	long long sum = 0;
	long long x_numer = 0;
	long long y_numer = 0;

	detector_recurse( Image, marked, low_val, sample_min, x, y, &total_count, &outstnd_count, &x_numer, &y_numer, &sum );

	unsigned long long denom = sum - ( total_count * sample_min );

	// Area = pi*r^2
	sfloat radius = sqrt( outstnd_count / (sfloat)M_PI );

	if ( denom > 0 && outstnd_count >= STAR_MIN_OUTSTND && radius <= MAX_RADIUS )
	{
		image_star star;
		star.centroid_x = x_numer / (sfloat)denom;
		star.centroid_y = IMAGE_HEIGHT-1 - ( y_numer / (sfloat)denom );
		star.identity = INDEX_INVALID;
		star.error = error_function( atan( radius / Image->FocalLength ) );
		star.radius = radius;

		// process image distortion, find angular distances

		GetSphericalFromImage( star.centroid_x - IMAGE_WIDTH/2.0f + 0.5f, star.centroid_y - IMAGE_HEIGHT/2.0f + 0.5f, Image->FocalLength, star.r_prime );

		ImageStars.push_back(star);
	}
}

bool compare_radius(image_star& left, image_star& right)
{
	return (left.radius > right.radius);
}

size_t DetectStars(vector<image_star>& ImageStars, detector_s* Detector, image_s* Image)
{
	Timer timer;
	timer.StartTimer();

	unsigned int num_image_pixels = Image->width*Image->height;
	unsigned int num_gl_sample_pixels = (Image->width/Detector->gl_sample_skip-1) * (Image->height/Detector->gl_sample_skip-1);

	vector<bool> marked(num_image_pixels, false);

	//vector<image_star> ImageStars;
	ImageStars.clear();
	//ImageStars.reserve(500);

	vector<short> LocalSample;
	LocalSample.reserve( (Detector->local_width/Detector->local_sample_skip-1) * (Detector->local_height/Detector->local_sample_skip-1) );

	unsigned long long sky_sum = 0;
	unsigned int sky_count = 0;

	debug_printf("Detecting stars...\n");
	int i, j, x, y;
	for ( j = 0; j < Image->height; j += Detector->local_height)
	{
		for ( i = 0; i < Image->width; i += Detector->local_width)
		{
			LocalSample.clear();
			for ( y = j + Detector->local_sample_skip; y < j + Detector->local_height && y < Image->height; y += Detector->local_sample_skip )
				for ( x = i + Detector->local_sample_skip; x < i + Detector->local_width && x < Image->width; x += Detector->local_sample_skip )
					LocalSample.push_back( get_pixel(Image, x, y) );

			sort( LocalSample.begin(), LocalSample.end() );

			size_t N = LocalSample.size();
			short sample_min = LocalSample[0];
			short Q1 = LocalSample[ N/4 ];
			short median = LocalSample[ N/2 ];
			short Q3 = LocalSample[ N*3/4 ];

			int low_val = (int)( Q3 + 1.0f * ( Q3 - Q1 ) );
			sfloat outst_val = 0.6f * ( Q3 - Q1 );

#if 0
			printf("(%02d, %02d) min = %d, Q1 = %d, median = %d, Q3 = %d, max = %d, outliers at %d\n",
					i, j, sample_min, Q1, median, Q3, LocalSample[ N-1 ], low_val);
#endif

			// get star pixels and find centroid
			for ( y = j; y < j + Detector->local_height && y < Image->height; y++ )
				for ( x = i; x < i + Detector->local_width && x < Image->width; x++ )
				{
					if ( !isOutstanding(Image, x, y, low_val) )
					{
						sky_sum += get_pixel(Image, x, y);
						sky_count += 1;
					}
					else if ( !is_marked(Image, marked, x, y) && mean_over_area(Image, x, y, 1) >= Q3 + outst_val )
					{
						CentroidAt(Image, marked, ImageStars, low_val, sample_min, x, y);
					}
				}
		}
	}

	size_t ret = ImageStars.size();

#ifdef MAX_IMAGE_STARS
	vector<image_star>::iterator nth = ( ImageStars.size() > MAX_IMAGE_STARS ? ImageStars.begin() + MAX_IMAGE_STARS : ImageStars.end() );
	nth_element( ImageStars.begin(), nth, ImageStars.end(), compare_radius );
	ImageStars.erase(nth, ImageStars.end());
#endif

	timer.StopTimer();

	Detector->mean_sky = (sfloat) ((double)sky_sum / sky_count);

#ifdef DEBUG_TEXT
	vector<image_star>::iterator it;
	for ( it = ImageStars.begin(); it != ImageStars.end(); it++ )
	{
		printf("Image Star %03d: (%.3f, %.3f) | Error =  %.15f | r' = (%f,%f,%f)\n",
				int(it - ImageStars.begin())+1, it->centroid_x, it->centroid_y, it->error, it->r_prime[0], it->r_prime[1], it->r_prime[2] );
	}
#endif

	debug_printf("ImageStars successfully found, size: %d | Time elapsed: %d ms\n", ImageStars.size(), timer.GetTime());
#ifdef PROMPT_USER
	printf("Press ENTER to continue.\n");
	getchar();
#endif

	//return ImageStars;
	return ret;
}
