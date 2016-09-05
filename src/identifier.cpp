/* identifier.cpp
Alexander Wickes and Gil Tabak
September 20, 2010
*/

#include "identifier.h"
#include <hash_map>
#include <algorithm>
#include "timer.h"

using namespace std;

typedef unsigned int vote_t;

static bool compare_catalog(catalog_pair first, catalog_pair second)
{
	return (first.distance < second.distance);
}

/*	Map to keep track of votes for candidates stars for a given image star. Key is the index of the candidate catalog star,
	value is number of votes for the candidate star. Candidates not found in map are added with 1 vote, those found have their
	votes incremented by 1 */
typedef hash_map<index_t, vote_t> candidate_map;

// After filling the map structure, the data is passed to a vector for easy sorting by number of votes
typedef pair<index_t, vote_t> candidate_pair;
typedef vector<candidate_pair> candidate_vector;

static bool compare_map( candidate_map::value_type& left, candidate_map::value_type& right)
{
	return (left.second > right.second);
}

static bool compare_vector( candidate_pair& left, candidate_pair& right )
{
	return (left.second > right.second);
}

// Add a vote for a specific catalog star to a specific image star's map
static inline void add_votes_helper(candidate_map* StarCand, index_t img_star, index_t ctg_star)
{
	candidate_map::iterator it = StarCand[img_star].lower_bound(ctg_star);
	if ( it != StarCand[img_star].end() && it->first == ctg_star )
		it->second += 1;
	else
		StarCand[img_star].insert( it, candidate_map::value_type(ctg_star, 1) );
}

// Add votes for both catalog stars to both image stars' maps
static inline void add_votes(candidate_map* StarCand, index_t img_star1, index_t img_star2, index_t ctg_star1, index_t ctg_star2)
{
	add_votes_helper(StarCand, img_star1, ctg_star1);
	add_votes_helper(StarCand, img_star1, ctg_star2);
	add_votes_helper(StarCand, img_star2, ctg_star1);
	add_votes_helper(StarCand, img_star2, ctg_star2);
}

// Global vector of image pairs, to avoid allocating memory repeatedly
vector<image_pair> ImagePairs;

void GetImagePairs(vector<image_star>& ImageStars)
{
	ImagePairs.clear();

	//vector<image_pair> ImagePairs;
	ImagePairs.reserve( ImageStars.size() * (ImageStars.size()-1) / 2 );


	int i, j;
	// Iterate over all distinct pairs of image stars
	for ( i = 0; i < ImageStars.size()-1; i++ )
	{
		for ( j = i+1; j < ImageStars.size(); j++ )
		{
			sfloat ang_dist = acos( dot_product(ImageStars[i].r_prime, ImageStars[j].r_prime) );
			// Error is sum of errors for both stars
			sfloat error = ImageStars[i].error + ImageStars[j].error;
			//error *= .5f;

			image_pair pair = { i, j, cos(ang_dist + error), cos(ang_dist - error) };
			ImagePairs.push_back(pair);
		}
	}

	//return ImagePairs;
}

void IdentifyImageStars(vector<image_star>& ImageStars, Catalog& theCatalog, prior_s* Prior)
{
	Timer timer;
	timer.StartTimer();

	vector<catalog_star>& CatalogStars = theCatalog.Stars;
	vector<catalog_pair>& CatalogPairs = theCatalog.SortedPairs;

	//vector<image_pair> ImagePairs = GetImagePairs( ImageStars );
	ImagePairs.reserve( 1000 );
	GetImagePairs( ImageStars );

	sfloat rPrior[3];
	sfloat prior_dot;
	if ( Prior != NULL )
	{
		rPrior[0] = cos(Prior->coords.DEC) * cos(Prior->coords.RA);
		rPrior[1] = cos(Prior->coords.DEC) * sin(Prior->coords.RA);
		rPrior[2] = sin(Prior->coords.DEC);
		prior_dot = cos(Prior->phi_max);
	}

	candidate_map* StarCandidates = new candidate_map[ImageStars.size()];

	debug_printf("Identifying stars...\n");

	//fill data structure
	vector<image_pair>::iterator it;
	for ( it = ImagePairs.begin(); it != ImagePairs.end(); it++ )
	{
#if 0
		printf("ImagePair(%04d, %04d) | R = [%.15f, %.15f] = cos[%.15f, %.15f]\n", 
					it->star1, it->star2, it->lower, it->upper, acos(it->upper), acos(it->lower) );
#endif

		catalog_pair tmp = {
			INDEX_INVALID, INDEX_INVALID, it->lower
		};

		// do STL binary search for iterator to beginning of range
		vector<catalog_pair>::iterator it_cat = lower_bound( CatalogPairs.begin(), CatalogPairs.end(), tmp, compare_catalog );
		//printf( "lower bound, Pair(%d,%d), Dist = %f\n", it_cat->star1, it_cat->star2, it_cat->distance );
		for ( ; it_cat != CatalogPairs.end() && it_cat->distance <= it->upper; it_cat++ )
		{
			if ( Prior != NULL && (	dot_product( rPrior, CatalogStars[it_cat->star1].r) < prior_dot ||
									dot_product( rPrior, CatalogStars[it_cat->star2].r) < prior_dot ) )
				continue;

			/*
			printf("	Adding votes: ImagePair(%03d, %03d), CatalogPair(%04d, %04d) | Distance: %.15f\n",
					it->star1, it->star2, it_cat->star1, it_cat->star2, it_cat->distance);
			*/

			// for each CatalogStars pair in range, add vote for each CatalogStars star to each image star
			add_votes(StarCandidates, it->star1, it->star2, it_cat->star1, it_cat->star2); // Commenting this line fixes errors
		}
	}

	debug_printf( "Successfully populated candidate stuff, time: %d ms\n", timer.GetTime() );
	Timer timer2;
	timer2.StartTimer();

#if 0
	int i;
	for ( i = 0; i < ImageStars.size(); i++ )
	{
		printf("Image Star %d\n", i);
		candidate_map::iterator it;
		for ( it = StarCandidates1[i].begin(); it != StarCandidates1[i].end(); it++ )
		{
			printf("	1) CatalogStars ID: %d | Votes: %d\n", it->first, it->second);
		}
		for ( it = StarCandidates2[i].begin(); it != StarCandidates2[i].end(); it++ )
		{
			printf("	2) CatalogStars ID: %d | Votes: %d\n", it->first, it->second);
		}
	}
#endif

#ifdef ONLY_ONE_STAGE
	int i;
	for ( i = 0; i < ImageStars.size(); i++ )
	{
#ifdef PRINT_TOP_VOTES
		printf("%03d) (x,y) = (%f, %f) | Candidates: %d\n", i, ImageStars[i].centroid_x, ImageStars[i].centroid_y, StarCandidates[i].size() );

		candidate_vector tmp_vec(StarCandidates[i].begin(), StarCandidates[i].end());
		sort( tmp_vec.begin(), tmp_vec.end(), compare_vector );
		candidate_vector::iterator it;
		for ( it = tmp_vec.begin(); it != tmp_vec.end() && int(it - tmp_vec.begin()) < 10; it++ )
			printf("	%03d Votes | CtgStar: %04d | (RA,DEC) = (%.15f, %.15f)\n",
					it->second, it->first, CatalogStars[it->first].RA, CatalogStars[it->first].Dec );

		it = tmp_vec.begin();
		if ( it != tmp_vec.end() && it->second >= 0.75f*ImageStars.size() )
			ImageStars[i].identity = it->first;
#else
		candidate_map::iterator it = min_element(StarCandidates[i].begin(), StarCandidates[i].end(), compare_map );
		if ( it != StarCandidates[i].end() && it->second >= 0.75f*ImageStars.size() )
		{
			ImageStars[i].identity = it->first;
#if 0
			printf("	%03d Votes | CtgStar: %04d | (RA,DEC) = (%.15f, %.15f)\n",
					it->second, it->first, CatalogStars[it->first].RA, CatalogStars[it->first].Dec );
#endif
		}
#endif
	}
#else
	candidate_vector* candidates_vec = new candidate_vector[ImageStars.size()];

	//int low_votes = int( LOW_VOTES_COEFF*ImageStars.size() );

#define NUM_CANDIDATES 5
#define NUM_PARTIAL_SORT 20

	int i, j;
	for ( i = 0; i < ImageStars.size(); i++ )
	{
		candidate_vector tmp_vec( StarCandidates[i].begin(), StarCandidates[i].end() );
		candidate_vector::iterator it_mid = ( tmp_vec.size() > NUM_PARTIAL_SORT ? tmp_vec.begin()+NUM_PARTIAL_SORT : tmp_vec.end() );
		partial_sort( tmp_vec.begin(), it_mid, tmp_vec.end(), compare_vector );

		int min_votes = 1;
		candidate_vector::iterator it;
		for (	it = tmp_vec.begin();
				it != it_mid && ( int(it - tmp_vec.begin()) < NUM_CANDIDATES || it->second >= min_votes );
				it++ )
		{
			candidates_vec[i].push_back( candidate_pair(it->first, 0) );
			if ( int(it - tmp_vec.begin()) == NUM_CANDIDATES-1 )
				min_votes = it->second;
		}
	}

	for ( i = 0; i < ImageStars.size()-1; i++ )
	{
		for ( j = i+1; j < ImageStars.size(); j++ )
		{
			candidate_vector::iterator it1, it2;
			for ( it1 = candidates_vec[i].begin(); it1 != candidates_vec[i].end(); it1++ )
			{
				for ( it2 = candidates_vec[j].begin(); it2 != candidates_vec[j].end(); it2++ )
				{
					index_t pair_index = i*ImageStars.size() + j - (i+1)*(i+2)/2;
					sfloat dot = dot_product(CatalogStars[it1->first].r, CatalogStars[it2->first].r);
					sfloat lower = ImagePairs[pair_index].lower;
					sfloat upper = ImagePairs[pair_index].upper;
					if ( dot >= lower && dot <= upper )
					{
						// both stars good, add vote to current map
						it1->second += 1;
						it2->second += 1;
					}
				}
			}
		}
	}

	//low_votes = int(0.8f*low_votes);

	int max_votes = 0;

	for ( i = 0; i < ImageStars.size(); i++ )
	{
		candidate_vector::iterator it_mid = (	candidates_vec[i].size() > NUM_PARTIAL_SORT ?
												candidates_vec[i].begin()+NUM_PARTIAL_SORT : candidates_vec[i].end() );
#ifndef PRINT_TOP_VOTES
		// assign identities and exit
		candidate_vector::iterator it = min_element( candidates_vec[i].begin(), candidates_vec[i].end(), compare_vector );
#else
		partial_sort( candidates_vec[i].begin(), it_mid, candidates_vec[i].end(), compare_vector );
		candidate_vector::iterator it = candidates_vec[i].begin();
#endif

		if ( it != it_mid && it->second > max_votes )
			max_votes = it->second;
	}

	if ( true ) //max_votes <= 4 )
	{
		debug_printf("Low votes detected, doing another stage\n");
		for ( i = 0; i < ImageStars.size(); i++ )
		{
			candidate_vector::iterator it_mid = (	candidates_vec[i].size() > NUM_PARTIAL_SORT ?
				candidates_vec[i].begin()+NUM_PARTIAL_SORT : candidates_vec[i].end() );

			candidate_vector::iterator it;
			for (	it = candidates_vec[i].begin();
					it != it_mid && (	( ImageStars.size() >= IMAGE_STARS_USE_MAX_CUTOFF && it->second >= int( MAX_VOTES_COEFF * max_votes ) ) ||
										( ImageStars.size() < IMAGE_STARS_USE_MAX_CUTOFF && it->second > 0 ) );
					it++ );

			candidates_vec[i].erase(it, candidates_vec[i].end());

			for ( it = candidates_vec[i].begin(); it != candidates_vec[i].end(); it++ )
				it->second = 0;
		}

		for ( i = 0; i < ImageStars.size()-1; i++ )
		{
			for ( j = i+1; j < ImageStars.size(); j++ )
			{
				candidate_vector::iterator it1, it2;
				for ( it1 = candidates_vec[i].begin(); it1 != candidates_vec[i].end(); it1++ )
				{
					for ( it2 = candidates_vec[j].begin(); it2 != candidates_vec[j].end(); it2++ )
					{
						index_t pair_index = i*ImageStars.size() + j - (i+1)*(i+2)/2;
						sfloat dot = dot_product(CatalogStars[it1->first].r, CatalogStars[it2->first].r);
						sfloat lower = ImagePairs[pair_index].lower;
						sfloat upper = ImagePairs[pair_index].upper;
						if ( dot >= lower && dot <= upper )
						{
							// both stars good, add vote to current map
							it1->second += 1;
							it2->second += 1;
						}
					}
				}
			}
		}

		max_votes = 0;

		for ( i = 0; i < ImageStars.size(); i++ )
		{
			candidate_vector::iterator it_mid = (	candidates_vec[i].size() > NUM_PARTIAL_SORT ?
				candidates_vec[i].begin()+NUM_PARTIAL_SORT : candidates_vec[i].end() );

#ifndef PRINT_TOP_VOTES
			// assign identities and exit
			candidate_vector::iterator it = min_element( candidates_vec[i].begin(), candidates_vec[i].end(), compare_vector );
#else
			partial_sort( candidates_vec[i].begin(), it_mid, candidates_vec[i].end(), compare_vector );
			candidate_vector::iterator it = candidates_vec[i].begin();
#endif

			if ( it != it_mid && it->second >= max_votes )
				max_votes = it->second;
		}
	}

	int low_votes = int( MAX_VOTES_COEFF * max_votes );

	debug_printf("Max Votes: %d | Low Votes: %d\n", max_votes, low_votes );

	for ( i = 0; i < ImageStars.size(); i++ )
	{
#ifndef PRINT_TOP_VOTES
		// assign identities and exit
		candidate_vector::iterator it = min_element( candidates_vec[i].begin(), candidates_vec[i].end(), compare_vector );
#else
		candidate_vector::iterator it = candidates_vec[i].begin();
#endif

		if ( low_votes >= 3 && it != candidates_vec[i].end() && it->second >= low_votes )
			ImageStars[i].identity = it->first;
		else
			ImageStars[i].identity = INDEX_INVALID;

		sfloat RA = 0.0f, DEC = 0.0f;
		if (ImageStars[i].identity != INDEX_INVALID && ImageStars[i].identity != INDEX_FALSE_STAR)
		{
			RA = CatalogStars[ImageStars[i].identity].RA;
			DEC = CatalogStars[ImageStars[i].identity].Dec;
		}

#ifdef DEBUG_TEXT
#ifdef PRINT_TOP_VOTES
		printf("(%02d) Candidates: %d | Identity: %d | Votes: %d | (X,Y) = (%.3f,%.3f) | (RA,DEC) = (%f,%f)\n",
			i+1, candidates_vec[i].size(), ImageStars[i].identity, it->second,
			ImageStars[i].centroid_x, ImageStars[i].centroid_y, RA, DEC );

		for ( j = 0; it != candidates_vec[i].end() && j < 20; it++, j++ )
		{
			printf("	%04d Votes | CtgStar: %04d | (RA,DEC) = (%f, %f)\n",
				it->second, it->first, CatalogStars[it->first].RA, CatalogStars[it->first].Dec );
		}
#endif /* PRINT_TOP_VOTES */
#endif /* DEBUG_TEXT */
	}

	delete[] candidates_vec;
#endif
	delete[] StarCandidates;

	timer2.StopTimer();
	timer.StopTimer();

	debug_printf("Successfully identified stars. | Time elapsed: %d ms | Total Time: %d ms\n", timer2.GetTime(), timer.GetTime());
#ifdef PROMPT_USER
	printf("Press ENTER to continue.\n");
	getchar();
#endif
}