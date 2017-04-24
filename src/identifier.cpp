/* identifier.cpp
Alexander Wickes and Gil Tabak
September 20, 2010
*/

#include "identifier.h"
#include <hash_map>
#include <algorithm>
#include "timer.h"

#define NUM_CANDIDATES		5		//
#define NUM_PARTIAL_SORT	20		//max number of candidates, modulo a tie in votes
#define NUM_TO_USE			15		//number of image stars to use

using namespace std;

typedef unsigned int vote_t;

/*	Map to keep track of votes for candidates stars for a given image star. Key is the index of the candidate catalog star,
	value is number of votes for the candidate star. Candidates not found in map are added with 1 vote, those found have their
	votes incremented by 1 */
typedef hash_map<index_t, vote_t> candidate_map;

// After filling the map structure, the data is passed to a vector for easy sorting by number of votes
// define a data type for candidates, consisting of index of catalog star and number of votes. same type as used internally by hash_map
typedef pair<index_t, vote_t> candidate_t;
// define a data type for vectors of candidates, for convenience
typedef vector<candidate_t> candidate_vector;
// Global vector of image pairs, to avoid allocating memory repeatedly
vector<image_pair> ImagePairs;

static bool compare_catalog(catalog_pair first, catalog_pair second)
{
	return (first.distance < second.distance);
}

static bool compare_map( candidate_map::value_type& left, candidate_map::value_type& right)
{
	return (left.second > right.second);
}

static bool compare_vector( candidate_t& left, candidate_t& right )
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

static inline void PurgeImageStars(vector<image_star>& ImageStars){
	//if the number of ImageStars is less then or equal to the given number, do nothing
    if (ImageStars.size() <= NUM_TO_USE)
        return;
    //sort the image stars by error (use reverse iterators to obtain descending order)
    sort(ImageStars.rbegin(), ImageStars.rend() );
   
    ImageStars.erase (ImageStars.begin()+NUM_TO_USE,ImageStars.end() );
}

// fills the ImagePairs vector with image_pair structs constructed from the image stars
void GetImagePairs(vector<image_star>& ImageStars){

	ImagePairs.clear();

	//vector<image_pair> ImagePairs;
	ImagePairs.reserve( ImageStars.size() * (ImageStars.size()-1) / 2 );

	// Iterate over all distinct pairs of image stars
	for (int i = 0; i < ImageStars.size()-1; i++ )
	{
		for (int j = i+1; j < ImageStars.size(); j++ )
		{
			sfloat ang_dist = acos( dot_product(ImageStars[i].r_prime, ImageStars[j].r_prime) );
			// Use sum of errors to make sure we include all stars within range.
			sfloat error = ImageStars[i].error + ImageStars[j].error;
			image_pair pair = { i, j, cos(ang_dist + error), cos(ang_dist - error) };
			ImagePairs.push_back(pair);
		}
	}
	//return ImagePairs;
}


// uncomment to print out debug information when running the method
//#define DEBUG_IDENTIFY_METHOD

// takes a vector of image stars, a catalog, and a prior (stores information about previously acquired attitude)
void IdentifyImageStars(vector<image_star>& ImageStars, Catalog& theCatalog, prior_s* Prior){

	// declare and start a timer for measuring execution speed
	Timer timer;
	timer.StartTimer();

	// restrict the number of stars used. Ordered based on signal (can modify in star.h)
	PurgeImageStars(ImageStars);

	// create local references to the catalog stars vector and catalog pairs vector. this is purely for convenience.
	vector<catalog_star>& CatalogStars = theCatalog.Stars;
	vector<catalog_pair>& CatalogPairs = theCatalog.SortedPairs;

	// requests that the vector capacity be enough to contain n=1000 elements
	GetImagePairs( ImageStars );

	// create a new array of candidate hash maps for each candidate star
	candidate_map* StarCandidates = new candidate_map[ImageStars.size()];

	debug_printf("Identifying stars...\n");


	// INITIAL VOTING STAGE: Populating votes
	// for each image pair, go through range of possible catalog pairs. for each catalog pair with both stars in range of
	// prior (or without a prior), add a vote for catalog stars to each image star in image pair

	// if prior data is valid, construct the prior vector and compute the cosine of phi_max for later use
	sfloat rPrior[3];
	sfloat prior_dot;
	if ( Prior != NULL )
	{
		rPrior[0] = cos(Prior->coords.DEC) * cos(Prior->coords.RA);
		rPrior[1] = cos(Prior->coords.DEC) * sin(Prior->coords.RA);
		rPrior[2] = sin(Prior->coords.DEC);
		prior_dot = cos(Prior->phi_max);
	}

	vector<image_pair>::iterator it;
	for ( it = ImagePairs.begin(); it != ImagePairs.end(); it++ )
	{
#ifdef DEBUG_IDENTIFY_METHOD
		printf("ImagePair(%04d, %04d) | R = [%.15f, %.15f] = cos[%.15f, %.15f]\n", 
					it->star1, it->star2, it->lower, it->upper, acos(it->upper), acos(it->lower) );
#endif

		// create a temporary catalog_pair with current ImagePair's lower bound distance, for the binary search that follows
		catalog_pair tmp = {
			INDEX_INVALID, INDEX_INVALID, it->lower
		};
		
		// do STL binary search for iterator pointing to beginning of range specified by the current ImagePair.
		// comparison function compare_catalog (specified at beginning of this file) just returns the boolean value (first.distance < second.distance)
		vector<catalog_pair>::iterator it_cat = lower_bound( CatalogPairs.begin(), CatalogPairs.end(), tmp, compare_catalog );

		/*
		// legacy debug printout
		printf( "lower bound, Pair(%d,%d), Dist = %f\n", it_cat->star1, it_cat->star2, it_cat->distance );
		*/

		// for loop iterator initialized above to catalog_pair iterator at beginning of range specified by the ImagePair.
		// stops at iterator at end of range specified by the ImagePair, or at last catalog_pair
		for (/*initialized above*/; it_cat != CatalogPairs.end() && it_cat->distance <= it->upper; it_cat++ )
		{
			if ( Prior != NULL && (	dot_product( rPrior, CatalogStars[it_cat->star1].r) < prior_dot ||
									dot_product( rPrior, CatalogStars[it_cat->star2].r) < prior_dot ) )
				// if either catalog star is out of range of prior, skip everything after this and go on to next loop iteration
				continue;

			/*
			// legacy debug printout
			printf("	Adding votes: ImagePair(%03d, %03d), CatalogPair(%04d, %04d) | Distance: %.15f\n",
					it->star1, it->star2, it_cat->star1, it_cat->star2, it_cat->distance);
			*/

			// add a vote for each catalog star in the catalog_pair to each image star in the image_pair
			add_votes(StarCandidates, it->star1, it->star2, it_cat->star1, it_cat->star2);
		}
	}

	// print out time taken up to populating votes in initial voting stage
	debug_printf( "Successfully populated candidate stuff, time: %d ms\n", timer.GetTime() );

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// INITIAL VOTING STAGE: 
	// create another timer for timing 
	Timer timer2;
	timer2.StartTimer();

	// create candidate_vector (vector of candidate_t objects) for each image star, to be used for NEXT voting stage
	candidate_vector* candidates_vec = new candidate_vector[ImageStars.size()];


	// declare loop counters
	int i, j;
	// iterate over each image star
	for ( i = 0; i < ImageStars.size(); i++ )
	{
		// create a temporary candidate_vector filled with all candidates from the corresponding hash_map
		// note that the OBJECT created by hash_map<index_t,vote_t> is exactly of type candidate_t, that is pair<index_t,vote_t>
		candidate_vector tmp_vec( StarCandidates[i].begin(), StarCandidates[i].end() );

		// create an iterator to the NUM_PARTIAL_SORTth element of tmp_vec if number of candidate pairs exceeds NUM_PARTIAL_SORT, else to the last element
		candidate_vector::iterator it_stop = ( tmp_vec.size() > NUM_PARTIAL_SORT ? tmp_vec.begin()+NUM_PARTIAL_SORT : tmp_vec.end() );

		// do a partial (quick) sort of tmp_vec so that the N = [it_mid - tmp_vec.begin()] candidates with the highest votes
		// are at the beginning (in descending order? should be but not sure anymore if that's actually true)
		partial_sort( tmp_vec.begin(), it_stop, tmp_vec.end(), compare_vector );

		// declare iteration variables
		int min_votes = 1;
		candidate_vector::iterator it;
		// iterates until the NUM_CANDIDATESth candidate, and then until the next candidate with properly fewer votes
		// this is so that all candidates with this lowest number of votes is included
		for (	it = tmp_vec.begin();
				it != it_stop && ( int(it - tmp_vec.begin()) < NUM_CANDIDATES || it->second >= min_votes );
				it++ )
		{
			// add a copy of the candidate star with no votes to the end of the candidates_vec vector for image star i
			candidates_vec[i].push_back( candidate_t(it->first, 0) );
			if ( int(it - tmp_vec.begin()) == NUM_CANDIDATES-1 )
				// just sets min_votes according to the behavior described immediately before the for loop
				min_votes = it->second;
		}
	}

	
	// VOTING STAGES 2 & 3
	// define variables for later use
	candidate_vector::iterator it1, it2;
	int low_votes = 0;
	int max_votes = 0;
	int stage;
	for (stage = 2; stage <= 3; stage++)
	{
		// only do this in stage 3. If you wanted to generalize this algorithm to larger data sets and more stages,
		// you could do this in each stage with the subset of candidates deleted depending functionally on the stage 
		if ( stage == 3 )
		{
			// for each image star, find the last candidate star leading it_stop with too few votes and erase that one and all that follow it
			for ( i = 0; i < ImageStars.size(); i++ )
			{
				// generate end iterator as before
				candidate_vector::iterator it_stop = (	candidates_vec[i].size() > NUM_PARTIAL_SORT ?
					candidates_vec[i].begin()+NUM_PARTIAL_SORT : candidates_vec[i].end() );

				// iterate until either it_stop, or there are too many image stars and we find one with < low_votes, or there are not enough image stars and we find one with 0 votes
				candidate_vector::iterator it;
				for (	it = candidates_vec[i].begin();
						it != it_stop && (	( ImageStars.size() >= IMAGE_STARS_USE_MAX_CUTOFF && it->second >= low_votes ) ||
											( ImageStars.size() < IMAGE_STARS_USE_MAX_CUTOFF && it->second > 0 ) );
						it++ );
				
				// throw away all candidates for the image star starting with the one just found and including all with fewer votes
				candidates_vec[i].erase(it, candidates_vec[i].end());

				// reset all votes to 0 for stage 3 voting
				for ( it = candidates_vec[i].begin(); it != candidates_vec[i].end(); it++ )
					it->second = 0;
			}
		}

		// iterate over distinct image star pairs i<j
		for ( i = 0; i < ImageStars.size()-1; i++ )
		{
			for ( j = i+1; j < ImageStars.size(); j++ )
			{
				// iterate over all possible candidate pairs from the image pair
				for ( it1 = candidates_vec[i].begin(); it1 != candidates_vec[i].end(); it1++ )
				{
					for ( it2 = candidates_vec[j].begin(); it2 != candidates_vec[j].end(); it2++ )
					{
						// generate the image pair index corresponding to i,j. the equation comes from the fact that the pairs 
						// are those in the upper triangle of the ij matrix (i<j) written one after the other, so we take take the
						// linearized ij index and subtract off the positions skipped in the lower triangle.
						// I may be wrong, but there seems to have been a bug here. was using (i+1)(i+2)/2 instead of i*(i+1)/2.
						// Gil: new values did not seem to work. Changed back to old values.
						index_t pair_index = i*ImageStars.size() + j - (i+1)*(i+2)/2; 
						// compute cosine of the angle between the two candidate stars, for use in comparison (below)
						sfloat dot = dot_product(CatalogStars[it1->first].r, CatalogStars[it2->first].r);
						// just renaming lower and upper bounds on the image pair dot product range, for convenience
						sfloat lower = ImagePairs[pair_index].lower;
						sfloat upper = ImagePairs[pair_index].upper;
						// see if the dot product is in the correct range
						if ( dot >= lower && dot <= upper )
						{
							// if so, add votes to current candidates
							it1->second += 1;
							it2->second += 1;
						}
					}
				}
			}
		}

		// find the highest number of votes out of all candidate stars for all image stars, put it in max_votes
		max_votes = 0;
		for ( i = 0; i < ImageStars.size(); i++ )
		{
			// generate end iterator like before
			candidate_vector::iterator it_stop = (	candidates_vec[i].size() > NUM_PARTIAL_SORT ?
													candidates_vec[i].begin()+NUM_PARTIAL_SORT : candidates_vec[i].end() );

			// partial/quick sort like before
			partial_sort( candidates_vec[i].begin(), it_stop, candidates_vec[i].end(), compare_vector );
		
			// if first iterator is different from end iterator, and its vote count is higher, update max_votes with its value
			candidate_vector::iterator it = candidates_vec[i].begin();
			if ( it != it_stop && it->second > max_votes )
				max_votes = it->second;
		}
	
		// vote bound for use in stage 3 and later
		low_votes = int( MAX_VOTES_COEFF * max_votes );

		debug_printf("Max Votes: %d | Low Votes: %d\n", max_votes, low_votes );
	}


	// for each image star, 
	// if there the highest vote count is too small, we don't have enough data so assign all image stars with invalid identities,
	// else if there are valid candidates with enough votes, assign the image star with identity of its candidate star with most votes
	for ( i = 0; i < ImageStars.size(); i++ )
	{
		// set iterator to the candidate we are interested in, the one that bubbled up to the top with the most votes after all stages
		candidate_vector::iterator it = candidates_vec[i].begin();

		// if the max vote count is sufficiently large, and there is a surviving candidate, and the best candidate has sufficiently many votes,
		// then assign it with that candidate's identity. otherwise, assign it with an invalid identity.
		if ( low_votes >= 3 && it != candidates_vec[i].end() && it->second >= low_votes )
			ImageStars[i].identity = it->first;
		else
			ImageStars[i].identity = INDEX_INVALID;

		// if the image star now has a valid identity, give it the corresponding catalog star's coordinates
		sfloat RA = 0.0f, DEC = 0.0f;
		if (ImageStars[i].identity != INDEX_INVALID && ImageStars[i].identity != INDEX_FALSE_STAR)
		{
			RA = CatalogStars[ImageStars[i].identity].RA;
			DEC = CatalogStars[ImageStars[i].identity].Dec;
		}
	}

	// clean up data structures. this is probably inefficient and could/should be optimized!
	delete[] candidates_vec;
	delete[] StarCandidates;
	
#ifdef TIMERS_ON
	// stop both timers
	timer2.StopTimer();
	timer.StopTimer();
	printf("Successfully identified stars. | Time for ID algorithm: %d ms | Total Time for ID method: %d ms\n", timer2.GetTime(), timer.GetTime());
#endif
#ifdef PROMPT_USER
	printf("Press ENTER to continue.\n");
	getchar();
#endif
}
