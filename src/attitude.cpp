/* attitude.cpp
Alexander Wickes and Gil Tabak
September 20, 2010
*/

#include "attitude.h"
#include <cstdio>
#include <cmath>
#include "timer.h"

using namespace std;

#define maxEigenError 0.000000001 							//max error allowed on the rotational quaternion

#define DELTA_ERROR 0.000000001

sfloat eigen_vector[4];

void quart(double a[5], double rr[4], double ri[4]); 		// quartic solver
void computeNco(double N[4][4], double Nco[4][4]); 			//computes the cofactor of the N matrix (puts it in Nco)
double reOrient(double eVector[4], double N[4][4]); 		//iterative method for fixing any errors caused by floating point

void WeighByError( vector<image_star>& ImageStars, vector<catalog_star>& Catalog, double R[3][3], double rOut[3] ); //not tested extensively.  Might not help much.

/*

---------------------------------------------------------------------------------------------------------------------------------------

GetAttitude takes a set of coordinates from the catalog (as catalog stars) and the corresponding coordinates of the image stars,
calculated from RA, DEC, and focal length.  The rotational matrix and center of the image coordinates are determined.

The first method used is described thoroughly here: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.65.971&rep=rep1&type=pdf
Find the quaternion which is an eigenvector to a certain matrix (corresponding to the largest eigenValue).
Follow the methods described to determine the matrix, compute its eigenValues by solving the
corresponding polynomial (a quartic with four real solutions), and determine the corresponding eigenVector
by the method described in Appendix 5 of the paper.

There were small errors in the solution, probably due to floating point. To make 
up for that, use an iterative method to sharpen the precision of the eigenVector determined.

Once this was all done, the quaternion was converted to a rotation matrix, and finally the coordinates of the center were determined.

---------------------------------------------------------------------------------------------------------------------------------------

The S matrix is a matrix of the sum of the dot products of the coordinates of stars
in the two reference systems used.

eValue is the largest eigenValue, and eVector is its corresponding eigenVector.

N is a matrix determined from S, as specified in the paper by Horn (Closed-form solution of absolute orientation using unit quaternions)
We determine the eigenvalue of N by solving its characteristic quartic a[4].

R is the rotation matrix
*/

double GetAttitude(	vector<image_star>& ImageStars, vector<catalog_star>& CatalogStars, 
					coordinates* BoresightOut, coordinates* CoordOut, sfloat RMat[3][3], double RQuat[4])
{
	Timer timer;
	timer.StartTimer();

	double R[3][3];				//the rotational matrix to be determined.
	double a[5], rr[4], ri[4];	//a[0] = highest coefficient for polynomial, a[4] = constant coefficient
								//rr are the real parts of the solutions, ri are the imaginary parts
								//there should be no imaginary solution due to the construction of the quartic
								//(eigenvalues of an orthogonal matrix)

	//calculate the S matrix
	double S[3][3] = {	{ 0, 0, 0 },
						{ 0, 0, 0 },
						{ 0, 0, 0 } };
	int valid_stars = 0;
	int i, j;
	vector<image_star>::iterator it;
	for ( it = ImageStars.begin(); it != ImageStars.end(); it++ )					//iterate over image stars
	{
		if ( it->identity != INDEX_INVALID && it->identity != INDEX_FALSE_STAR )	//valid star
		{
			for ( i = 0; i < 3; i++ )
				for ( j = 0; j < 3; j++ )
					S[i][j] += CatalogStars[it->identity].r[i] * it->r_prime[j];	//include its info in the S matrix
			valid_stars += 1;
		}

	}

	//find the a[2] coefficient
	a[2] = 0;
	for (i = 0; i <3; i++)
		for (j = 0; j <3; j++)
		{
			a[2] += -2 * S[i][j] * S[i][j];
			//printf("S%d%d = %f\n",i,j,S[i][j]);
		}

	//the other a[] coefficients
	a[0] = 1;
	a[1] = 0;
	a[3] = 8 * (S[0][0]*S[1][2]*S[2][1]+S[1][1]*S[2][0]*S[0][2]
				+ S[2][2]*S[0][1]*S[1][0] - S[0][0]*S[1][1]*S[2][2]
				- S[1][2]*S[2][0]*S[0][1] - S[2][1]*S[1][0]*S[0][2]);

	//compute the N matrix from S
	double N[4][4];
	N[0][0] = S[0][0] + S[1][1] + S[2][2];
	N[1][0] = N[0][1] = S[1][2] - S[2][1];
	N[2][0] = N[0][2] = S[2][0] - S[0][2];
	N[3][0] = N[0][3] = S[0][1] - S[1][0];
	N[1][1] = S[0][0] - S[1][1] - S[2][2];
	N[2][1] = N[1][2] = S[0][1] + S[1][0];
	N[3][1] = N[1][3] = S[2][0] + S[0][2];
	N[2][2] = -S[0][0] + S[1][1] - S[2][2];
	N[3][2] = N[2][3] = S[1][2] + S[2][1];
	N[3][3] = -S[0][0] - S[1][1] + S[2][2];

	//printf("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f",N[0][0],N[0][1],N[0][2],N[0][3],N[1][0],N[1][1],N[1][2],N[1][3],N[2][0],N[2][1],N[2][2],N[2][3],N[3][0],N[3][1],N[3][2],N[3][3]);

	//compute the final a[] coefficient
	a[4] = N[0][0] * ( N[1][1]*N[2][2]*N[3][3] + N[2][1]*N[3][2]*N[3][1]
					 + N[3][1]*N[2][1]*N[3][2] - N[1][1]*N[3][2]*N[3][2]
					 - N[2][1]*N[2][1]*N[3][3] - N[3][1]*N[2][2]*N[3][1] )
		 - N[1][0] * ( N[2][1]*N[3][2]*N[3][0] + N[3][1]*N[2][0]*N[3][2]
					 + N[1][0]*N[2][2]*N[3][3] - N[3][2]*N[3][2]*N[1][0]
					 - N[3][3]*N[2][0]*N[2][1] - N[3][0]*N[2][2]*N[3][1] )
		 + N[2][0] * ( N[3][1]*N[2][0]*N[3][1] + N[1][0]*N[2][1]*N[3][3]
					 + N[1][1]*N[3][2]*N[3][0] - N[3][3]*N[2][0]*N[1][1]
					 - N[3][0]*N[2][1]*N[3][1] - N[3][1]*N[3][2]*N[1][0] )
		 - N[3][0] * ( N[1][0]*N[2][1]*N[3][2] + N[1][1]*N[2][2]*N[3][0]
					 + N[2][1]*N[2][0]*N[3][1] - N[3][0]*N[2][1]*N[2][1]
					 - N[3][1]*N[2][2]*N[1][0] - N[3][2]*N[2][0]*N[1][1] );

	//printf("%fx^4+%fx^3+%fx^2+%fx+%f = 0\n",a[0],a[1],a[2],a[3],a[4]);
	
	//solve the quartic
	quart(a, rr, ri);

	//determine the largest eigenValue
	double eValue = MAX(rr[0], MAX(rr[1], MAX(rr[2], rr[3])));

	//subtract the largest eigenValue from the main diagonal of the N matrix
	for (i = 0; i < 4; i++)
		N[i][i] -= eValue;

	//compute the adjacent matrix co-factor of N
	double Nco[4][4];
	computeNco(N, Nco);

	//compute the eigenVector of N from the adjacent matrix co-factor
	double eVector[4] = {0, 0, 0, 0};
	for ( i = 0; i < 4; i++ )
		for ( j = 0; j < 4; j++ )
			eVector[i] += Nco[i][j];

	//fix the N matrix to how it was before
	for ( i = 0; i < 4; i++ )
		N[i][i] += eValue;

	//this fixes any floating-point errors in the eigenVector using a recursive method
	eValue = reOrient(eVector, N);

	for ( i = 0; i < 4; i++ )
		eigen_vector[i] = eVector[i];

	//calculate the error from the eigenvalue
	double rms_error = sqrt( 2 * (1.0 - eValue / valid_stars) );

	debug_printf( "eValue: %.15f | Error: %.15f\n", eValue, rms_error );

	debug_printf("eigenVector: %f %f %f %f\n", eVector[0], eVector[1], eVector[2], eVector[3] );

	debug_printf("eigenvalues: %f %f %f %f\n", rr[0], rr[1], rr[2], rr[3] );

	//compute the rotation matrix from the quaternion
	GetMatrixFromQuat(R, eVector);

	double rCenter[3], rBoresightVector[3];
	double Rt[3][3];
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			Rt[i][j] = R[j][i];						//compute the transpose of R, called Rt
#ifdef WEIGH_STAR_ERROR	
	WeighByError( ImageStars, CatalogStars, Rt, rCenter );
#else
	for ( i = 0; i < 3; i++ )
		rCenter[i] = R[0][i];						//compute the coordiantes of the center of the image
	for ( i = 0; i < 3; i++ )
		rBoresightVector[i] = R[1][i];				//the boresight vector is used to determine the boresight angle (rotation about the centeral axis).
	for ( it = ImageStars.begin(); it != ImageStars.end(); it++ )
		if (it->identity != INDEX_INVALID)
			matrix_mult(Rt, it->r_prime, it->r);	//
#endif

	//compute the RA and DEC of the center
	double tmpRA, tmpDec;
	tmpRA = atan2( rCenter[1], rCenter[0] );
	if ( tmpRA < 0 )
		tmpRA += 2 * M_PI;
	tmpDec = atan( rCenter[2] / sqrt( SQUARE(rCenter[0]) + SQUARE(rCenter[1]) ) );

	//print the RA and DEC.
	debug_printf(	"(RA,DEC) = (%.15f, %.15f)\n", tmpRA, tmpDec);
	debug_printf(	"RA: %02d:%02d:%f DEC: %02d:%02d:%f\n",
					int( tmpRA/(2*M_PI)*24 ), int( tmpRA/(2*M_PI)*24*60 ) % 60, fmod( tmpRA/(2*M_PI)*24*60*60, 60 ),
					int( tmpDec*(180/M_PI) ), int( abs(tmpDec)*(180/M_PI)*60 ) % 60, fmod( abs(tmpDec)*(180/M_PI)*60*60, 60) );

	if ( CoordOut != NULL )
	{
		CoordOut->RA = (sfloat)tmpRA;
		CoordOut->DEC = (sfloat)tmpDec;
	}

	//compute the RA and DEC of the boresight vector
	tmpRA = atan2( rBoresightVector[1], rBoresightVector[0] );
	if ( tmpRA < 0 )
		tmpRA += 2 * M_PI;
	tmpDec = atan( rBoresightVector[2] / sqrt( SQUARE(rBoresightVector[0]) + SQUARE(rBoresightVector[1]) ) );

	if ( BoresightOut != NULL )
	{
		BoresightOut->RA = (sfloat)tmpRA;
		BoresightOut->DEC = (sfloat)tmpDec;
	}

	if ( RMat != NULL )
	{
		//print the R matrix
		for ( i = 0; i < 3; i++ )
			for ( j = 0; j < 3; j++ )
				RMat[i][j] = (sfloat)R[i][j];
	}

	if ( RQuat != NULL )
	{
		for ( i = 0; i < 3; i++ )
			RQuat[i] = eVector[i];
	}

	timer.StopTimer();

	debug_printf(	"Successfully estimated attitude. | Time elapsed: %d ms\n", timer.GetTime());

	return rms_error;
}

//used by quart
void quadra(double a1, double a2, double a3, double *rr1, double *rr2, double *ri1, double *ri2)
{

	double rad, srad;

	rad = a2 * a2 - 4 * a1 * a3;
	if ( rad >= 0 )
	{
		srad = sqrt( rad );
		*rr1 = (-a2 - srad) / (2 * a1);
		*rr2 = (-a2 + srad) / (2 * a1);
		*ri1 = 0;
		*ri2 = 0;
		return;
	}
	else
	{
		srad = sqrt( -rad );
		*rr1 = -a2 / (2 * a1);
		*rr2 = * rr1;
		*ri1 = srad/ (2 * a1);
		*ri2 = - *ri1;
		return;
	}

}

//used by quart
void cubic(double a[4], double rr[3], double ri[3])
{
	int i;
	double a0, a1, a2, a3;
	double g, h, y1, sh, xk, theta, xy1, xy2, xy3;
	double y2, z1, z2, z3, z4;

	for (i = 0; i < 3; i ++)
	{
		rr[i] = 0;
		ri[i] = 0;
	}

	a0 = a[0];
	a1 = a[1] / 3;
	a2 = a[2] / 3;
	a3 = a[3];

	g = (a0 * a0) * a3 - 3 * a0 * a1 * a2 + 2 * pow(a1, 3);
	h = a0 * a2 - a1 * a1;
	y1 = g * g + 4 * pow(h, 3);

	if (y1 < 0)
	{
		sh = sqrt(-h);
		xk = 2 * sh;
		theta = acos(g / (2 * h * sh)) / 3;
		xy1 = 2 * sh * cos(theta);
		xy2 = 2 * sh * cos(theta + (2 * M_PI / 3));
		xy3 = 2 * sh * cos(theta + (4 * M_PI / 3));
		rr[0] = (xy1 - a1) / a0;
		rr[1] = (xy2 - a1) / a0;
		rr[2] = (xy3 - a1) / a0;
		return;
	}
	else
	{
		y2 = sqrt(y1);
		z1 = (g + y2) / 2;
		z2 = (g - y2) / 2;
		if (z1 < 0)
		{
			z3 = pow(-z1, 1/3);
			z3 = -z3;
		}
		else 
			z3 = pow(z1, 1/3);
		if (z2 < 0)
		{
			z4 = pow(-z2, 1/3);
			z4 = - z4;
		}
		else
			z4 = pow(z2, 1/3);

		rr[0] = -(a1 + z3 + z4) / a0;
		rr[1] = (-2 * a1 + z3 + z4) / (2 * a0);
		ri[1] = sqrt(3.0f) * (z4 - z3) / (2 * a0);
		rr[2] = rr[1];
		ri[2] = -ri[1];

		return;

	}
}

//compute the roots for the quartic. There should be 4 real solutions due to the nature
//of the matrix (orthogonal).
void quart(double a[5], double rr[4], double ri[4])
{
	int i;
	double aa[5], b[4], rrc[3], ric[3];
	double x, y, z, c1, c2, c3, qr1, qr2, qi1, qi2;

	aa[0] = a[0];
	for (i = 1; i < 5; i ++)
		aa[i] = a[i] / a[0];

	b[0] = 1;
	b[1] = -aa[2];
	b[2] = aa[3] * aa[1] - 4 * aa[4];
	b[3] = aa[4] * (4 * aa[2] - aa[1] * aa[1]) - aa[3] * aa[3];

	cubic(b, rrc, ric);

	if (fabs(ric[1]) < 1.0e-6 ) // ric[1] == 0
	{
		x = MAX(rrc[0], MAX(rrc[1], rrc[2]));
		rrc[0] = x;
	}
	x = rrc[0] / 2;
	if ((x * x - aa[4]) > 0)
	{
		y = sqrt(x * x - aa[4]);
		z = - (aa[3] - aa[1] * x) / (2 * y);
	}
	else
	{
		y = 0;
		z = sqrt( pow(aa[1] / 2, 2) + 2 * x - aa[2] );
	}

	c1 = 1;
	c2 = aa[1] / 2 + z;
	c3 = x + y;

	quadra(c1, c2, c3, &qr1, &qr2, &qi1, &qi2);

	rr[0] = qr1;
	rr[1] = qr2;
	ri[0] = qi1;
	ri[1] = qi2;

	c1 = 1;
	c2 = aa[1] / 2 - z;
	c3 = x - y;

	quadra(c1, c2, c3, &qr1, &qr2, &qi1, &qi2);
	rr[2] = qr1;
	rr[3] = qr2;
	ri[2] = qi1;
	ri[3] = qi2;

	return;

}

//gives the co-factor of N at coordinate i,j
double coFactor(double N[4][4], int i, int j)
{
	double M[3][3]; //the matrix whose determinant is found
	int xNew,yNew; //coordinates for M
	bool odd;
	if ((i+j)%2 != 0)
		odd = true;
	else
		odd = false; //multiply by -1 if i+j is odd

	for (int x = 0; x < 4; x++)
		for (int y = 0; y < 4; y++) //x and y are the coordinates of N
			if (x != i && y != j) //skip the parts where x == i or y == j
			{
				if (x > i) xNew = x - 1; else xNew = x; 
				if (y > j) yNew = y - 1; else yNew = y; //find the coordinates of M
				M[xNew][yNew] = N[x][y]; //take values from N
				if (odd) M[xNew][yNew] *= -1; //change sign if necessary
			}

	return	  M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1]
			- M[2][0]*M[1][1]*M[0][2] - M[2][1]*M[1][2]*M[0][0] - M[2][2]*M[1][0]*M[0][1]; //determinant of M
}

//compute the adjacent co-factor matrix of N
void computeNco(double N[4][4], double Nco[4][4])
{
	for (int i = 0; i <4; i++)
	{
		//printf("\n");
		for (int j = 0; j <4; j++)
		{
			Nco[i][j] = coFactor(N, i,j);
			//printf("%f ", Nco[i][j]);
		}
	}
}

static inline double magnitude4(double vector[4])
{
	return sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]+vector[3]*vector[3]);
}

//normalized eVector
static inline void normalize4(double vector[4])
{
	double mag = magnitude4(vector);
	if ( abs(mag) < 0.0000001 )
		return;

	int i;
	for ( i = 0; i < 4; i++ )
		vector[i] /= mag;
}

//this is the iterative method that fixes errors resulting from floating-point crap.
//sometimes this is called the power method for finding eigenvectors. Basically,
//multiply the matrix by the eVector, normalize, repeat. This gives the eVector
//corresponding to the greatest eValue.
double reOrient(double eVector[4], double N[4][4])
{
	double eigenError = 1, tempVector[4]; //eigenError is a measurement of error used.
	double eValue;
	int count = 0;
	while (eigenError > maxEigenError) //exit once the error is small enough
	{
		normalize4(eVector);
		for (int i = 0; i < 4; i++)
			tempVector[i] = eVector[i]; //load the temp vector
		for (int i = 0; i < 4; i++)
			eVector[i] = N[i][0]*tempVector[0] + N[i][1]*tempVector[1] + N[i][2]*tempVector[2] + N[i][3]*tempVector[3]; //multiply
		eValue = magnitude4(eVector);
		//printf("New Eigen Value: %.17f\n", eValue );
		normalize4(eVector);
		eigenError = 1 - (eVector[0]*tempVector[0] + eVector[1]*tempVector[1] + eVector[2]*tempVector[2] + eVector[3]*tempVector[3]); //compute error
		count += 1;
		if (count > 1000)
		{
			debug_printf("Over 1000 iterations... Not converging fast enough D=\n");
			break;
		}
	}
	if (count <= 9000)
		debug_printf("stopped after %d iterations beyond the linear solution\n", count);
	normalize4(eVector);

	return eValue;
}

void rotateImplicit(double v[3], double u[3])												// rotates v about u by theta = asin(||u||)
{
	double w[3];
	double mag = normalize(u, w);
	double theta = asin(mag);

	rotate(v, w, theta, v);
}

/*weights the stars by their given errors to obtain a better prediction of the rotational matrix once an initial prediction has been made:
The cross product of the actual coordinates and the expected coordinates of each star is taken.  Individually,
for each particular star, this information is sufficient to determine a rotation on the expected coordinates
that will make them the same as the actual coordinates.  However, if we wish to use this information for multiple
stars, we must compromise on which rotation to use.  So we take a weighted average of the cross products, and 
rotate *all* of the stars' expected coordinates by the resulting axis (the angle by which we rotate depends
on the magnitude of this vector).*/
void WeighByError( vector<image_star>& ImageStars, vector<catalog_star>& Catalog, double Rt[3][3], double rOut[3] )
{
	for (int i = 0; i < 3; i++)
		rOut[i] = Rt[i][0];

	//r is the image coordinate, transformed to sky by R
	//find r = Rt (rP). Now r should be close the catalog coordiantes.  In fact, the difference
	//between them is the current error determined by the quaternion algorithm.
	vector<image_star>::iterator it;
	for ( it = ImageStars.begin(); it != ImageStars.end(); it++ )
		if ( it->identity != INDEX_INVALID )
			matrix_mult( Rt, it->r_prime, it->r );

	//print stuff out
	/*for (int i = 0; i < SIZE; i++)
	{
		printf("r:			%f %f %f\n: ", stars[i].r[0], stars[i].r[1], stars[i].r[2]);
		printf("rP_transformed:	%f %f %f\n\n", stars[i].rP_transformed[0], stars[i].rP_transformed[1], stars[i].rP_transformed[2]);
	}*/

	double total_weight = 0;													//used to scale the errors

	for ( it = ImageStars.begin(); it != ImageStars.end(); it++ )
		if ( it->identity != INDEX_INVALID )
			total_weight += 1.0 / SQUARE(it->error);							//compute the total weight to be used

	double error = 0.0,															//used to determine when to escape the loop
	old_error = -1.0,
	cross[3],																	//used in the calculation
	cross_avg[3],																//weighted average of the crosses
	delta_rP_trans[3];															//the change in rP_transform for each star
																				//as determined by cross products
																				//(first order approx.)

	int count = 0;																//break out of loop if gone for over 9000

	while ( abs(error - old_error) > DELTA_ERROR)
	{
		cross_avg[0] = 0;														//initialize these every loop
		cross_avg[1] = 0;
		cross_avg[2] = 0;

		for ( it = ImageStars.begin(); it != ImageStars.end(); it++ )
		{
			if ( it->identity != INDEX_INVALID )
			{
				cross_product( it->r, Catalog[it->identity].r, cross );			//find the cross product of r and the catalog coordinates
				for ( int j = 0; j < 3; j++)
				{
					cross_avg[j] += cross[j] / SQUARE(it->error);				//determined the weighted average of crosses
				}
			}
		}
		
		for (int i = 0; i <3; i++)
			cross_avg[i] /= total_weight;										//take into account the scaling

		old_error = error;														//remember the old error
		error = 0;

		for ( it = ImageStars.begin(); it != ImageStars.end(); it++ )
		{
			if ( it->identity != INDEX_INVALID )
			{
				rotateImplicit( it->r, cross_avg );								//most important part -- rotate the rP_trans

				cross_product( cross_avg, it->r, delta_rP_trans );				//used to determine error

				for (int j = 0; j < 3; j++)
					error += delta_rP_trans[j] * delta_rP_trans[j] * 10000;		//determine the new error from first order apprx
			}
		}

		rotateImplicit( rOut, cross_avg );										//rotate the center of the image by the appropriate amount.

		//printf("old error: %.10f\n",old_error);
		//printf("error = %.10f\n", error);

		count++;																//update counter

		if (count > 9000)
		{
			debug_printf("iterated over 9000 times!!!");						//break out of taking too long
			break;
		}
	}

	debug_printf( "Error Weighing Iterations: %d\n", count );

#ifdef DEBUG_TEXT
	for ( it = ImageStars.begin(); it != ImageStars.end(); it++ )
		if ( it->identity != INDEX_INVALID )
			debug_printf( "star %d | r = (%.10f, %.10f, %.10f)\n", int(it-ImageStars.begin()), it->r[0], it->r[1], it->r[2] );
#endif
}



