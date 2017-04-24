/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
*/

#include "stdafx.h"

#include <cstdio>
#include <cstdarg>
#include <direct.h>
#include <cmath>
#include "ini_reader.h"
#include "detector.h"
#include "catalog.h"
#include "identifier.h"
#include "attitude.h"
#include "timer.h"
#include "atmosphere.h"
#include "star.h"
#include "fitsio.h"

#include <fstream>
#include <string>
#include <vector>

using namespace std;

#define CCD_FRAME_RAW			0
#define CCD_FRAME_BMP			1
#define CCD_MODE_NORMAL			0
#define CCD_MODE_TRIGGER		1

#define FRAME_TYPE				CCD_FRAME_RAW
#define CAMERA_MODE(trigger)	( (trigger) ? CCD_MODE_TRIGGER : CCD_MODE_NORMAL )

#define CAMERA_WIDTH	1280
#define CAMERA_HEIGHT	960

#define FONT_SIZE 16

	//#define PROCESS_IMAGE

const char* settings_file = "settings.txt";
bool fault = false;

	// globals for reading from a directory
string line;
vector<string> names;
const char* namesFile = "names.txt";
char namesFileDir[100];	// set at runtime
int namesIterator = 0;

	// Static framebuffer for storing camera image
void* pixel_data = NULL;

	// Global catalog structure, so it can be allocated in main and used in callback
Catalog glCatalog;

	// Global vector for image stars, to be used repeatedly (avoid constant memory allocation)
vector<image_star> ImageStars;

	// Global variables, modified in callback and viewed in main
	// Transpose (inverse) of rotation matrix
double Image_Rt[3][3] = {	{ 1, 0, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 } };
	// RA,DEC of image
coordinates Image_Coords = {0, 0};
	//Azi, Ele of image
double Image_Azi = 0.0;
double Image_Ele = 0.0;
double Image_Bore = 0.0;
double Image_LST = 0.0;
double Image_HA = 0.0;

// used for focal length adjustment in CaliMode
double running_foc = 0;
int num_im_cali = 0;

int Ele_Corr = 0;
int num_detected = 0;
	// Local time of image
SYSTEMTIME Image_LocalTime;
	// Whether attitude was correctly determined
bool Image_Correct_ID = false;
	// Mean Sky Level
float Mean_Sky_Level = 0.0f;
	// Whether there is a valid prior for use when identifying stars in new image
bool Prior_Valid = false;
	// Contains coords of a previous image to make it easier to identify another image
prior_s Prior;

char dir_name[512];
char fits_dir_name[512];
char output_name[512];
char output_foc[512];

	// Global structure for the dark image, which is an image taken by the camera without the lens on to be subtracted from actual images
image_s DarkImage = {
	0,
	0,
	0,
	NULL
};

	// Settings loaded from file about how the program should run
struct settings_s {
	bool TriggerMode;
	bool DumpRaws;
	bool DumpData;
	bool ProcessImages;
	float FocalLength;
	int CaptureDelay;	// milliseconds
	int BitMode;		// 8 or 12
	int BinMode;		// 1, 2, 3, 4
	int ImageWidth;
	int ImageHeight;
	int ExposureTime;	// milliseconds (1 - 200,000)
	int RedGain;		//
	int GreenGain;		// dB (6 - 41)
	int BlueGain;		//
	float DarkCoeff;
	float Latitude;		// Degrees
	float Longitude;	// 
	float Altitude;		// meters
	bool ReadRaw;
	bool ReadFits;
	bool DumpFits;
	bool CaliMode;		// Auto-calibration. adjusts focal length.
	int CaliTries;
	float FocalError;
	float PriorRA;		// used as prior when CaliMode is on
	float PriorDEC;
	float PriorError;
	float OutRad;
	float InRad;
	bool LostInSpace;
} Settings;

	// CCD Camera settings
struct camera_s {
	int DeviceID;
	char ModuleNo[32];
	char SerialNo[32];
	unsigned long long Timer;
	unsigned long long PrevTimeStamp;
	bool Initialized;
	bool EngineStarted;
} Camera;

#ifdef SDL_ENABLED
	// Easy function to draw text using SDL
void DrawSDLText(TTF_Font* font, SDL_Surface* dst, int x, int y, char* format, ... )
{
	char buffer[256];
	va_list args;
	va_start(args, format);
	vsnprintf(buffer, 256, format, args);
	va_end(args);
	
	SDL_Color textColor = { 255, 255, 255 };
	SDL_Surface* msg = TTF_RenderText_Blended( font, buffer, textColor );
	
	SDL_Rect offset;
	offset.x = (short)x;
	offset.y = (short)y;
	SDL_BlitSurface( msg, NULL, dst, &offset );
	
	SDL_FreeSurface(msg);
}
#endif

	// clamps x < a to x = a, x > b to x = b
int clamp(int x, int a, int b)
{
	if (b < a)
	{
		int tmp = a;
		a = b;
		b = tmp;
	}
	if (x < a)
		return a;
	if (x > b)
		return b;
	return x;
}

	// Function for loading program settings
int CameraLoadSettings( const char* filename )
{
		
	printf( "Welcome to the latest version of the Star Tracker, 5/26/13.\n" );
	printf( "This program is the result of combined efforts by Gil Tabak, Alexander Wickes, and Karanbir Toor, with the instruction of Professor Philip Lubin.\n" );
	printf( "This program has been tested on a Mightex CGE-B013-U Camera.\n" );
	printf( "Please see the README file for further instructions and a description of features offered.\n" );
	printf( "Please also note the Settings.txt file, which the program uses as input.\n" );

	FILE* file = fopen( filename, "r" );
	if ( file == NULL )
		printf( "Error opening settings file.\n" );
	
	Settings.TriggerMode =	IniGetBool( file, "TriggerMode", false );

	Settings.DumpRaws =	IniGetBool( file, "DumpRaws", false );
	Settings.DumpFits =		IniGetBool( file, "DumpFits", false );
	Settings.DumpData =		IniGetBool( file, "DumpData", false );

	Settings.ProcessImages =IniGetBool( file, "ProcessImages", false );
	Settings.BitMode =		( IniGetBool( file, "BitMode12", true ) ? 12 : 8 );
	Settings.FocalLength =	IniGetFloat( file, "FocalLength", 4293.0f );
	Settings.CaptureDelay =	IniGetInt( file, "CaptureDelay", 5000 );
	Settings.BinMode =		clamp( IniGetInt( file, "BinMode", 1 ), 1, 4 );
#if 0
	Settings.ImageWidth =	IniGetInt( file, "ImageWidth", 1280 );
	Settings.ImageHeight =	IniGetInt( file, "ImageHeight", 960 );
#else
	Settings.ImageWidth =	CAMERA_WIDTH/Settings.BinMode - (CAMERA_WIDTH % Settings.BinMode);
	Settings.ImageHeight =	CAMERA_HEIGHT/Settings.BinMode - (CAMERA_HEIGHT % Settings.BinMode);
#endif
	Settings.ExposureTime =	IniGetInt( file, "ExposureTime", 20 );
	Settings.RedGain =		clamp( IniGetInt( file, "RedGain", 14 ), 6, 41 );
	Settings.GreenGain =	clamp( IniGetInt( file, "GreenGain", 14 ), 6, 41 );
	Settings.BlueGain =		clamp( IniGetInt( file, "BlueGain", 14 ), 6, 41 );
	Settings.DarkCoeff =	IniGetFloat( file, "DarkCoeff", 1.0f );
	Settings.Latitude =		IniGetFloat( file, "Latitude", 0.0f );
	Settings.Longitude =	IniGetFloat( file, "Longitude", 0.0f );
	Settings.Altitude =		IniGetFloat( file, "Altitude", 0.0f );
	Settings.ReadRaw =	IniGetBool( file, "ReadRaw", false );
	Settings.ReadFits =	IniGetBool( file, "ReadFits", false );
	Settings.CaliMode =		IniGetBool( file, "CaliMode", false);
	Settings.CaliTries =	IniGetInt( file, "CaliTries", 10 );
	Settings.FocalError =	IniGetFloat( file, "FocalError", 0.01);
	Settings.PriorRA =		IniGetFloat( file, "PriorRA", 999);
	Settings.PriorDEC =		IniGetFloat( file, "PriorDEC", 999);
	Settings.PriorError =	IniGetFloat( file, "PriorError", 0.4);
	Settings.OutRad =		IniGetFloat( file, "OutRad", 3.0);
	Settings.InRad =		IniGetFloat( file, "InRad", 1.3);
	Settings.LostInSpace =	IniGetFloat( file, "LostInSpace", FALSE);
	
	if ( file != NULL )
		fclose(file);
	
	printf( "TriggerMode: %s\n", (Settings.TriggerMode ? "True" : "False") );
	printf( "DumpRaws: %s\n", (Settings.DumpRaws ? "True" : "False") );
	printf( "DumpData: %s\n", (Settings.DumpData ? "True" : "False") );
	printf( "ProcessImages: %s\n", (Settings.ProcessImages ? "True" : "False") );
	printf( "BitMode: %d\n", Settings.BitMode );
	printf( "FocalLength: %f\n", Settings.FocalLength );
	printf( "CaptureDelay: %d\n", Settings.CaptureDelay );
	printf( "BinMode: %d\n", Settings.BinMode );
	printf( "ImageWidth: %d\n", Settings.ImageWidth );
	printf( "ImageHeight: %d\n", Settings.ImageHeight );
	printf( "ExposureTime: %d\n", Settings.ExposureTime );
	printf( "RedGain: %d\n", Settings.RedGain );
	printf( "GreenGain: %d\n", Settings.GreenGain );
	printf( "BlueGain: %d\n", Settings.BlueGain );
	printf( "DarkCoeff: %f\n", Settings.DarkCoeff );
	printf( "Latitude: %f\n", Settings.Latitude );
	printf( "Longitude: %f\n", Settings.Longitude );
	printf( "Altitude: %f\n", Settings.Altitude );
	//printf( "DumpFits: %s\n", Settings.DumpFits );
	printf( "ReadRaw: %s\n", (Settings.ReadRaw ? "True" : "False") );
	printf( "ReadFits: %s\n", (Settings.ReadFits ? "True" : "False") );
	// The following print statements made the code crash. Not sure why.
	/*printf( "CaliMode: %s\n", Settings.CaliMode );
	printf( "CaliTries: %d\n", Settings.CaliTries );
	printf( "FocalError: %f\n", Settings.FocalError );
	printf( "PriorRA: %f\n", Settings.PriorRA );
	printf( "PriorDEC: %f\n", Settings.PriorDEC );*/
	
	return 1;
}

#ifdef SDL_ENABLED

SDL_Surface* CCD_surface = NULL;
SDL_Surface* screen = NULL;
TTF_Font* font = NULL;
int yText = 0;

void GUISetCursor( int y )
{
	yText = y;
}

void GUIFlipScreen()
{
	SDL_BlitSurface( CCD_surface, NULL, screen, NULL );
	GUISetCursor(0);
}

void GUIPrintText( const char* txt, ... )
{
	char buffer[4*1024];
	va_list argptr;
	va_start(argptr, txt);
	vsnprintf(buffer, 4*1024, txt, argptr);
	
	DrawSDLText( font, screen, 0, yText, buffer );
	
	va_end(argptr);
	
	yText += FONT_SIZE;
}

#endif

void printTo_FIT_file(FILE * pFile, unsigned char* data){
	string temp = "";
    // Header is 469 bytes long.
	temp += "SIMPLE  =     T";
	for (int i = 0; i<80-15; i++)
		temp += " ";
	temp += "BITPIX  =    16";
	for (int i = 0; i<80-15; i++)
		temp += " ";
	temp += "NAXIS   =     2";
	for (int i = 0; i<80-15; i++)
		temp += " ";
	temp += "NAXIS1  =   1280";
	for (int i = 0; i<80-16; i++)
		temp += " ";
	temp += "NAXIS2  =   960";
	for (int i = 0; i<80-15; i++)
		temp += " ";
	temp += "END";
	for (int i = 0; i<80-15; i++)
		temp += " ";
	fwrite(temp.c_str(), temp.length(), 1, pFile);
	
    // Now I add the padding because the header must be 2880 bytes long. 
	unsigned char c = 0;
	unsigned char* buffer = new unsigned char[2880-468];
	for (int i = 0; i<2880-468; i++) {
		buffer[i] = c;
	}
	fwrite(buffer, 2880-468, 1, pFile);
		//cout << "wrote header size: " << temp.length() + 2880-468 << endl;
	
	fwrite (data , 2*CAMERA_WIDTH*Settings.ImageHeight , 1 , pFile );
	
	unsigned char* buffer2 = new unsigned char[2880 - ((2*CAMERA_WIDTH*Settings.ImageHeight) % 2880)];
	c = 0;
	for (int i=0; i< 2880 - ((2*CAMERA_WIDTH*Settings.ImageHeight) % 2880) ; i++) {
		buffer2[i] = c;
	}
	
	fwrite(buffer2, 2880 - ((2*CAMERA_WIDTH*Settings.ImageHeight) % 2880), 1, pFile);
	
	delete buffer2;
	delete buffer;
	fclose (pFile);
}

unsigned char* readFitsFileIntoBuffer(const char* f) {
	FILE * pFile;
	long lSize;
	unsigned char * buffer;
	size_t result;
	
	pFile = fopen ( f , "rb" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	
		// obtain file size:
	rewind (pFile);
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	fseek ( pFile , 2880 , SEEK_SET );
	
		// allocate memory to contain the whole file:
	buffer = new unsigned char[CAMERA_WIDTH*Settings.ImageHeight*2];
		//if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
	
		// copy the file into the buffer:
	result = fread (buffer,CAMERA_WIDTH*Settings.ImageHeight*2, 1,pFile);
		//if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
	
	/* the whole file is now loaded in the memory buffer. */
	
		// terminate
	fclose (pFile);
	
	return buffer;
}

unsigned char* extractDataFromRawFile(std::string fileName){
  FILE * pFile;
  long lSize;
  unsigned char * buffer;
  size_t result;
	
  pFile = fopen ( "./data/image_20130518_190308696.raw" , "rb" );
  if (pFile==NULL) {fputs ("File error\n",stderr);}
  
		// obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);
	
		// allocate memory to contain the whole file:
  buffer = (unsigned char*) malloc (sizeof(unsigned char)*lSize);
  if (buffer == NULL)  {printf("Memory error\n");}
	
		// copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize)  {printf("Reading error\n");}
	
		// terminate
  fclose (pFile);
  printf("---------------------------\n");
  printf("%d\n",lSize);
  printf("---------------------------\n");
  return (unsigned char*)buffer;
}

void putFileNamesIntoVector(const char* dir){
	FILE * pFile;
  long lSize;
  char * buffer;
  size_t result;
	
	char nameFilePath[100];
	sprintf(nameFilePath, "./%s/names.txt",dir);
  pFile = fopen ( nameFilePath , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr);}
	
		// obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);
	
		// allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {printf("Memory error\n");}
	
		// copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {printf("Reading error\n");}
	
  /* the whole file is now loaded in the memory buffer. */
	
	// push individual file names into the names vector
  int i=0;
  int j;
  while(i<lSize){
	  for(j=i; buffer[j] != '\n'; j++)
		  j;
	  names.push_back(string(&buffer[i],j-i));
	  i=j+1;
  }

  std::vector<string>::size_type sz = names.size();
  printf("read the following names from names.txt\n");
  for(unsigned int i=0; i<sz; i++)
	  printf("%s\n",names[i]);
		// terminate
  fclose (pFile);
  free (buffer);
}

int to_int(char c)
{
     switch(c)
     {   //  makes no assumptions about character order
          case '0': return 0 ;
          case '1': return 1 ;
          case '2': return 2 ;
          case '3': return 3 ;
          case '4': return 4 ;
          case '5': return 5 ;
          case '6': return 6 ;
          case '7': return 7 ;
          case '8': return 8 ;
          case '9': return 9 ;
          default: return -1 ;
     }
}

//gets date from image name, and puts it in Image_LocalTime struct.
void parseFileNameForDate(char* filename){

	filename += 20;
	for (; *filename != '_'; filename++); //loop until '_' after image
	filename++; // ptr == '_', skip

	Image_LocalTime.wYear = to_int(*filename++)*1000;
	Image_LocalTime.wYear += to_int(*filename++)*100;	
	Image_LocalTime.wYear += to_int(*filename++)*10;
	Image_LocalTime.wYear += to_int(*filename++);
	
	Image_LocalTime.wMonth = to_int(*filename++)*10;
	Image_LocalTime.wMonth += to_int(*filename++);
	
	Image_LocalTime.wDay = to_int(*filename++)*10 ;
	Image_LocalTime.wDay += to_int(*filename++);

	*filename++; // ptr == '_', skip
	
	Image_LocalTime.wHour = to_int(*filename++)*10;
	Image_LocalTime.wHour += to_int(*filename++);
	
	Image_LocalTime.wMinute = to_int(*filename++)*10;
	Image_LocalTime.wMinute += to_int(*filename++);
	
	Image_LocalTime.wSecond = to_int(*filename++)*10 ;
	Image_LocalTime.wSecond += to_int(*filename++);

	Image_LocalTime.wMilliseconds = to_int(*filename++)*100;
	Image_LocalTime.wMilliseconds += to_int(*filename++)*10;
	Image_LocalTime.wMilliseconds += to_int(*filename++);

//	printf("Image Time IN PARSE METHOD is %04d%02d%02d_%02d%02d%02d%03d\n", Image_LocalTime.wYear, Image_LocalTime.wMonth, Image_LocalTime.wDay,
//	Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds );

	return;
}

////////////////////////////////////////
//	Methods Used In Calibration Mode  //
////////////////////////////////////////

void updatePriorFocalLength()
{
	Prior.focal_length = Settings.FocalLength;
}

/*
if calibration mode is on, try different focal lengths until one works.
Does not save the focal length that worked, only whether there was a correct ID. 
Specify bounds and increments in settings.txt file.
Use focal_length as the starting point, and move outwards.
*/
void findValidSolution(	vector<image_star>& ImageStars, vector<catalog_star>& CatalogStars, 
					coordinates* BoresightCoords, coordinates* Coords, sfloat R[3][3], double RQuat[4], double &err)
{
	float currentFocal = Settings.FocalLength;
	float startingFocal = Settings.FocalLength;//Prior_Valid? Prior.focal_length: Settings.FocalLength;

	//max error, in pixels
	float maxError = Settings.FocalError * Settings.FocalLength;
	// maybe try:
	Image_Correct_ID = false;
	for(  int f = 0; f < Settings.CaliTries && !Image_Correct_ID; f++ ){

		//add or subtract, moving outwards from the best guess.
		if (f%2 == 0)
			currentFocal = startingFocal + ( (float) f) * maxError / Settings.CaliTries;
		else currentFocal = startingFocal - ( (float) f) * maxError / Settings.CaliTries;

		//adjust image stars r_prime vector using GetSphericalFromImage with different focal lengths:
		for (int i = 0; i < ImageStars.size(); i++)
			GetSphericalFromImage(	ImageStars[i].centroid_x - IMAGE_WIDTH/2.0f + 0.5f, 
			ImageStars[i].centroid_y - IMAGE_HEIGHT/2.0f + 0.5f,
			currentFocal, ImageStars[i].r_prime );

		// enough stars, continue with identification
		IdentifyImageStars(ImageStars, glCatalog, ( Prior_Valid ? &Prior : NULL ));

		double rms_error = GetAttitude(ImageStars, CatalogStars, BoresightCoords, Coords, R, RQuat, err);
		Image_Correct_ID = ( rms_error <= 0.01 );
	}
	printf("focal length that gives solutions is %f\n", currentFocal);
}

/*
Finds the best focal length by minimizing the sum of errors
of the identified image stars as a function of focal length

If both caliMode and dumpData are on, then record total error
versus focal length in foc_... file
*/

void findBestFocalLength(	vector<image_star>& ImageStars, vector<catalog_star>& CatalogStars, 
					coordinates* BoresightCoords, coordinates* Coords, sfloat R[3][3], double RQuat[4], double &err)
{
	FILE* file;
	if (Settings.CaliMode && Settings.DumpData)
		file = fopen(output_foc, "a");
	else file = NULL;

	double currentFocal = Settings.FocalLength - Settings.FocalError * Settings.FocalLength;
	double min_err = err; // we'll find the min err. Set to the one retrieved previously, which is >= min_err.
	double min_focal = currentFocal;

	vector<double> rms_error_trials;
	rms_error_trials.clear();
	rms_error_trials.reserve(Settings.CaliTries);
	for (int f = 0; f < Settings.CaliTries; f++){

		//adjust image stars r_prime vector using GetSphericalFromImage with different focal lengths:
		for (int i = 0; i < ImageStars.size(); i++)
			GetSphericalFromImage(	ImageStars[i].centroid_x - IMAGE_WIDTH/2.0f + 0.5f, 
			ImageStars[i].centroid_y - IMAGE_HEIGHT/2.0f + 0.5f,
			currentFocal, ImageStars[i].r_prime );

		currentFocal += 2*Settings.FocalError * Settings.FocalLength / Settings.CaliTries;
		GetAttitude(ImageStars, glCatalog.Stars, BoresightCoords, Coords, R, RQuat, err) ;
		rms_error_trials.push_back(err);

		if ( file != NULL ){
			fprintf(file, "%.3f ", currentFocal);
		}

		if (err < min_err){
			min_err = err;
			min_focal = currentFocal;
		}
	}
	if ( file != NULL ){
		fprintf(file, "\n");
	}

	//print min focal length and err.
	//printf("min foc: %f | min err: %f\n",min_focal, min_err);

	currentFocal = Settings.FocalLength - Settings.FocalError * Settings.FocalLength;
	//print all focals and their error
	// for (int f = 0; f < Settings.CaliTries; f++){
	// currentFocal += Settings.FocalError * Settings.FocalLength / Settings.CaliTries;
	// printf("Foc Len: %f | err: %f\n",currentFocal, rms_error_trials[f]);
	// }

	if ( file != NULL )
		for (int i = 0; i < rms_error_trials.size(); i++){
			fprintf(file, "%.7f ", rms_error_trials[i]);
		}
	if ( file != NULL ){
		fprintf(file, "\n");
		}

		num_im_cali ++;
		running_foc = (running_foc*(num_im_cali-1) + min_focal) / num_im_cali;
		if (num_im_cali%20 == 9)
			Settings.FocalLength = running_foc;

		//print calibrated focal length and err.
		printf("running calibrated focal length: %.2f\n",running_foc);
		if (file != NULL)
			fclose(file);
} 

//////////////////////////////////////
//  End Calibration Mode Methods  ////
//////////////////////////////////////

//////////////////////////////////////
//  Dump Images Methods  /////////////
//////////////////////////////////////

void saveRawImages(unsigned char* BytePtr)
{
	// Save image to file w/ timestamp
	char buffer[512];

	// get system timestamp (millisecond accuracy)
	int nameSize = sprintf( buffer, "./%s/image_%04d%02d%02d_%02d%02d%02d%03d.raw", dir_name, Image_LocalTime.wYear, Image_LocalTime.wMonth, Image_LocalTime.wDay,
		Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds );

	// create a file containing the names of all the files in the directoy
	char nameFileBuffer[512];
	int lastIndex = sprintf(nameFileBuffer, "./%s/names.txt",dir_name);
	FILE* namesFile = fopen(nameFileBuffer,"a");
	fwrite(buffer, nameSize, 1, namesFile );
	fwrite("\n", 1, 1, namesFile);
	fclose(namesFile);
	// -------------------------------------------------------------------


	FILE* file = fopen(buffer, "wb");
	if ( file == NULL )
	{
		printf( "Error: Could not open file.\n" );
	}
	else
	{
		int pixel_bytes = ( Settings.BitMode == 8 ? 1 : 2 );
		fwrite( BytePtr, CAMERA_WIDTH*Settings.ImageHeight*pixel_bytes, 1, file );
		fclose(file);
		//printf( "Wrote to file successfully.\n" );
	}
}

void saveFITSImages(unsigned char* BytePtr)
{
		unsigned short* BytePtr16 = (unsigned short*)BytePtr;
		fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
		int status, ii, jj;
		long  fpixel, nelements, exposure;

		/* initialize FITS image parameters */
		int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
		long naxis    =   2;  /* 2-dimensional image                            */    
		long naxes[2] = { CAMERA_WIDTH , Settings.ImageHeight }; 

		status = 0;         /* initialize status before calling fitsio routines */

		// Save image to file w/ timestamp
		char buffer[512];

		// get system timestamp (millisecond accuracy)
		int nameSize = sprintf( buffer, "./%s/image_%04d%02d%02d_%02d%02d%02d%03d.fits", fits_dir_name, Image_LocalTime.wYear, Image_LocalTime.wMonth, Image_LocalTime.wDay,
			Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds );

		if (fits_create_file(&fptr, buffer, &status)) /* create new FITS file */
			printf("Error: Couldn't create a FITS file");           /* call printerror if error occurs */

			/* write the required keywords for the primary array image.     */
			/* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
			/* a FITS image with BITPIX = 16 (signed short integers) with   */
			/* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
			/* FITS uses to store unsigned integers.  Note that the BSCALE  */
			/* and BZERO keywords will be automatically written by cfitsio  */
			/* in this case.                                                */

		if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
			printf( "Error: Couldn't create the image" );

		fpixel = 1;                               /* first pixel to write      */
		nelements = naxes[0] * naxes[1];          /* number of pixels to write */

		//Writes the image to the file
		if ( fits_write_img(fptr, TUSHORT, fpixel, nelements, BytePtr16, &status) )
			printf("Error: Couldn't write the image");

		if ( fits_close_file(fptr, &status) )                /* close the file */
			printf("Error: Couldn't close the file");
		
		// create a file containing the names of all the files in the directoy
		char nameFileBuffer[512];
		int lastIndex = sprintf(nameFileBuffer, "./%s/names.txt",fits_dir_name);
		FILE* namesFile = fopen(nameFileBuffer,"a");
		fwrite(buffer, nameSize, 1, namesFile );
		fwrite("\n", 1, 1, namesFile);
		fclose(namesFile);
		// -------------------------------------------------------------------
}

//////////////////////////////////////
//  Dump Images Methods  /////////////
//////////////////////////////////////

// Callback function for handling 

// broke apart FrameCallBack. This function is called either through the camera FrameCallBack function,
// or directly if Settings.readRaw is true.
void FrameCallBackHelper(unsigned char* BytePtr )
{
#ifdef TIMERS_ON 
	Timer FrameCallBackHelperTimer;
	FrameCallBackHelperTimer.StartTimer();
#endif

	if ( Settings.CaliMode )
		void updatePriorFocalLength();
	else running_foc = Settings.FocalLength;
	if ( Settings.PriorRA != 999 && Settings.PriorDEC !=999 && Prior_Valid == false && !Settings.LostInSpace){
		Prior.coords.RA = Settings.PriorRA;
		Prior.coords.DEC = Settings.PriorDEC;
		Prior_Valid = true;
	}
		SYSTEMTIME localtime;
		GetLocalTime( &localtime );
		Image_LocalTime = localtime;

	// return when passed NULL, unless reading from file (when NULL is used)
	if ( BytePtr == NULL && !Settings.ReadRaw && !Settings.ReadFits)
		return;

	// If in 12 bit mode, fix the messed up bit ordering.
	// No need to run this if reading images from file
	if ( Settings.BitMode == 12 && !Settings.ReadRaw && !Settings.ReadFits )
	{
		int i;
		unsigned short* BytePtr16 = (unsigned short*)BytePtr;
		for ( i = 0; i < CAMERA_WIDTH*Settings.ImageHeight; i++ )
			BytePtr16[i] = (BytePtr[2*i] << 4) | BytePtr[2*i+1];
	}

	if(Settings.ReadRaw){

#ifdef TIMERS_ON
		Timer ReadRawTimer;
		ReadRawTimer.StartTimer();
#endif

		if (namesIterator < names.size()){
			char fileName[512];
			//char hardCode[] = "./data_20130518_193209/image_20130518_193209735.raw";
			//sprintf(fileName, "%s", names[namesIterator].data());
			printf("Name     : %s\n",names[namesIterator]);
			//printf("char name: %s\n",fileName);
			BytePtr = new unsigned char[CAMERA_WIDTH*Settings.ImageHeight*2];

			int pos = 0;
			for(int i=0; i<names[namesIterator].length()-1; i++){
				fileName[i] = names[namesIterator].data()[i];
			}

			fileName[names[namesIterator].length()-1] = 0;
			
			//if reading from image, read the time off the image. Overwrites Image_LocalTime.
			 parseFileNameForDate(fileName);
			//printf("Image Time is %04d%02d%02d_%02d%02d%02d%03d\n", Image_LocalTime.wYear, Image_LocalTime.wMonth, Image_LocalTime.wDay,
		//	Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds );


			//FILE* inputFile = fopen("./data_20130518_193209/image_20130518_193209735.raw", "r");
			FILE* inputFile = fopen(fileName, "rb");
			if ( inputFile == NULL ){
				printf( "Error: Could not open file.\n" );
			} else {
				printf("Raw file Open.\n");
				fread (BytePtr,CAMERA_WIDTH*Settings.ImageHeight*2,1,inputFile);
				fclose(inputFile);
			}
			//initialize the catalog.
			if (namesIterator == 0)
				glCatalog.Initialize( "catalog.txt", &Image_LocalTime, Settings.FocalLength );

			namesIterator++; //next time, access the next image
			//return;			
#ifdef TIMERS_ON
	ReadRawTimer.StopTimer();
	printf("Successfully Read File. | Time elapsed: %d ms \n", ReadRawTimer.GetTime());
#endif
		} else
			return;
	}

	if(Settings.ReadFits){

#ifdef TIMERS_ON
		Timer ReadFitsTimer;
		ReadFitsTimer.StartTimer();
#endif

		if (namesIterator < names.size()){
			char fileName[512];
			//char hardCode[] = "./data_20130518_193209/image_20130518_193209735.raw";
			//sprintf(fileName, "%s", names[namesIterator].data());
			printf("Name     : %s\n",names[namesIterator]);
			//printf("char name: %s\n",fileName);
			BytePtr = new unsigned char[CAMERA_WIDTH*Settings.ImageHeight*2];
			unsigned short* BytePtr16 = (unsigned short*)BytePtr;

			for(int i=0; i<names[namesIterator].length()-1; i++){
				fileName[i] = names[namesIterator].data()[i];
			}
			fileName[names[namesIterator].length()-1] = 0;

			
			//if reading from image, read the time off the image. Overwrites Image_LocalTime.
			parseFileNameForDate(fileName);

			////////////////////////
			//FILE* inputFile = fopen("./data_20130518_193209/image_20130518_193209735.raw", "r");
			fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
			int status, ii, jj;
			long  fpixel, nelements, exposure;

			/* initialize FITS image parameters */
			int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
			long naxis    =   2;  /* 2-dimensional image                            */    
			long naxes[2] = { CAMERA_WIDTH , Settings.ImageHeight }; 

			status = 0;         /* initialize status before calling fitsio routines */
			printf(">>>>>>>>>>>>>>>>>>>>\n");
			if (fits_open_file(&fptr, fileName, READONLY, &status)) /* open FITS file */
				printf("Error: Couldn't open a FITS file\n");           /* call printerror if error occurs */

			int anul;
			if (fits_read_img(fptr, USHORT_IMG, 1, (long long)CAMERA_WIDTH*Settings.ImageHeight, NULL, (void*)BytePtr16, &anul, &status)) /* read data from FITS file */
				//printf("Error: Couldn't read a FITS file\n");           /* call printerror if error occurs */

			if ( fits_close_file(fptr, &status) )                /* close the file */
				printf("Error: Couldn't close the file\n");
			/*FILE* inputFile = fopen(fileName, "rb");
			if ( inputFile == NULL ){
				printf( "Error: Could not open file.\n" );
			} else {
				printf("Raw file Open.\n");
				// 2880 is the size of the fits header
				fseek(inputFile, 2880, SEEK_SET);
				fread (BytePtr,1,CAMERA_WIDTH*Settings.ImageHeight*2,inputFile);
				fclose(inputFile);
			}*/
			//////////////
			printf("<<<<<<<<<<<<<<<<<<<<<<\n");

			//initialize the catalog.
			if (namesIterator == 0)
				glCatalog.Initialize( "catalog.txt", &Image_LocalTime, Settings.FocalLength );

			namesIterator++; //next time, access the next image
			//return;			
#ifdef TIMERS_ON
	ReadFitsTimer.StopTimer();
	printf("Successfully Read File. | Time elapsed: %d ms \n", ReadFitsTimer.GetTime());
#endif
		} else
			return;
	}

	if (Settings.DumpRaws)
		saveRawImages(BytePtr);

	if (Settings.DumpFits)
		saveFITSImages(BytePtr);

#ifdef TIMERS_ON
	Timer SumTimer;
	SumTimer.StartTimer();
#endif

	int bytes_pp = (Settings.BitMode == 12 ? 2 : 1);
	unsigned short* BytePtr16 = (unsigned short*)BytePtr;
	// If Binning is enabled, average pixels horizontally (camera only bins vertically)
	if ( Settings.BinMode > 1 )
	{
		long long sum;
		int x, y, i;
		for ( y = 0; y < Settings.ImageHeight; y++ )
		{
			for ( x = 0; x < Settings.ImageWidth; x++ )
			{
				sum = 0;
				for ( i = 0; i < Settings.BinMode; i++ )
				{
					if ( Settings.BitMode == 12 )
						sum += BytePtr16[y*CAMERA_WIDTH + x*Settings.BinMode+i];
					else
						sum += BytePtr[y*CAMERA_WIDTH + x*Settings.BinMode+i];
				}
				if ( Settings.BitMode == 12 )
					BytePtr16[y*Settings.ImageWidth + x] = sum / Settings.BinMode;
				else
					BytePtr[y*Settings.ImageWidth + x] = sum / Settings.BinMode;
			}
		}
	}
#ifdef TIMERS_ON
	SumTimer.StopTimer();
	printf("Successfully computed sum| Time elapsed: %d ms \n", SumTimer.GetTime());
#endif

	
#ifdef TIMERS_ON
	Timer ProcessingTimer;
	ProcessingTimer.StartTimer();
#endif
	// Copy temporary framebuffer to static framebuffer for display in GUI
	memcpy( pixel_data, BytePtr, Settings.ImageWidth*Settings.ImageHeight*bytes_pp );

	if ( Settings.ProcessImages )
	{

#ifdef TIMERS_ON
	Timer DarkImageTimer;
	DarkImageTimer.StartTimer();
#endif
		if ( DarkImage.data != NULL )
		{
			// Subtract the dark image.
			int i;
			short* BytePtr16s = (short*)BytePtr;
			unsigned short* BytePtr16 = (unsigned short*)BytePtr;
			unsigned short* DarkPtr16 = (unsigned short*)DarkImage.data;
			for ( i = 0; i < Settings.ImageWidth*Settings.ImageHeight; i++ )
				BytePtr16s[i] = short( BytePtr16[i] - Settings.DarkCoeff * DarkPtr16[i] );
		}

#ifdef TIMERS_ON
	DarkImageTimer.StopTimer();
	printf("Successfully subtracted dark. | Time elapsed: %d ms \n", DarkImageTimer.GetTime());
#endif

		if ( glCatalog.initialized )
		{
			detector_s Detector = {
				GL_SAMPLE_SKIP,
				LOCAL_WIDTH,
				LOCAL_HEIGHT,
				LOCAL_SAMPLE_SKIP,
				STAR_MIN_OUTSTND,
				0.0
			};
			image_s Image = {
				Settings.ImageWidth,
				Settings.ImageHeight,
				Settings.BitMode,
				Settings.FocalLength,
				BytePtr,
				NULL
			};

			num_detected = DetectStars(ImageStars, &Detector, &Image, Settings.OutRad, Settings.InRad);
			//printf( "	Number of Stars Detected: %d\n", num_detected );
			Mean_Sky_Level = Detector.mean_sky < 100000? (float)Detector.mean_sky: 0;


			if (ImageStars.size() >= 3)
			{

				coordinates Coords, BoresightCoords;
				sfloat R[3][3];
				sfloat RQuat[4];
				double err = 0;

				if (Settings.CaliMode){
					printf("finding valid solution in calibration mode\n" );
					findValidSolution(ImageStars, glCatalog.Stars, &BoresightCoords, &Coords, R, RQuat, err);
				}
				// else calibration mode is NOT on
				else{ 
					printf("finding valid solution -- no calibration mode\n" );
					// enough stars, continue with identification
					IdentifyImageStars(ImageStars, glCatalog, ( Prior_Valid ? &Prior : NULL ));
					double rms_error = GetAttitude(ImageStars, glCatalog.Stars, &BoresightCoords, &Coords, R, RQuat, err);
					Image_Correct_ID = ( rms_error <= 0.01 );
				}

				if ( Image_Correct_ID )
				{
					
					if (Settings.CaliMode)
						findBestFocalLength(ImageStars, glCatalog.Stars, &BoresightCoords, &Coords, R, RQuat, err);

					Prior.coords = Coords;
					Prior.tickCount = GetTickCount();
					if (!Settings.LostInSpace)
						Prior_Valid = true;
					//Prior.focal_length = currentFocal;

					Image_Coords = Coords;

					double JDN;
					// Get Az, El
					if (Settings.ReadRaw || Settings.ReadFits)
						JDN = getJulianDate( &Image_LocalTime );
					else{
						SYSTEMTIME systime;
						GetSystemTime( &systime );
						JDN = getJulianDate( &systime );
					}
					double LAT = Settings.Latitude * M_PI / 180;
					double LONG = Settings.Longitude * M_PI / 180;
					//printf("Julian: %f\n", JDN);
					Image_LST = getLST( LONG, JDN );
					Image_HA = getHA1( Image_LST, Coords.RA );
					Image_Azi = getAzi( Image_HA, LAT, Coords.DEC );
					Image_Ele = getEle( Image_HA, LAT, Coords.DEC );
					//Ele_Corr = correctEle( );

					//find the Boresight Angle from BoresightCoords
					double HA, AZI;
					HA = getHA1( Image_LST, BoresightCoords.RA );
					AZI = getAzi( HA, Settings.Latitude, BoresightCoords.DEC );
					Image_Bore = acos ( cos(Image_Azi) * sin(AZI) - sin(Image_Azi) * cos(AZI) );

					int i, j;
					for ( i = 0; i < 3; i++ )
						for ( j = 0; j < 3; j++ )
							Image_Rt[i][j] = R[j][i];
				}
				// no solution:
				else if ( Prior_Valid && GetTickCount() - Prior.tickCount >= PRIOR_TIMEOUT*1000 )
				{
					Prior_Valid = false;
				}

				if ( Settings.DumpData )
				{
					// Write data to output.txt
					FILE* file = fopen(output_name, "a");
					if ( file != NULL )
					{
						// change to CSV output
						if ( Image_Correct_ID )
						{
							coords_discrete coords;
							GetDiscreteCoords(&Coords, &coords);
							fprintf(file, "%02d/%02d/%04d "
								"%02d:%02d:%02d.%03d,"
								"%.4f, %.4f, %.1f,"
								"%.6f, %.6f,"
								"%02d:%02d:%ffits,"
								"%02d:%02d:%f,"
								"%.4f,%.4f,%d,"
								"%.1f,%d,"
								"%f,%f,%f,%f,"
								"%f\n",
								Image_LocalTime.wMonth, Image_LocalTime.wDay, Image_LocalTime.wYear,
								Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds,
								Settings.Latitude, Settings.Longitude, Settings.Altitude,
								Coords.RA*180/M_PI, Coords.DEC*180/M_PI,
								coords.RA_hr, coords.RA_min, coords.RA_sec,
								coords.DEC_deg, coords.DEC_min, coords.DEC_sec,
								Image_Azi*180/M_PI, Image_Ele*180/M_PI, Ele_Corr,
								Detector.mean_sky, num_detected,
								RQuat[0], RQuat[1], RQuat[2], RQuat[3],running_foc);
							//			Prior.focal_length);
						}
						else
						{
							fprintf(file, "%02d/%02d/%04d "
								"%02d:%02d:%02d.%03d,"
								"%.4f, %.4f, %.1f,"
								"ERR,ERR,ERR,ERR,"
								"ERR,ERR,ERR,"
								"%.1f,%d,"
								"ERR,ERR,ERR,ERR,ERR\n",
								Image_LocalTime.wMonth, Image_LocalTime.wDay, Image_LocalTime.wYear,
								Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds,
								Settings.Latitude, Settings.Longitude, Settings.Altitude,
								Detector.mean_sky, num_detected );
						}
						fclose(file);
					}
				}
			}
			else
			{
				printf("	Not enough stars to continue.\n");
			}
		}
	} // Settings.ProcessImages
#ifndef SDL_ENABLED

#ifdef TIMERS_ON
	Timer PrintingTimer;
	PrintingTimer.StartTimer();
#endif

	if ( Image_Correct_ID )
		printf( "(RA, DEC) = (%.6f, %.4f) | (AZI, ELE) = (%.4f, %.4f) | Boresight: %.4f | Stars: %d | Mean: %.1f\n"
		"EleCorr: %d | LST: %f | HA: %f | (LAT, LONG) = (%.4f,%.4f) | ALT = %.1f | Focal Length =%.2f\n",
		Image_Coords.RA*12/M_PI, Image_Coords.DEC*180/M_PI,
		Image_Azi*180/M_PI, Image_Ele*180/M_PI, Image_Bore*180/M_PI,
		num_detected, Mean_Sky_Level,
		Ele_Corr, Image_LST, Image_HA,
		Settings.Latitude, Settings.Longitude, Settings.Altitude,
		running_foc);
	else
		printf( "(RA, DEC) = ERR | (AZI, ELE) = ERR | Boresight: ERR | Stars: %d | Mean: %.1f\n"
		"EleCorr: ERR | LST: ERR | HA: ERR | (LAT, LONG) = (%.4f,%.4f) | ALT = %.1f | Focal Length =%.2f\n",
		num_detected, Mean_Sky_Level,
		Settings.Latitude, Settings.Longitude, Settings.Altitude, running_foc );
#ifdef TIMERS_ON
	PrintingTimer.StopTimer();
	printf("Successfully Printed values| Time elapsed: %d ms \n", PrintingTimer.GetTime());
#endif

#endif //SDL_ENABLED

	//when reading from file, we filled up BytePtr in FrameCallBackHelper.
	if ( (Settings.ReadRaw || Settings.ReadFits) && BytePtr != NULL)
		delete[] BytePtr;

#ifdef TIMERS_ON
	FrameCallBackHelperTimer.StopTimer();
	printf("Completed entire frame. | Time elapsed: %d ms\n", FrameCallBackHelperTimer.GetTime() );
#endif

}

// Called when the camera is working. Passes to FrameCallBackHelper if device matches and not in read mode.
// When ReadRaw is on, the FrameCallBackHelper is called directly.
void FrameCallBack(TProcessedDataProperty* Attributes, unsigned char* BytePtr ){
	if ( Attributes->CameraID == Camera.DeviceID && !Settings.ReadRaw && !Settings.ReadFits) // The working camera, and NOT reading from files.
		FrameCallBackHelper( BytePtr );
}

	// Called when there is a problem with the camera
void CameraFaultCallBack( int ImageType )
{
	printf( "Error: Camera fault.\n" );
	fault = true;
}

int main(int argc, char* argv[])
{
	MSG msg; // Required for camera to interface with windows
	
	int ret;
	
	ImageStars.reserve( 500 );
	
	printf( "Loading camera settings from file: %s\n", settings_file );
	CameraLoadSettings(settings_file);


	if(Settings.ReadRaw || Settings.ReadFits){
		printf("Input directory you wish to read from.\n");
		scanf("%s",namesFileDir);
		printf("you entered: %s\n",namesFileDir);
		putFileNamesIntoVector(namesFileDir);
		//extractDataFromRawFile("./data_20130518_175503/image_20130518_175504026.raw");
	}
	
		// Setup Priors
	Prior.coords.RA = 0.0f;
	Prior.coords.DEC = 0.0f;
	Prior.tickCount = GetTickCount();
	
	if ( Settings.BitMode != 12 )
	{
		printf( "Error: Program only supports 12 bit images.\n" );
		return 0;
	}
	
	//Prior.phi_max = MAX_PRIOR_DIST;
	Prior.phi_max = Settings.PriorError;
	
		// Load the dark image
	DarkImage.width = Settings.ImageWidth;
	DarkImage.height = Settings.ImageHeight;
	DarkImage.bitmode = Settings.BitMode;
	DarkImage.data = new unsigned short[CAMERA_WIDTH*CAMERA_HEIGHT];
	
	if (DarkImage.data != NULL)
	{
		FILE* file = fopen("dark_image.raw", "rb");
		if (file != NULL)
		{
			fread(DarkImage.data, CAMERA_WIDTH*CAMERA_HEIGHT*sizeof(unsigned short), 1, file);
			fclose(file);
			if ( Settings.BinMode > 1 )
			{
					// If binning is enabled, simulate binning on the dark image
				unsigned short* DarkPtr = (unsigned short*)DarkImage.data;
				unsigned int sum;
				int x, y, i, j;
				for ( y = 0; y < Settings.ImageHeight; y++ )
				{
					for ( x = 0; x < Settings.ImageWidth; x++ )
					{
						sum = 0;
						for ( i = 0; i < Settings.BinMode; i++ )
							for ( j = 0; j < Settings.BinMode; j++ )
								sum += DarkPtr[(y+j)*CAMERA_WIDTH + (x+i)];
						DarkPtr[y*Settings.ImageWidth+x] = sum / Settings.BinMode;
					}
				}
			}
			
		}
		else
		{
			delete[] DarkImage.data;
			DarkImage.data = NULL;
		}
	}
	
		// Allocate a static framebuffer for the camera image
	pixel_data = new unsigned short[Settings.ImageWidth*Settings.ImageHeight];
	if ( pixel_data != NULL )
		memset( pixel_data, 0, Settings.ImageWidth*Settings.ImageHeight*sizeof(unsigned short) );
	
#ifdef SDL_ENABLED
		// Initialize SDL
	ret = SDL_Init( SDL_INIT_VIDEO );
	printf( "SDL_Init returns %d\n", ret );
	screen = SDL_SetVideoMode( Settings.ImageWidth, Settings.ImageHeight, 32, SDL_HWSURFACE | SDL_DOUBLEBUF );
	SDL_WM_SetCaption( "CCD Camera", NULL );
	SDL_Event event;
	
	ret = TTF_Init();
	printf( "TTF_Init returned %d\n", ret );
	font = TTF_OpenFont( "consola.ttf", FONT_SIZE );
		//printf("Font at 0x%08X\n", (unsigned int)font);
	
	int bytes_pp = (Settings.BitMode == 12 ? 2 : 1);
	CCD_surface = SDL_CreateRGBSurfaceFrom( pixel_data, Settings.ImageWidth, Settings.ImageHeight,
																				 bytes_pp*8, Settings.ImageWidth*bytes_pp,
																				 0x0ff0, 0x0ff0, 0x0ff0, 0);
		//printf( "CCD_surface at 0x%08X\n", (unsigned int)CCD_surface );
#endif
	
	if ( Settings.ProcessImages && !Settings.ReadRaw && !Settings.ReadFits)
	{
			// Load Star Catalog from current time, if not reading from file
		SYSTEMTIME systime;
		GetSystemTime( &systime );
		
		glCatalog.Initialize( "catalog.txt", &systime, Settings.FocalLength );
	}
	
	Camera.Timer = 0;
	Camera.Initialized = false;
	Camera.EngineStarted = false;

	// Initialize camera engine, if NOT reading from file.
	if (!Settings.ReadRaw && !Settings.ReadFits){
		printf( "Initializing device.\n" );
		ret = BUFCCDUSB_InitDevice();
		printf( "There are %d devices.\n", ret );
		if ( ret <= 0 )
		{
			// don't give error if ReadRaw
			if (Settings.ReadRaw)
				printf( "No devices, but ReadRaw is on.\n" );
			// if ReadRaw is off, give an error
			else if (Settings.ReadFits) {
				printf( "No devices, but ReadFits is on.\n" );
			} else {
				printf( "Error: No devices.\n" );
				goto exit;
			}
		}
		Camera.Initialized = true;

		if (!Settings.ReadRaw && !Settings.ReadFits)
			printf( "Assuming device number 1.\n" );
		Camera.DeviceID = 1;
		ret = BUFCCDUSB_GetModuleNoSerialNo( Camera.DeviceID, Camera.ModuleNo, Camera.SerialNo );
		//check camera details. Skip is using ReadRaw
		if ( ret == 1 && !Settings.ReadRaw && !Settings.ReadFits)
		{
			printf( "Module number: %s\n", Camera.ModuleNo );
			printf( "Serial number: %s\n", Camera.SerialNo );
		}
		else if (!Settings.ReadRaw && !Settings.ReadFits) 
		{
			printf( "Error: Couldn't retrieve module or serial number.\n" );
			goto exit;
		}

		if (!Settings.ReadRaw && !Settings.ReadFits)
			printf( "Starting camera engine.\n" );
		ret = BUFCCDUSB_AddDeviceToWorkingSet( Camera.DeviceID );
		ret = BUFCCDUSB_StartCameraEngine( NULL, Settings.BitMode );
		if ( ret == 2 )
		{
			printf("Error: Tried setting to 12 bit mode, but not supported -- setting to 8 bit instead.\n");
			Settings.BitMode = 8;
		}
		else if ( ret == -1 )
		{
			// don't throw error if reading raw
			if (Settings.ReadRaw || Settings.ReadFits)
				printf("Camera engine did not start.\n");
			// otherwise throw error and exit.
			else {
				printf("Error: Camera engine did not start.\n");
				goto exit;
			}
		}
		Camera.EngineStarted = true;
		printf( "BUFCCDUSB_StartCameraEngine returns %d\n", ret );

		printf( "Changing camera settings.\n");
		ret = BUFCCDUSB_SetCameraWorkMode( Camera.DeviceID, CAMERA_MODE(Settings.TriggerMode) );
		ret = BUFCCDUSB_SetGains( Camera.DeviceID, Settings.RedGain, Settings.GreenGain, Settings.BlueGain );
		ret = BUFCCDUSB_SetCustomizedResolution( Camera.DeviceID, CAMERA_WIDTH, Settings.ImageHeight,
			( Settings.BinMode == 1 ? 0 : 0x81 + (Settings.BinMode-2) ), 4 );
		// exposure input is in 50 microsecond units
		ret = BUFCCDUSB_SetExposureTime( Camera.DeviceID, Settings.ExposureTime*1000/50 );

		printf( "Installing hookers, start grabbing.\n" );
		ret = BUFCCDUSB_InstallFrameHooker( FRAME_TYPE, FrameCallBack );
		printf( "BUFCCDUSB_InstallFrameHooker returns %d\n", ret );
		ret = BUFCCDUSB_InstallUSBDeviceHooker( CameraFaultCallBack );

		if ( Settings.TriggerMode == true )
			BUFCCDUSB_StartFrameGrab( GRAB_FRAME_FOREVER );
	} // end !Settings.ReadRaw

	short mouse_x = 0, mouse_y = 0;
	coordinates MouseCoords;
	
		// create file structure names, add format to top of data file
	SYSTEMTIME systime;
	GetLocalTime( &systime );

	if (Settings.DumpRaws || Settings.DumpData){
	sprintf( dir_name, "data_%04d%02d%02d_%02d%02d%02d", systime.wYear, systime.wMonth, systime.wDay,
					systime.wHour, systime.wMinute, systime.wSecond );
	mkdir( dir_name );
	}

	if (Settings.DumpFits){
	sprintf( fits_dir_name, "fits_data_%04d%02d%02d_%02d%02d%02d", systime.wYear, systime.wMonth, systime.wDay,
					systime.wHour, systime.wMinute, systime.wSecond );
	mkdir( fits_dir_name );
	}

	if (Settings.DumpData) {
		sprintf( output_name, "./%s/log_%04d%02d%02d_%02d%02d%02d.csv", dir_name, systime.wYear, systime.wMonth, systime.wDay,
			systime.wHour, systime.wMinute, systime.wSecond );
	}

	sprintf( output_foc, "./%s/foc_%04d%02d%02d_%02d%02d%02d.csv", dir_name, systime.wYear, systime.wMonth, systime.wDay,
					systime.wHour, systime.wMinute, systime.wSecond );
	
	if (Settings.DumpData) { 
		FILE* file = fopen(output_name, "w");
		if ( file != NULL ) {
			fprintf(file, "MM/DD/YYYY Time,LAT (deg),LONG (deg),ALT (m),"
						  "RA (deg),DEC (deg),RA (hr:min:sec),DEC_deg (deg:min:sec),"
						  "AZI (deg),ELE (deg),EleCorr (arcsec),"
						  "MeanSky,NumDetected,"
						  "RQuat0,RQuat1,RQuat2,RQuat3,Focal Length\n");
			fclose(file);
		}
	}
	
	printf( "Entering main loop.\n" );
	while ( !fault )
	{
			// Handle SDL drawing
#ifdef SDL_ENABLED
		/*
		 
		 printf( "(RA, DEC) = (%.6f, %.4f) | (AZI, ELE) = (%.4f, %.4f) | Stars: %d | Mean: %.1f\n"
		 "EleCorr: %d | LST: %f | HA: %f | (LAT,LONG) = (%.4f,%.4f) | ALT = %.1f\n",
		 */
		GUIFlipScreen();
		GUIPrintText( "%02d/%02d/%04d %02d:%02d:%02d.%03d", Image_LocalTime.wMonth, Image_LocalTime.wDay, Image_LocalTime.wYear,
								 Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds );
		GUIPrintText( "Mean Sky Level: %.1f", Mean_Sky_Level );
		GUIPrintText( "Stars Detected: %d", num_detected );
		GUIPrintText( "Latitude:  %.4f", Settings.Latitude );
		GUIPrintText( "Longitude: %.4f", Settings.Longitude );
		GUIPrintText( "Altitude:  %.1f", Settings.Altitude );
		if ( Image_Correct_ID )
		{
			GUIPrintText( "Ele Corr: %d", Ele_Corr );
			GUIPrintText( "Boresight: %.4f", Image_Bore*180/M_PI );
			GUIPrintText( "AZI: %.4f", Image_Azi*180/M_PI );
			GUIPrintText( "ELE: %.4f", Image_Ele*180/M_PI );
			GUIPrintText( "RA: %.6f", Image_Coords.RA*12/M_PI );
			GUIPrintText( "DEC: %.4f", Image_Coords.DEC*180/M_PI );
			//GUIPrintText( "Focal Length Used: %.1f", Prior.focal_length );
		}
		else
		{
			GUIPrintText( "Ele Corr: ERR" );
			GUIPrintText( "Boresight: ERR" );
			GUIPrintText( "AZI: ERR" );
			GUIPrintText( "ELE: ERR" );
			GUIPrintText( "RA:  ERR" );
			GUIPrintText( "DEC: ERR" );
			//GUIPrintText( "Focal Length Used: ERR" );
		}
		GUIPrintText( "" );
		GUIPrintText( "Mouse: (%d,%d)", mouse_x, mouse_y );
		if ( Image_Correct_ID )
		{
			GUIPrintText( "MouseAZI: %.4f", 0.0 );
			GUIPrintText( "MouseELE: %.4f", 0.0 );
			GUIPrintText( "MouseRA:  %.6f", MouseCoords.RA*12/M_PI );
			GUIPrintText( "MouseDEC: %.4f", MouseCoords.DEC*180/M_PI );
		}
		else
		{
			GUIPrintText( "MouseAZI: ERR" );
			GUIPrintText( "MouseELE: ERR" );
			GUIPrintText( "MouseRA:  ERR" );
			GUIPrintText( "MouseDEC: ERR" );
		}
		SDL_Flip(screen);
			// Handle SDL events
		while( SDL_PollEvent( &event ) )
		{
			switch (event.type)
			{
				case SDL_MOUSEMOTION:
						// Get Mouse coordinates and find RA,DEC
					mouse_x = event.motion.x;
					mouse_y = Settings.ImageHeight-1 - event.motion.y;
					double rP[3], r[3];
					GetSphericalFromImage( mouse_x - Settings.ImageWidth/2.0f + 0.5f, mouse_y - Settings.ImageHeight/2.0f + 0.5f, Settings.FocalLength, rP );
					matrix_mult( Image_Rt, rP, r );
					GetCoordsFromSpherical( r, &MouseCoords );
					break;
				case SDL_QUIT:
					goto exit;
				default:
					break;
			}
		}
#endif
		// If there is no camera, call FrameCallBackHelper directly. The BytePtr is created inside, so pass NULL.
		if (Settings.ReadRaw || Settings.ReadFits)
			FrameCallBackHelper(NULL);
		else{
			// The following is to let camera engine to be active..it needs message loop.
			if ( PeekMessage( &msg, NULL, 0, 0, PM_REMOVE ) )
			{
				if( msg.message == WM_QUIT )
				{
					goto exit;
				}
				else if ( msg.message == WM_TIMER )
				{
					TranslateMessage(&msg);
					DispatchMessage(&msg);
				}
			}
			if ( Settings.TriggerMode == false )
			{
				// If trigger mode is off, program is on a user-set timer.
				unsigned long long ticks = GetTickCount();	// milliseconds
				//printf( "TickCount: %llu\n", ticks );

				if ( Camera.Timer < 0 )
				{
					Camera.Timer = 0;
				}
				else
				{
					Camera.Timer += ticks - Camera.PrevTimeStamp;
				}
				Camera.PrevTimeStamp = ticks;

				if ( Camera.Timer >= Settings.CaptureDelay )
				{
					// reset timer
					Camera.Timer %= Settings.CaptureDelay;
					BUFCCDUSB_StartFrameGrab( 1 );
				}
			}
		} //if not reading from image
	} //while (!fault)
	
exit:
	printf( "Cleaning up.\n" );
	
#ifdef SDL_ENABLED
	TTF_CloseFont( font );
	TTF_Quit();
	SDL_Quit();
#endif
	
	BUFCCDUSB_InstallFrameHooker( FRAME_TYPE, NULL );
	BUFCCDUSB_StopFrameGrab();
	if ( Camera.EngineStarted )
		BUFCCDUSB_StopCameraEngine();
	if ( Camera.Initialized )
		BUFCCDUSB_UnInitDevice();
	
#ifdef SDL_ENABLED
	if ( CCD_surface != NULL )
		SDL_FreeSurface( CCD_surface );
#endif
	
	if (pixel_data != NULL)
		delete[] pixel_data;
	if (DarkImage.data != NULL)
		delete[] DarkImage.data;
	
	printf( "Press ENTER to exit.\n" );
	getchar();
	return 0;
}	
