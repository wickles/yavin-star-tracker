/* Star_Tracker.cpp
Star Tracker by Alexander Wickes and Gil Tabak
September 20, 2010

This program connects to a Mightex CCD Camera, model CGE-B013-U, in order to capture
images. These images are processed (see detector.cpp) in order to detect possible
stars and acquire their local spherical coordinates, in the image's reference frame.
A star catalog is pre-loaded (see catalog.cpp) and a geometric voting algorithm is
used to identify the image stars (see identifier.cpp). We then use a quaternion-based
least squares regression method to find the coordinates of the center of the image
(see attitude.cpp).

Settings are loaded from settings.txt. SDL (Simple DirectMedia Layer) is used for the
GUI. 
*/

/* TODO
	- Compute mean for dark and image, scale dark accordingly.
	*/

#include "stdafx.h"

#include <cstdio>
#include <cstdarg>
#include <direct.h>
#include "ini_reader.h"
#include "detector.h"
#include "catalog.h"
#include "identifier.h"
#include "attitude.h"
#include "timer.h"
#include "atmosphere.h"

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
char output_name[512];

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
	bool DumpImages;
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
	float Altitude;	// ???
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
	FILE* file = fopen( filename, "r" );
	if ( file == NULL )
		printf( "Error opening settings file.\n" );

	Settings.TriggerMode =	IniGetBool( file, "TriggerMode", false );
	Settings.DumpImages =	IniGetBool( file, "DumpImages", false );
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
	Settings.Altitude =	IniGetFloat( file, "Altitude", 0.0f );

	if ( file != NULL )
		fclose(file);
	
	/*
	printf( "TriggerMode: %s\n", (Settings.TriggerMode ? "True" : "False") );
	printf( "DumpImages: %s\n", (Settings.DumpImages ? "True" : "False") );
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
	*/

	return 1;
}

// Callback function for handling 
void FrameCallBack( TProcessedDataProperty* Attributes, unsigned char* BytePtr )
{
	if ( Attributes->CameraID == Camera.DeviceID ) // The working camera.
	{
		SYSTEMTIME systime;
		GetLocalTime( &systime );
		/*
		printf( "Grabbing frame: %04d.%02d.%02d : %02d:%02d:%02d.%03d\n", systime.wYear, systime.wMonth, systime.wDay,
							systime.wHour, systime.wMinute, systime.wSecond, systime.wMilliseconds );
							*/
		Image_LocalTime = systime;

		//printf( "BytePtr at %08X\n", (unsigned int)BytePtr);

		if ( BytePtr == NULL )
			return;

		// If in 12 bit mode, fix the messed up bit ordering.
		if ( Settings.BitMode == 12 )
		{
			int i;
			unsigned short* BytePtr16 = (unsigned short*)BytePtr;
			for ( i = 0; i < CAMERA_WIDTH*Settings.ImageHeight; i++ )
				BytePtr16[i] = (BytePtr[2*i] << 4) | BytePtr[2*i+1];
		}

		if (Settings.DumpImages)
		{
			// Save image to file w/ timestamp
			char buffer[512];

			// get system timestamp (millisecond accuracy)
			sprintf( buffer, "./%s/image_%04d%02d%02d_%02d%02d%02d%03d.raw", dir_name, systime.wYear, systime.wMonth, systime.wDay,
							systime.wHour, systime.wMinute, systime.wSecond, systime.wMilliseconds );

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
		
		int bytes_pp = (Settings.BitMode == 12 ? 2 : 1);
		unsigned short* BytePtr16 = (unsigned short*)BytePtr;
		// If Binning is enabled, average pixels horizontally (camera only bins vertically)
		if ( Settings.BinMode > 1 )
		{
			unsigned int sum;
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
		// Copy temporary framebuffer to static framebuffer for display in GUI
		memcpy( pixel_data, BytePtr, Settings.ImageWidth*Settings.ImageHeight*bytes_pp );

		if ( Settings.ProcessImages )
		{
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

			if ( glCatalog.initialized )
			{
				detector_s Detector = {
					GL_SAMPLE_SKIP,
					LOCAL_WIDTH,
					LOCAL_HEIGHT,
					LOCAL_SAMPLE_SKIP,
					STAR_MIN_OUTSTND
				};
				image_s Image = {
					Settings.ImageWidth,
					Settings.ImageHeight,
					Settings.BitMode,
					Settings.FocalLength,
					BytePtr,
					NULL
				};

				int num_stars = DetectStars(ImageStars, &Detector, &Image);
				//printf( "	Number of Stars Detected: %d\n", num_stars );
				Mean_Sky_Level = Detector.mean_sky;

				if (ImageStars.size() >= 3)
				{
					// enough stars, continue with identification
					IdentifyImageStars(ImageStars, glCatalog, ( Prior_Valid ? &Prior : NULL ));

					coordinates Coords;
					sfloat R[3][3];
					double rms_error = GetAttitude(ImageStars, glCatalog.Stars, &Coords, R);
					Image_Correct_ID = ( rms_error <= 0.01 );

					if ( Image_Correct_ID )
					{
						Prior.coords = Coords;
						Prior.tickCount = GetTickCount();
						Prior_Valid = true;

						Image_Coords = Coords;

						// Get Az, El
						double relJDN = getJulianDate( &systime );
						double LST, HA;
						LST = getLST( Coords.RA, Coords.DEC, relJDN );
						HA = getHA1( LST, Coords.RA );
						Image_Azi = getAzi( HA, Settings.Latitude, Coords.DEC );
						Image_Ele = getEle( HA, Settings.Latitude, Coords.DEC );


						int i, j;
						for ( i = 0; i < 3; i++ )
							for ( j = 0; j < 3; j++ )
								Image_Rt[i][j] = R[j][i];
					}
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
								//fprintf(file, "%02d/%02d/%04d %02d:%02d:%02d.%03d | (RA, DEC) = (%.15f, %.15f) = (%02d:%02d:%f, %02d:%02d:%f) | Rt = { { %f, %f, %f }, { %f, %f, %f }, { %f, %f, %f } }\n",
								fprintf(file, "%02d/%02d/%04d \
											   %02d:%02d:%02d.%03d,\
											   %.15f,%.15f,\
											   %02d:%02d:%f,\
											   %02d:%02d:%f,\
											   %f,%d,\
											   %f,%f,%f,%f,%f,%f,%f,%f,%f\n",
										systime.wMonth, systime.wDay, systime.wYear,
										systime.wHour, systime.wMinute, systime.wSecond, systime.wMilliseconds,
										Coords.RA, Coords.DEC,
										coords.RA_hr, coords.RA_min, coords.RA_sec,
										coords.DEC_deg, coords.DEC_min, coords.DEC_sec,
										Detector.mean_sky, num_stars,
										Image_Rt[0][0], Image_Rt[0][1], Image_Rt[0][2],
										Image_Rt[1][0], Image_Rt[1][1], Image_Rt[1][2],
										Image_Rt[2][0], Image_Rt[2][1], Image_Rt[2][2]);
							}
							else
							{
								fprintf(file, "%02d/%02d/%04d \
											   %02d:%02d:%02d.%03d,\
											   FAIL,FAIL,FAIL,FAIL,\
											   %f,%d,\
											   FAIL,FAIL,FAIL,FAIL,FAIL,FAIL,FAIL,FAIL,FAIL\n",
										systime.wMonth, systime.wDay, systime.wYear,
										systime.wHour, systime.wMinute, systime.wSecond, systime.wMilliseconds,
										Detector.mean_sky, num_stars );
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
	}
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

	// Setup Priors
	Prior.coords.RA = 0.0f;
	Prior.coords.DEC = 0.0f;
	Prior.tickCount = GetTickCount();
	Prior_Valid = true;

	if ( Settings.BitMode != 12 )
	{
		printf( "Error: Program only supports 12 bit images.\n" );
		return 0;
	}

	Prior.phi_max = MAX_PRIOR_DIST;

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
	SDL_Surface* screen = SDL_SetVideoMode( Settings.ImageWidth, Settings.ImageHeight, 32, SDL_HWSURFACE | SDL_DOUBLEBUF );
	SDL_WM_SetCaption( "CCD Camera", NULL );
	SDL_Event event;

	ret = TTF_Init();
	printf( "TTF_Init returned %d\n", ret );
	TTF_Font* font = TTF_OpenFont( "consola.ttf", FONT_SIZE );
	//printf("Font at 0x%08X\n", (unsigned int)font);

	int bytes_pp = (Settings.BitMode == 12 ? 2 : 1);
	SDL_Surface* CCD_surface = SDL_CreateRGBSurfaceFrom( pixel_data, Settings.ImageWidth, Settings.ImageHeight,
													bytes_pp*8, Settings.ImageWidth*bytes_pp,
													0x0ff0, 0x0ff0, 0x0ff0, 0);
	//printf( "CCD_surface at 0x%08X\n", (unsigned int)CCD_surface );
#endif
	
	if ( Settings.ProcessImages )
	{
		// Load Star Catalog
		SYSTEMTIME systime;
		GetSystemTime( &systime );

		glCatalog.Initialize( "catalog.txt", &systime, Settings.FocalLength );
	}

	Camera.Timer = 0;
	Camera.Initialized = false;
	Camera.EngineStarted = false;

	// Initialize camera engine
	printf( "Initializing device.\n" );
	ret = BUFCCDUSB_InitDevice();
	printf( "There are %d devices.\n", ret );
	if ( ret <= 0 )
	{
		printf( "Error: No devices.\n" );
		goto exit;
	}
	Camera.Initialized = true;

	printf( "Assuming device number 1.\n" );
	Camera.DeviceID = 1;
	ret = BUFCCDUSB_GetModuleNoSerialNo( Camera.DeviceID, Camera.ModuleNo, Camera.SerialNo );
	if ( ret == 1 )
	{
		printf( "Module number: %s\n", Camera.ModuleNo );
		printf( "Serial number: %s\n", Camera.SerialNo );
	}
	else
	{
		printf( "Error: Couldn't retrieve module or serial number.\n" );
		goto exit;
	}

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
		printf("Error: Camera engine could not start.\n");
		goto exit;
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

	short mouse_x = 0, mouse_y = 0;
	coordinates MouseCoords;

	// create file structure names, add format to top of data file
	SYSTEMTIME systime;
	GetLocalTime( &systime );
	sprintf( dir_name, "data_%04d%02d%02d_%02d%02d%02d", systime.wYear, systime.wMonth, systime.wDay,
			systime.wHour, systime.wMinute, systime.wSecond );
	mkdir( dir_name );
	sprintf( output_name, "./%s/log_%04d%02d%02d_%02d%02d%02d.csv", dir_name, systime.wYear, systime.wMonth, systime.wDay,
			systime.wHour, systime.wMinute, systime.wSecond );
	FILE* file = fopen(output_name, "w");
	if ( file != NULL )
	{
		fprintf(file, "Month/Day/Year Hour:Minute:Second:Millisecond,RA,DEC,RA_hr,DEC_deg,Mean_Sky,Num_Detected,Rt00,Rt01,Rt02,Rt10,Rt11,Rt12,Rt20,Rt21,Rt22\n");
		fclose(file);
	}

	printf( "Entering main loop.\n" );
	while ( !fault )
	{
		// Handle SDL drawing
#ifdef SDL_ENABLED
		SDL_BlitSurface( CCD_surface, NULL, screen, NULL );
		int yText = 0;
		DrawSDLText( font, screen, 0, yText, "%02d/%02d/%04d %02d:%02d:%02d.%03d", Image_LocalTime.wMonth, Image_LocalTime.wDay, Image_LocalTime.wYear,
											Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds );
		yText += FONT_SIZE;
		DrawSDLText( font, screen, 0, yText, "Mean Sky Level: %f", Mean_Sky_Level );
		yText += FONT_SIZE;
		if ( Image_Correct_ID )
		{
			DrawSDLText( font, screen, 0, yText, "Camera (RA, DEC) = (%.10f, %.10f)", Image_Coords.RA, Image_Coords.DEC );
			yText += FONT_SIZE;
			DrawSDLText( font, screen, 0, yText, "Camera (AZI, ELE) = (%.5f, %.5f)", Image_Azi, Image_Ele );
		}
		else
			DrawSDLText( font, screen, 0, yText, "Camera (RA, DEC) = INVALID" );
			yText += FONT_SIZE;
			DrawSDLText( font, screen, 0, yText, "Camera (AZI, ELE) = INVALID" );
		yText += FONT_SIZE;
		if ( Image_Correct_ID )
			DrawSDLText( font, screen, 0, yText, "Mouse (%d, %d), (RA, DEC) = (%.10f, %.10f)", mouse_x, mouse_y, MouseCoords.RA, MouseCoords.DEC );
		else
			DrawSDLText( font, screen, 0, yText, "Mouse (%d, %d), (RA, DEC) = INVALID", mouse_x, mouse_y );
		yText += FONT_SIZE;
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
#else
		/*
		printf( "%02d/%02d/%04d %02d:%02d:%02d.%03d", Image_LocalTime.wMonth, Image_LocalTime.wDay, Image_LocalTime.wYear,
											Image_LocalTime.wHour, Image_LocalTime.wMinute, Image_LocalTime.wSecond, Image_LocalTime.wMilliseconds );
											*/
		if ( Image_Correct_ID )
			printf( "MeanSky: %f | (RA, DEC) = (%.10f, %.10f) | (AZI, ELE) = (%.5f, %.5f)\n",
					Mean_Sky_Level, Image_Coords.RA, Image_Coords.DEC, Image_Azi, Image_Ele );
		else
			printf( "MeanSky: %f | (RA, DEC) = INVALID | (AZI, ELE) = INVALID\n", Mean_Sky_Level );
#endif

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
	}

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
