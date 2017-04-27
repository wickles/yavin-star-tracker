This program encompasses the latest working version of our static Star Tracker.
Here, we summarize the directions for its application. Originally developed at
the University of California, Santa Barbara, by Gil Tabak, Alexander Wickes, and
Karanbir Toor, under the supervision of professor Philip Lubin. We thank our
associate Gary Hughes for useful discussions. This work is supported in part by
the NASA California Space Grant and by the Institute for Terahertz Science and
Technology (ITST) at UC Santa Barbara.

The program was designed and tested with a Mightex CCD camera (model CGE-B013-U)
and the Mightex SDK, but could be extended to work with other CCD cameras.
Captured images are processed for detection of possible star candidates. The
star candidates are then fed to a geometric voting algorithm that attempts to
identify them. Finally a quaternion-based least-squares fitting method
determines the celestial coordinates of the image center.

Before running the code, you might want to look at the settings file, which
controls various options. This includes how to use the program to read/write
the images used.

Binaries depend on MSVC runtime components.
Redistributables available from Microsoft:
MSVC-2010-x86: https://www.microsoft.com/en-us/download/details.aspx?id=5555
MSVC-2012: https://www.microsoft.com/en-us/download/details.aspx?id=30679

Mightex SDK and documentation available from:
http://www.mightexsystems.com/family_info.php?categories_id=141


////////////////////////////////////////////////////////////////////////

settings.txt file:

In here, the camera-based settings mode are included. The comments in the
settings file are mostly self-explanatory.

The focal length can be estimated from the specifications of the chip and lens
used. We computed this number more precisely for our lens using empirical
measurements.

We set the gains to their max level to maximize the signal for the camera. We
recommend leaving this option as is.

The altitude is used for the atmospheric correction, which is applied to the
azimuth and elevation.

DarkCoeff is the coefficient used to subtract the dark image. This reduces dark
current and other effects due to chip imperfections. We used 0.75, which seemed
to reduce the background artifacts for ~150ms. In general, DarkCoeff may need
to be adjusted for longer or shorter exposures in order to obtain correct. We
remark the relationship may not be linear between DarkCoeff and the exposure
time in general. The dark image used is called dark_image.raw. We used an
average of 30 different dark images. If this image is not present, a dark image
will not be used. Another dark image should be used for a different camera
chip, since the dark currents and other artifacts will be different. A
different dark image may need to be used in the case when the temperature
changes.

////////////////////////////////////////////////////////////////////////

Saving Images and algorithm output:

The program will generate an excel file log_YYYYMMDD_HHMMSS inside of a
directory data_YYYYMMDD_HHMMSS. This file keeps track of the algorithm output.
The letters in the name represent the date and time.

There is also an option to write RAW images as the camera runs (set DumpImages =
TRUE). Warning: this might take a lot of memory on your hard drive, as each
image consumes several megabytes. The images are placed in the directory
data_YYYYMMDD_HHMMSS. Each image is named image_YYYYMMDD_HHMMSS.raw. The program
will also generate a names.txt file with the image names, which is used by the
program to read the images, see below.

There was an issue with how the Mightext camera read the data, for which we had
to account by shifting bits in the data directly. In the images we write, this
issue is accounted for. However, we do NOT yet subtract the dark image. This
will have to be done if the files will be processed.

The RAW images in the directory may be used by the image reading feature. The
RAW images may also be opened in AstroArt.

To correctly open RAW files in AstroArt:

Go to File->Configure RAW, and make the following modifications:
Under Custom RAW, select Grayscale
Under Bits, select 16 Int
Under Word, select PC
Under Sign, select Unsigned,
Under File Extension, write *.raw
Under Header size, set 0
Under X,Y set 1280, 960.

Saving FITS Images:

There is also an option to write FITS images as the camera runs (set DumpFits =
TRUE). Warning: this might take a lot of memory on your hard drive, as each
image consumes several megabytes. The images are placed in the directory
data_YYYYMMDD_HHMMSS. Each image is named image_YYYYMMDD_HHMMSS.fits. The
program will also generate a names.txt file with the image names, which is used
by the program to read the images, see below.

Reading Images:

To do this, set ReadRaw = TRUE in the settings.txt file, and make sure to run
Star_Tracker_silent, which works via a command prompt. Star_Tracker_silent will
ask for a directory name. Input the name of the directory that the program
generated. It is important that the directory include the names.txt files to
know which images to read (and their order). This feature is only available in
the silent version, see below.

Also the file names must be written the same way the program generates them,
image_YYYYMMDD_HHMMSS. The code figures out the time of the image by parsing
the name of the file.

The program will subtract the dark image and execute the algorithm. It will also
output solutions into a new directory.

Hint: Instead of typing the directory name, it may be easier to turn on
QuickEdit mode. With QuickEdit mode on, you simply right click to paste text.

Reading FITS:

You can read in FITS files using this same procedure as reading raw images,
except set ReadFits = TRUE. The names.txt file and the image file names must
also be configured the same way as for RAW images, defined above. This feature
also only works in the silent version of the star tracker, see below.

Turning on QuickEdit mode:

Click the icon in the upper-left corner of the Command Prompt window, and then
click Properties. On the Options tab, click to select the QuickEdit Mode and
Insert Mode check boxes. Click OK.

////////////////////////////////////////////////////////////////////////

Calibration Mode:

A recent addition is the option to ask the computer to actively seek a focal
length. This is done in a two step process:

(1) find a focal length for which a solution is obtained

(2) find the best focal length once a series of image stars have been
identified. This is done by minimizing the sum of the errors of the image stars'
locations.

Both stages begin their search at the focal length specified in the settings.txt
file.

The maximum interval considered is determined by the quantity FocalError, which
is found in the settings.txt file. This quantity is the PERCENTAGE of the
quantity FocalLength.

The interval size of searching is determined by CaliTries, which is also
specified in the setttings.txt file.

If DumpData is TRUE, the program will also generate a file foc_YYYYMMDD_HHMMSS
in the data_YYYYMMDD_HHMMSS directory. This file is used to show the
relationship between the (sum of the) errors and the focal length. Each analyzed
solution has two lines: the focal lengths and the errors, in that order.

If DumpData is TRUE, a running average of the focal length will also be recorded
in the log_YYYYMMDD_HHMMSS file.


////////////////////////////////////////////////////////////////////////


Running the program on the sky:

Remember to adjust you camera lens by focusing, opening the shutter
sufficiently, and holding the camera steady. The application works in two
modes, silent and release. Also be wary of light pollution, as the star tracker
has inconsistent performance if there is a bright light source nearby such a
lamp.

Versions:

The star tracker comes in two versions, silent and release. The release version
offers a graphical user interface to view what the camera sees and displays the
algorithm information on screen. This version also read particular attitude
values by hovering the curser over different parts of the graphical user
interface.

The star tracker also come in a silent version. This version is less intensive
on the hardware and only uses the command line to communicate with the user.
This is the only version that supports the read files (RAW and FITS see above)
option.


Code and Class Structure::

Star_Tracker.cpp: This file does the necessary preprocessing, such loading the
setting from the settings.txt file, initializing the catalog, connecting the
camera, and obtaining the user input, if relevant. After all this, the code in
this file calls the various stages of the algorithm, then prints the results
appropriately.

ini_reader.h/.cpp: These parse and read from the settings.txt file.

Star.h/.cpp: These files define most of the constants used in the codes
algorithm, define many of the data structures used in various stages of the
algorithm, and implement many of the mathematical methods needed to perform the
star tracker algorithm.

Catalog.h/.cpp/.txt: These files define the catalog used by other parts of the
code to create paris of stars and identify image stars.

Detector.h/.cpp: These files correspond the stage of the algorithm that detects
the images stars from the star field image, centroids them and approximates the
associated error.

Identifier.h/.cpp: These files use multiple stages of the geometric voting
algorithm to map catalog stars to image stars.

Attitude.h/.cpp: These files take in identified image stars and return the
attitude of the camera based on the images.

Atmosphere.h.cpp: These files TODO.

////////////////////////////////////////////////////////////////////////

TODO:
- Compute mean for dark and image, scale dark accordingly.
