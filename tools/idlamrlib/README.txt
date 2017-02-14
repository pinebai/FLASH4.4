README and quick start guide for IDLAMRVis

Table of contents:

0.  Conventions
1.  Getting IDLAMRVis running
    1.1 Restarting in case of an error or crash
    1.2 Quick start with a specific plotfile
    1.3 Quick start with an AMR object or derived data
2.  The Root IDLAMRVis program
3.  Column-density visualization
    3.1 Rotation sliders
    3.2 Comments
4.  Slice visualization
    4.1 Slice Scalar Range Window
    4.2	Slice Geometry Window	
    4.3   Rendering Area   
5.  Line visualization
6.  Command-line access to visualization tools
    6.1 Column-density (the nice-but-slow type)
	6.1.1 Column-density aligned to an axis	 
	6.1.2.1 Lowest-level tool: COLUMN_ANY	
	6.1.2.2 Higher-level tool: COL_ANY_FRAME
    6.2 Slices
	6.2.1 Slices aligned with an axis
	6.2.2 Slices along any angle
    6.3 Lines
7.  Making movies
    7.1 Rendering your frames
    7.2 Encoding the movie
	7.2.1 Movies in time



0.  Conventions
Primary mouse button refers to the left button on a right-handed mouse, or the right button on a left-handed mouse. It is also refered to as MB0.
Second mouse button refers to the middle mouse button.  If you don't have a 3-button mouse, this won't work for you. It is also refered to as MB1
Ternary mouse button refers to the right mouse button on a right-handed mouse, or the left mouse button on a left-handed mouse.  It is also refered to as MB2.
AMR Object refers to a variable in IDL of an AMR component.  They are created by GET_AMRCOMPONENT and can be manipulated by a number of functions in the idlamrlib distribution.

1.  Getting IDLAMRVis running

Checkout the idlamrlib package from CVS.
In your system shell type:
  cd idlamrlib
  idl
In IDL type:
  .com idlamrvis  (this load a bunch of functions)
  idlamrvis       (this starts the program)
IDLAMRVis will start up with a file selection dialog.
Select the plotfile you want to look at and then click Okay
IDLAMRVis will load the 'density' component selected plotfile.
To use different components, see section 1.2
If you click cancel or select something other than a plotfile when you click OK, it will
cause an error and you will have to restart idlamrvis:

1.1   Restarting due to error or crash
To restart idlamrvis in case of an error:
In the IDLAMRVis window click the 'X' in the title bar
then in IDL type idlamrvis

1.2    Quick start with a specific plotfile
In stead of typing idlamrvis to start, in IDL type
idlamrvis,filename='PLOTFILE'
Example:  idlamrvis,filename='/gpfs/ux453850/plt0195'
IDLAMRVis will start up and bypass the file selection dialog
The selected component will be density unless you specify otherwise
To select a componenet other than density, type
idlamrvis,filename='PLOTFILE',variable='VAR'
Example:  idlamrvis,filename='/gpfs/ux453850/plt0195',variable='pressure'

1.3	  Quick start with an AMR object or derived data
If you have an AMR object loaded into IDL's memory already, you can load it into IDLAMRVis.
This is also the way to vidualize data from some derived quantity.  Load your AMR object and
perform any operations on it you'd like, then import the AMR object into IDLAMRVis as follows:
.com idlamrvis
filename = '/gpfs/ux453850/plt0195'
component = 'density'
amr = get_amrcomponent(filename, component)
; Your data manipulation goes here
idlamrvis,amr=amr

2.	The Root IDLAMRVis Program
Once you have loaded data into IDLAMRVis, you are presented with the root of the GUI.
Here are the features:
Close:  Closes IDLAMRVis root window and all windows that it created
Column Density:  Opens the Column Density visualization tool.  This may take several minutes to start up.  If you use this tool with a variable other than 'density', you will get strange results.  (ex. What does the column-density of a velocity fiels look like?)
Slice:  Opens the Slice visualization tool.  Fairly Quick start up.
Status Field:  Says what the main program is doing.
Window Manager close button:  If IDLAMRVis has an error, the Close button will not work.  Close the program like you would any other misbehaving X-Windows program.

3.  Column-density visualization
This tool shows you Column-densities at any angle.  It may take a couple minutes to start up.  It renders in 2 modes, quick-and-dirty, and nice-but-slow.
Once it finishes loading, you will see several controls:
Close:  Closes the Column-density tool only
ReCalc Voxel:  Does nothing.  Some day it will rebuild the voxel for doing quick-and-dirty renderings at other zooms and centerings
Render:  Produces a nice-but-slow rendering of column-density based on the current view parameters.
3.1	 Rotation sliders
This program defines the orientation in terms of 3 rotations about fixed axes.  The view plane is defined at being perpendicular to the Z-axis.  The rotation is defined as follows:
First, rotate from X to Y by an angle Azimuth
Next, rotate from Y to Z by an angle Altitude
Finally, rotate from X to Y by an angle Pitch
3.2	 Comments
The Column-density tool was the first one I wrote for this visualization package.  As such, it displays my most primitive understanding of how to build a GUI in IDL.  This is not particularly user friendly or pretty and I intend to improve it.  If you need to view column-densities at a higher zoom, you have to do that manually at the command line.  See section 6 for command-line access to vis tools.

4.  Slice Visualization
This tool allows viewing slices of AMR data at arbitrary angles and zooms, as well as the ability to potput that data to a PPM or JPEG file, and also to view data along a 1-D line through the slice.
Two bugs should be mentioned.  The slice tool occasionaly produces jagged seams along the edge of a FAB.  This has been a tricky bug to track down.   Also, if you slice a AMR component that has negative values, It currently crashes because it tries to take the log of the data.  I will add a check to ban Log mode on negative data.
This tool is fairly well developed.  Here are the features in detail:
Close:  Closes the slice tool and any windows it spawned
Write JPEG:  Opens a file selection dialog.  Pick a filename and it will write a JPEG inage of the current window contents.
Write PPM:  Same as Write JPEG, but writes a PPM image, which is a larger file but is a lossless bitmap (higher quality).

4.1   Slice Scalar Range Window
This window houses controls for how the rendered data should be displayed.  It's features are:
Min & Max sliders: Allows you to select the max and min of the colormap for the data being displayed.  The absolute limits of the slider are the global max and min of the data set.  
Min & Max input boxes: Allows precise control of the max and min of the colormap.  Hit ENTER to set the new value.
Number of coutours:  Sets the number of coutour lines to be applied to the data plot.  Set 0 to turn them off.  The major ticks of the colormap bar are equal to the values of the contour lines.  Using contour lines may exxagerate the bug that causes seams along FAB boundaries because the coutour lines catch the seams.
Log/Linear button:  To switch between viewing the data as log or linear.  The large, bold word is the current mode. 

4.2	   Slice Geometry Window
This window houses controls for the projection of the slice.  There is the zoom control, the center control, and the angle control.
Log of Zoom:  The slider control the zoom of the slice.  The slider operates exponentially so as to give good control at all scales.  The total zoom range is from 0.10 to 4R, where R is the length ratio of the level-0 cell to the highest-level cell.  There is a text input box to set the zoom value exactly.  Higher zooms take longer to renger due to inefficencies in the rendering algorithm.
X,Y & Z Centers:  These control the center point of the rotations.  They can be dragged over the entire domain of the simulation, or set by the text boxes.  If your slice seems coarse, try setting to center to a point covered by high-level data.
Angle Sliders:  These control the angle of the slice plane.  For a more detailed explanation of the rotations used, see section 3.1.

4.3   Rendering Area
The area the lice is drawn into isn't just for pretty pictures.  A number of commands are available with the mouse in the rendering area.
Primary mouse button (MB0):
In the rendering area, drag with MB0 to recenter the slice.  This does not affect the rotation, only the offset of the viewing window relative to the center of the rotation.  When you release MB0, the view will re-render.
Double-click MB0 to make the view center on the point double-clicked.
Secondary mouse button (MB1):
In the rendering area, drag MB1 to view a 1-D sample of the data along the straight line from the start of the drag to the end of the drag.
Ternary mouse button (MB2):
Dragging MB2 will one day allow rotation of the slice plane, but I need to work out some of the math first.

5. Line Visualization
This is a new feature which allows viewing the value of AMR data on an arbitrarily-oriented 1-D line.  It is accessed by dragging the secondary mouse button over the rendering window of the Slice tool.  The line visualization feature soule be considered bare-bones at this time.

6.  Command-line access to visualization tools
All types of visualization in the IDLAMRVis suite are rendered by tools that can be called at the IDL prompt, if you know the syntax.  
6.1 Column-density (the nice-but-slow type):
6.1.1 Column-density aligned to an axis
For column-density aligned to an axis, the fastest tool is Mark's old column_amr.pro.  
The syntax is: column_amr,amr,plane 
COLUMN_AMR takes as input an AMR Object (from GET_AMRCOMPONENT), and the view axis.  Many other options are available and can be fairly-well explained from the header of column_amr.pro
6.1.2 Column-density along any angle
For column-density along an arbitrary angle, the tools are COL_ANY_FRAME and COLUMN_ANY.  
Both take an AMR Tree as input.  An AMR tree is a hierarchical representation of AMR date which can be made by calling AMR_TREE with the syntax:
tree = amr_tree(amr)
6.1.2.1 Lowest-level tool: COLUMN_ANY
COLUMN_ANY is the lower-level tool which simply returns a 2-D array representing the column-density.  The syntax of COLUMN_ANY is:
image = column_any(TREE,ANGLES)
ANGLES is a dblarr(3) of [azimuth,altitude,pitch] as described previously.  ANGLES are measured in radians
Other options of interest:
xrange,yrange:  These control the area that will show up in the image.  Each are dblarr(2) variables
zrange: This dblarr(2) controls the depth of the column.  If you make the span on zrange really small, you will get bad output.
subsample: This makes column_any generate a larger image and then average it down to a smaller one.  It can give prettier output, but can also make column_any run MUCH more slowly.
imgdim:  This controls the size of the output image.  It is an intarr(3).  An example would be [512,512,1].  It can be non-square, but imgdim and [xrange,yrange] should have the same aspect ratio or your output will be garbage.
6.1.2.2 Higher-level tool: COL_ANY_FRAME
COL_ANY_FRAME is intended to be a full-featured replacement for the older COLUMN_AMR tool.  
The minimal syntax is:
col_any_frame(tree=TREE, angles=ANGLES)
where TREE and ANGLES are as explained in 6.1.2.1
COL_ANY_FRAME supports the options xrange,yrange,zrange in the same way as COLUMN_ANY.
Other options of COL_ANY_FRAME are:
colorbar: Set true to draw the colorbar
scalerrange: Set tomanually specify the max and min of the colortable, ex: scalerrange=[1.22e-20,8.51e-12]
outM: What type of output do you want?  Available are 'X' for X Windows, 'PS' for postscript, 'PPM' for a PPM image, and 'JPEG' for a JPEG image
filename: The filename to write to if outM is set to 'PS','PPM',or 'JPEG'
aftersmooth: If you want to call IDL's SMOOTH function on the image after rendering, aftersmooth specifies the width of the smoothing window.  Acceptable values are odd numbers >= 3.

6.2 Slices
6.2.1 Slices aligned with an axis
For a slice oalogn an axis, use Mark's raster_amr.  The syntax is:
raster_amr,AMR,PLANE,SLICE
where AMR is an AMR Object, PLANE is the depth of the slice plane, and SLICE is the number of the axis the slice is perpendicular to, (0,1 or 2)
6.2.2 Slices along any angle
For slices along an arbitrary angle, call SLICE_ANY.  The syntax is:
image = slice_any(TREE, ANGLES)
Where TREE is from TREE = amr_tree(AMR) and angles follows the [azimuth,altitude,pitch] convention
Other options are:
target: An [x,y,z] specifying the genter of the slice-plane's rotation
xrange,yrange: The form [xmin,xmax],[ymin,ymax], specifying the offset and size of the region of the slice to render
imgdim: Of the form [x,y] specifying the pixel size of the output image.
SLICE_ANY is not as full-featured as RASTER_AMR.  If only produces raw images.  For pretty finished output, use the Slice tool of IDLAMRVis, or bug the author to make a more full-featured command-line version

6.3 Lines
Abritrary 1-D lines through AMR data are rather new and imcomplete.  The command-line tool is AMR_LINE.  The syntax is:
line = amr_line(amr=AMR, vector=VECTOR, point=POINT, range=RANGE)
AMR is an AMR object.
VECTOR is a normalized (or not) vector which the line is paralell to
POINT is a point in the AMR domain that the line goes through
RANGE is how far the line extends forwards and back from POINT.
As an example, if VECTOR=[1,0,0], POINT=[0,0,0] and RANGE=[-1,3] then the line would extand from [-1,0,0] to [3,0,0]
Other options are:
maxlevel: Specifies the maximum level to sample at
numpoints: Hom many points along the line to sample.  (Clamped at 64k points for speed reasons).
amr_line returns a 1-D dblarr().  It does not render.  

7.  Making movies
The basic way to make a movie of AMR data is to render a lot of frames (perhaps 1000).  Then take those frames to something called an Encoder, which turns them into a standard movie format.  If you are making a movie in time, you will probably want to record what time in the simulation (i.e. Millions of years) the frame came from.  If you don't have the time index of all the frames, a movie in time will not look good because the time step is jumping around.
7.1 Rendering your frames
Use any rendering tool you like to make column-densities slices, whatever.  The thing to remember is that you want to output in either PPM or JPEG format.  PostScript is nice for journals, but useless for movie making.  The way I made frames for movies is with a movie script.  One example is the routine RASTER_ZOOM_MOVIE.  Also see plotting.txt for lots of examples.  In plotting.txt, I made arrays out of the parameters for the movie and then iterate over the arrays.  Of particular interest is the main loop instruction: "for i=cpuid,999,ncpu do begin"  Rendering 1000 frames can take a long time, but is paralellizable.  I couldn't find a way to get IDL to use multiple CPUs on the same machine to speed up unless I ran seperate instances of IDL.  When running on a machine like meso.berkeley.edu, I would set ncpu=6 to use 3/4 of the machine, and run 6 copies of IDL, each one with a different value of cpuid (0 through 5).  This way the frames rendered independantly in paralell.  If I had split up the 1000 frames into 6 linear chunks, I would get worse "SpeedUp" because different parts of the same movie often render at different speeds.
7.2 Encoding the movie
Once you have your 1000 frames, you need to put them into a standard movie format.  I have been using the MPEG2 format because while it is not the most efficient, it is very widely compatible and has a lot of free tools available.  I used the free encoder ppmtompeg, part of the Berkeley MPEG Tools.  If it's not on the astron system, you can probably install it without too much difficulty.  Making the inputs file for it is a little tricky.  I've put a sample in movie_inputs.txt  This version uses JPEG frames.  For a version using PPM frames, it's a little simpler.  See 'man ppmtompeg'.
7.2.1 Movies in time
Making a movie in time is the same as making a movie in position with one caveat:  Your frames are probably not all spaced evenly in time.  You need them to be evenly spaced ot else the movie will speed up and slow down and be rather confusing.  To solve just this problem I wrote a little C program called smartframelinker.  In smartframelinker.tar are all the source files.  For instruction on how to use, put it in a new directory, untar it, make it, then do "smartframelinker --help".  The gist of it that if you say "smartframelinker --f1000  /directory_with_the_frames/* /directory_with_the_plotfiles/*" it will produce 1000 symbolic links to the frames, with the number of links to a frame proportional to the length of the timestep of the plotfile that produced that frame.  You can then feed these links to you movie-encoding program and it will make a movie with a constant time scale.  Another point:  smartframelinker does not require whole plotfiles, only the directory and Header.  Smartframelinker can also be used to stretch out a movie if you don't have enough frames and your encoding program will not let you use a very low FPS number.

