/*
 * Sameed Siddiqui
 * EE 113DA Winter-Spring 2016
 *
 * Fracture detection project
 * 		This program is designed to give as a proof-of concept that fractures can be
 * 		detected using some combination of techniques. We illustrate a very rough
 * 		fracture detection program. First, we implement a canny edge detector, which
 * 		detects edges with higher accuracy than the sobel filter. Next we use OpenCV's
 * 		probabilistic hough transform to detect lines in the filtered image. Noting
 * 		that un-broken bones generally look like parallel cylinders, we define a
 * 		broken bone based on the intersections of cylinders.
 *
 *		NOTE: This program was built specifically for the TI LCDK c6748 DSP chip.
 *			  Due to some quirks in this chip, as many variables as possible had to be
 *			  declared in a global scope, and the code itself could not be split into 
 *			  multiple .c files.

/* ==========================================================================
/*                          INCLUDE FILES
/* ========================================================================== */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <xdc/std.h>
#include <xdc/runtime/Error.h>
#include <xdc/runtime/System.h>
#include <ti/sysbios/BIOS.h>
#include <ti/sysbios/knl/Clock.h>
#include <ti/sysbios/knl/Task.h>
#include <ti/sysbios/knl/Semaphore.h>
#include <ti/sysbios/family/c64p/Hwi.h>



#include <stdio.h>
#include <math.h>
#include "string.h"
#include <stdlib.h>
#include "psc.h"
#include "vpif.h"
#include "raster.h"
#include "interrupt.h"
#include "lcdkC6748.h"
#include "soc_C6748.h"
#include "hw_psc_C6748.h"
#include "adv7343.h"
#include "tvp5147.h"
#include "cdce913.h"
#include "codecif.h"
#include "i2cgpio.h"
#include "cv.h"
#include "cxtypes.h"
#include "facedetect.h"
#include "L138_LCDK_switch_led.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*
* NOTE: This program was built specifically for the TI LCDK c6748 DSP chip.
*		Due to some quirks in this chip, as many variables as possible had to be
*	    declared in a global scope, and the code itself could not be split into
*	    multiple.c files.
*/

#define height_rows 150
#define width_cols 250

int demonstration_number = 1; 	// change necessary parameters for the different trials
								// trial 1: fracture
								// trial 2: no fracture


////////////////////////////////
//// Variables for Canny filter
////////////////////////////////
int upperThreshold = 120;	// Gradient strength nessicary to start edge
int lowerThreshold = 30;	// Minimum gradient strength to continue edge
int num_hysterisis_loops = 8; // Number of loops to run hysterisis for


unsigned char greyscale_image[height_rows][width_cols];
// #pragma DATA_SECTION(greyscale_image, ".EXT_RAM")

int gradients[height_rows][width_cols];
// #pragma DATA_SECTION(gradients, ".EXT_RAM")

int edge_dir[height_rows][width_cols];
// #pragma DATA_SECTION(edge_dir, ".EXT_RAM")

unsigned char* temp_image;
int next_edge_direction[2][2];
double gaussianMask[5][5];	// Gaussian mask
int GxMask[3][3];		// Sobel mask in the x direction
int GyMask[3][3];		// Sobel mask in the y direction
int binarize_upper_limit = 240;
int repeat_count = 0;
double temp_value_filter = 0;
double gx_temp;
double gy_temp;
int thisAngle;
int newAngle;
int temp_var_max_gradient = 0;
//// END VARIABLES FOR CANNY FILTER


////////////////////////////////
//// general globals
////////////////////////////////
int k = 0;
int i = 0;
int j = 0;

int row_counter;
int col_counter;
int counter;
//// END general globals


//////////////////////////////
// VARIABLES FOR Visual studio
//////////////////////////////
char str[10];
char buf[10];
char* fgets_input;
unsigned char bitmap[height_rows*width_cols];
FILE *ptr_file;
//// END VARIABLES FOR Visual studio


//////////////////////////////
//// Variables for Hough
//////////////////////////////

IplImage *ss_image=0; // declare an IplImage structure pointer
CvSeq* lines = 0;
CvPoint* line1;
CvPoint* line2;
float ss_m1; //slope
float ss_b1; // intercept
float ss_m2; //slope
float ss_b2; // intercept
float ss_x_intersection;
float ss_y_intersection;
float ss_intersection_distance;
float ss_max_intersection_square_dist = 100;
//// END Variables for Hough



void canny(unsigned char* bitmap)
// Canny edge detection
// this is better edge detection than a simple Sobel Filter.
// 1. First you start off with Gaussian filter to remove xray noise
// 2. Use the Sobel filter to detect all edges in the photo.
// 3. Find the edge magnitudes and find the edge directions
//    -> each edge direction/angle can be binned into 0, 45, 90, or 135 degrees
// 4. Keep only those edges with magnitudes above a certain threshold
// 5. Add to your images those edges which are above a lower threshold but are
//    in the same direction as and adjacent to the stronger edges.
// 6. binarize the image.
{
	// INSTRUCTIONS
	// For LCDK:
	// (1) uncomment "m_free(bitmap)". This frees up memory
	//		(a) However, this is only really effective if you dynamically allocate
	//			the other image matrices (e.g. greyscale_image) because otherwise
	//			freeing up memory wont do you any real good because at some point
	//			in your program you'll have used bitmap, greyscale_image, and
	//			temp_image
	// Note that some part of code is adapted from http://www.pages.drexel.edu/~nk752/cannyTut2.html

	// Define Gaussian masks
	gaussianMask[0][0] = 0.003765;		gaussianMask[0][1] = 0.015019;		gaussianMask[0][2] = 0.023792;		gaussianMask[0][3] = 0.015019;		gaussianMask[0][4] = 0.003765;
	gaussianMask[1][0] = 0.015019;		gaussianMask[1][1] = 0.059912;		gaussianMask[1][2] = 0.094907;		gaussianMask[1][3] = 0.059912;		gaussianMask[1][4] = 0.015019;
	gaussianMask[2][0] = 0.023792;		gaussianMask[2][1] = 0.094907;		gaussianMask[2][2] = 0.150342;		gaussianMask[2][3] = 0.094907;		gaussianMask[2][4] = 0.023792;
	gaussianMask[3][0] = 0.015019;		gaussianMask[3][1] = 0.015019;		gaussianMask[3][2] = 0.094907;		gaussianMask[3][3] = 0.059912;		gaussianMask[3][4] = 0.015019;
	gaussianMask[4][0] = 0.003765;		gaussianMask[4][1] = 0.015019;		gaussianMask[4][2] = 0.023792;		gaussianMask[4][3] = 0.015019;		gaussianMask[4][4] = 0.003765;

	// Define Sobel Masks
	GxMask[0][0] = -1; GxMask[0][1] = 0; GxMask[0][2] = 1;
	GxMask[1][0] = -2; GxMask[1][1] = 0; GxMask[1][2] = 2;
	GxMask[2][0] = -1; GxMask[2][1] = 0; GxMask[2][2] = 1;

	GyMask[0][0] = 1; GyMask[0][1] = 2; GyMask[0][2] = 1;
	GyMask[1][0] = 0; GyMask[1][1] = 0; GyMask[1][2] = 0;
	GyMask[2][0] = -1; GyMask[2][1] = -2; GyMask[2][2] = -1;



	/* convert our image to a 2-D array, grayscale. */
	i = 0; // treat i as a row counter.
	j = 0; // treat j as a column counter.
	for (k = 0; k < height_rows*width_cols; k++){
		greyscale_image[i][j] = bitmap[k];
		j++;
		if (j == width_cols)
		{
			j = 0;
			i++;
		}
	}
	//m_free(bitmap);


	// Create a copy of the image.
	for (i = 0; i <height_rows; i++)
	{
		for (j = 0; j < width_cols; j++)
			//temp_copy[i][j] = 0;
			gradients[i][j] = 0;
	}

	i = 0;
	j = 0;



	/* Do Gaussian filter, not counting the two pixel borders. We can do this because we can assume that these pixels should be black anyways. */

	for (i = 2; i < height_rows - 2; i++)
	{
		for (j = 2; j < width_cols - 2; j++)
		{
			temp_value_filter = 0;

			for (row_counter = -2; row_counter < 3; row_counter++)
			{
				for (col_counter = -2; col_counter < 3; col_counter++)
				{
					temp_value_filter = temp_value_filter + gaussianMask[2 + row_counter][2 + col_counter] * greyscale_image[i + row_counter][j + col_counter];
				}
			}

			//temp_copy[i][j] = (int)temp_value_filter;
			gradients[i][j] = (int)temp_value_filter;
		}

	}
	// we're doing this to get rid of temp_copy and to therefore save array space

	for (i = 0; i <height_rows; i++)
	{
		for (j = 0; j < width_cols; j++){
			greyscale_image[i][j] = gradients[i][j];
		}
	}

	for (i = 0; i <height_rows; i++)
	{
		for (j = 0; j < width_cols; j++){
			//temp_copy[i][j] = 0;
			gradients[i][j] = 0;
		}
	}


	/* END gaussian filter */

	/*Do Sobel filter ON OUR TEMP_COPY. Because it's about edges, we'll be more careful than in Gaussian to make sure that the borders are included.*/

	for (i = 1; i < height_rows - 1; i++)
	{
		for (j = 1; j < width_cols - 1; j++)
		{

			gx_temp = 0;
			gy_temp = 0;

			gx_temp = gx_temp + GxMask[0][0] * greyscale_image[i - 1][j - 1];
			gy_temp = gy_temp + GyMask[0][0] * greyscale_image[i - 1][j - 1];

			gx_temp = gx_temp + GxMask[1][0] * greyscale_image[i][j - 1];
			gy_temp = gy_temp + GyMask[1][0] * greyscale_image[i][j - 1];

			gx_temp = gx_temp + GxMask[2][0] * greyscale_image[i + 1][j - 1];
			gy_temp = gy_temp + GyMask[2][0] * greyscale_image[i + 1][j - 1];

			gx_temp = gx_temp + GxMask[0][1] * greyscale_image[i - 1][j];
			gy_temp = gy_temp + GyMask[0][1] * greyscale_image[i - 1][j];

			gx_temp = gx_temp + GxMask[2][1] * greyscale_image[i + 1][j];
			gy_temp = gy_temp + GyMask[2][1] * greyscale_image[i + 1][j];

			gx_temp = gx_temp + GxMask[0][2] * greyscale_image[i - 1][j + 1];
			gy_temp = gy_temp + GyMask[0][2] * greyscale_image[i - 1][j + 1];

			gx_temp = gx_temp + GxMask[1][2] * greyscale_image[i][j + 1];
			gy_temp = gy_temp + GyMask[1][2] * greyscale_image[i][j + 1];

			gx_temp = gx_temp + GxMask[2][2] * greyscale_image[i + 1][j + 1];
			gy_temp = gy_temp + GyMask[2][2] * greyscale_image[i + 1][j + 1];

			gradients[i][j] = sqrt(gx_temp*gx_temp + gy_temp*gy_temp);

			// For debugging purposes
			if (gradients[i][j] > temp_var_max_gradient) {
				temp_var_max_gradient = gradients[i][j];
			}

			thisAngle = (atan2(gx_temp, gy_temp) / 3.14159) * 180.0;		// Calculate actual direction of edge

																							/* Convert actual edge direction to approximate value */
			if (((thisAngle < 22.5) && (thisAngle > -22.5)) || (thisAngle > 157.5) || (thisAngle < -157.5))
				newAngle = 0;
			if (((thisAngle > 22.5) && (thisAngle < 67.5)) || ((thisAngle < -112.5) && (thisAngle > -157.5)))
				newAngle = 45;
			if (((thisAngle > 67.5) && (thisAngle < 112.5)) || ((thisAngle < -67.5) && (thisAngle > -112.5)))
				newAngle = 90;
			if (((thisAngle > 112.5) && (thisAngle < 157.5)) || ((thisAngle < -22.5) && (thisAngle > -67.5)))
				newAngle = 135;

			edge_dir[i][j] = newAngle;	// Store the approximate edge direction of each pixel in one array
		}

	}

	//for debugging purposes
	printf("The maximum gradient value encountered is: %d\n", temp_var_max_gradient);

	/*END sobel filter */


	// Create a 1-pixel wide border on the edges of gradients. This is because we never actually calculate the gradients at the edges.
	for (i = 0; i < height_rows; i++)
	{
		gradients[i][0] = 0;
		gradients[i][width_cols - 1] = 0;
		edge_dir[i][0] = 0;
		edge_dir[i][width_cols - 1] = 0;
	}
	for (j = 0; j < height_rows; j++)
	{
		gradients[0][j] = 0;
		gradients[height_rows - 1][j] = 0;
		edge_dir[0][j] = 0;
		edge_dir[height_rows - 1][j] = 0;
	}


	/*  Begin edge filtering.*/
	// we will re-fill our greyscale_image matrix with the canny-edge-detection results.

	for (i = 0; i < height_rows; i++)
	{
		for (j = 0; j < width_cols; j++)
		{
			greyscale_image[i][j] = 0;
		}
	}

	 // Non-maximal surpression
	    // we're going to surpress all edges that are non-maximal.

	    for (i = 1; i < height_rows - 1; i++)
	    {
	        for (j = 1; j < width_cols - 1; j++)
	        {
	            if (edge_dir[i][j] == 0) { // e.g. horizontal change --> vertical line
	                                        // For each edge, there can be two edge-pixels adjacent to it in a line.
	                                        // next_edge_direction[0] will be the next pixel, next_edge_direction[1] will be the previous pixel
	                                        // next_edge_direction[0][0] is the x-value increment of the next pixel, and [0][1] is the y-value increment
	                                        //
	                                        // Note that change in x-value means the COLUMN vector is updated (e.g. j), whereas a change in y means ROW (i) is updated!
	                next_edge_direction[0][0] = 0;
	                next_edge_direction[0][1] = 1;
	                next_edge_direction[1][0] = 0;
	                next_edge_direction[1][1] = -1;
	            }
	            if (edge_dir[i][j] == 45) { // e.g. the CHANGE is 45 degrees, so the line should be 135 degrees
	                next_edge_direction[0][0] = -1;
	                next_edge_direction[0][1] = 1;
	                next_edge_direction[1][0] = 1;
	                next_edge_direction[1][1] = -1;
	            }
	            if (edge_dir[i][j] == 90) { // e.g. the change is verticle, so we want a horizontal line.
	                next_edge_direction[0][0] = 1;
	                next_edge_direction[0][1] = 0;
	                next_edge_direction[1][0] = -1;
	                next_edge_direction[1][1] = 0;
	            }
	            if (edge_dir[i][j] == 135) { // e.g. the change is at 135 degrees, so we want a 45 degree line
	                next_edge_direction[0][0] = -1;
	                next_edge_direction[0][1] = -1;
	                next_edge_direction[1][0] = 1;
	                next_edge_direction[1][1] = 1;
	            }

	            if (gradients[i + next_edge_direction[0][1]][j + next_edge_direction[0][0]] > gradients[i][j])
	            {
	                gradients[i][j] = 0;
	            }
	            if (gradients[i + next_edge_direction[1][1]][j + next_edge_direction[1][0]] > gradients[i][j])
	            {
	                gradients[i][j] = 0;
	            }

	        }
	    }


	    // we're now going to populate our final edge-only picture
	    // DESIGN CHOICE:
	    // DESIGN CHOICE
	    // Basically, for hysterisis to run properly, you have to run it multiple times. Because something that was previously in a lower threshold could now
	    // be a part of the edge, where before its threshold was considered to be too low. Because now it's a part of the edge, you increase its gradient value
	    // to make it above the threshold. This means you have to then check the pixel adjacent to THAT pixel too. You can do this two ways: (1) via some sort
	    // of external function which is called as needed (you must be careful to make it not run infinitely) [hardest to do but probably the best] or just have
	    // a while loop which runs the above stuff like 5 times in a row. Which means your edge will be cut off if you have 6 stuff in a row below the high threshold but
	    // will be okay if you have e.g. 1 above low threshold, 2 above low threshold, 1 above high threshold, 4 above low threshold. Actually, try the while
	    // loop implementation first. It's easiest to impliment and more flexible for now
	    //
	    for (k = 0; k < num_hysterisis_loops; k++)
	    {
	        for (i = 1; i < height_rows - 1; i++)
	        {
	            for (j = 1; j < width_cols - 1; j++)
	            {

	                if (gradients[i][j] > upperThreshold)
	                {
	                    greyscale_image[i][j] = binarize_upper_limit;
	                    // we're going to binarize the image - anything that's an edge will be given one value and anything that's not will be given another.

	                    if (edge_dir[i][j] == 0) { // e.g. horizontal change -> implies vertical line.

	                        // For each edge, there can be two edge-pixels adjacent to it in a line.
	                        // next_edge_direction[0] will be the next pixel, next_edge_direction[1] will be the previous pixel
	                        // next_edge_direction[0][0] is the x-value increment of the next pixel, and [0][1] is the y-value increment
	                        //
	                        // Note that change in x-value means the COLUMN vector is updated (e.g. j), whereas a change in y means ROW (i) is updated!
	                        next_edge_direction[0][0] = 1;
	                        next_edge_direction[0][1] = 0;
	                        next_edge_direction[1][0] = -1;
	                        next_edge_direction[1][1] = 0;
	                    }
	                    if (edge_dir[i][j] == 45) {
	                        next_edge_direction[0][0] = 1;
	                        next_edge_direction[0][1] = 1;
	                        next_edge_direction[1][0] = -1;
	                        next_edge_direction[1][1] = -1;
	                    }
	                    if (edge_dir[i][j] == 90) {
	                        next_edge_direction[0][0] = 0;
	                        next_edge_direction[0][1] = 1;
	                        next_edge_direction[1][0] = 0;
	                        next_edge_direction[1][1] = -1;
	                    }
	                    if (edge_dir[i][j] == 135) {
	                        next_edge_direction[0][0] = -1;
	                        next_edge_direction[0][1] = 1;
	                        next_edge_direction[1][0] = 1;
	                        next_edge_direction[1][1] = -1;
	                    }

	                    // time for us to extend our edge by hysterisis.

	                    //let's analyze the next pixel
	                    if (gradients[i + next_edge_direction[0][1]][j + next_edge_direction[0][0]] > lowerThreshold)
	                    {
	                        // this way, in our next pass of our hysterisis function, this index will be recognized as well.
	                        gradients[i + next_edge_direction[0][1]][j + next_edge_direction[0][0]] = upperThreshold + 1;
	                        greyscale_image[i + next_edge_direction[0][1]][j + next_edge_direction[0][0]] = binarize_upper_limit;
	                    }
	                    if (gradients[i + next_edge_direction[1][1]][j + next_edge_direction[1][0]] > lowerThreshold)
	                    {
	                        // this way, in our next pass of our hysterisis function, this index will be recognized as well.
	                        gradients[i + next_edge_direction[1][1]][j + next_edge_direction[1][0]] = upperThreshold + 1;
	                        greyscale_image[i + next_edge_direction[1][1]][j + next_edge_direction[1][0]] = binarize_upper_limit;
	                    }

	                } // end if (gradients[i][j] > upperThreshold)
	            } // end for (j = 1; j < width_cols - 1; j++)
	        } // end for (i = 1; i < height_rows - 1; i++)
	    } // end for (k = 0; k < num_hysteresis_loops; k++)

	  // Create a 5-pixel wide border on the edges of gradients. This is because sometimes the images have weird border edges because of pre-processing.
	for (i = 0; i < height_rows; i++)
	{
		greyscale_image[i][0] = 0;
		greyscale_image[i][width_cols - 1] = 0;
		greyscale_image[i][1] = 0;
		greyscale_image[i][width_cols - 2] = 0;
		greyscale_image[i][2] = 0;
		greyscale_image[i][width_cols - 3] = 0;

	}
	for (j = 0; j < height_rows; j++)
	{
		greyscale_image[0][j] = 0;
		greyscale_image[height_rows - 1][j] = 0;
		greyscale_image[1][j] = 0;
		greyscale_image[height_rows - 2][j] = 0;
		greyscale_image[2][j] = 0;
		greyscale_image[height_rows - 3][j] = 0;
	}

	printf("We out here\n");

}

void read_file_visual_studio()
{

	// TO DO: CREATE the smaller_test_image file in MATLAb
	if (demonstration_number == 1)
		ptr_file = fopen("C:\\Program Files (x86)\\Texas Instruments\\c6sdk_02_00_00_00\\demos\\Fracture_detect\\Documents\\ulnarfracture150x200.txt", "r");
	else
		ptr_file = fopen("C:\\Program Files (x86)\\Texas Instruments\\c6sdk_02_00_00_00\\demos\\Fracture_detect\\Documents\\ulnarNOfracture150x250.txt", "r");

	// TO DO
	// TO DO
	if (!ptr_file)
	{
		printf("Error opening file!\n");
		return;
	}
	i = 0;
	while (1)
	{
		if (fgets(buf, 10, ptr_file) != NULL) {
			bitmap[i] = atoi(buf);
			i++;
		}
		else
			break;
	}
	fclose(ptr_file);
	//printf("i = %d\n", i);
	//printf("bitmap[0] should be 0, and is = %d\nbitmap[315344] should be 129, and is = %d\nbitmap[42000] should be 0, and is = %d\n", bitmap[0], bitmap[315344], bitmap[41999]);

}

int main()
{

	if (demonstration_number == 1) // case with fracture.
	{
		upperThreshold = 85;
		lowerThreshold = 50;
	}
	else {
		upperThreshold = 120;
		lowerThreshold = 30;
	}

	//mem_init();

	//
	// NOTE: there are two methods of reading in your image:
	// (1)	save as a text file, line by line in MATLAB
	//		use getline to somehow create a variable line by line
	// (2)  copy/paste the raw data to C.

	//bitmap = imread("C:\\Users\\EE113D\\MyProjects\\Fergilicious_SH\\Mini Project 1\\Test5.bmp");

	//temp_image = canny(bitmap);
	//m_free(bitmap);

	//bitmap = temp_image;


	//////////////////////////////////
	///// READ FILE FOR VISUAL STUDIO
	//////////////////////////////////
	printf("begin reading \n");
	read_file_visual_studio();
	//// END READ FILE FOR VISUAL STUDIO
	printf("end reading \n");

	printf("started canny \n");
	canny(bitmap);
	printf("Finished canny \n");


	ss_image=cvCreateImage(cvSize(width_cols, height_rows), IPL_DEPTH_8U, 1 );// create the image structure for your image.

	row_counter = 0;
	col_counter = 0;
	for (counter = 0; counter < ss_image->widthStep*height_rows; counter++)
	{
		if (col_counter  < width_cols){
			ss_image->imageData[counter] = greyscale_image[row_counter][col_counter];
		}
		else{
			ss_image->imageData[counter] = 0;
		}

		col_counter++;

		if (col_counter == ss_image->widthStep)
		{
			col_counter = 0;
			row_counter++;
		}
	}

	printf("ss_image->imageSize: %d \n, ss_image->nChannels: %d \n, ss_image->depth: %d \n, ss_image->widthStep: %d \n ss_image->width: %d \n, ss_image->height: %d \n", ss_image->imageSize, ss_image->nChannels, ss_image->depth, ss_image->widthStep, ss_image->width, ss_image->height );


	/////////////////////////////////////
	///// BEGIN HOUGH TRANSFORM ANALYSIS
	/////////////////////////////////////
	CvMemStorage* storage = cvCreateMemStorage(0);
	printf("Begin Hough transform-based analysis /n");
	// Previous: lines = cvHoughLines2(ss_image, storage, CV_HOUGH_PROBABILISTIC, 4, CV_PI/180, 20, 20, 15 );
	lines = cvHoughLines2(ss_image, storage, CV_HOUGH_PROBABILISTIC, 4, CV_PI/90, 20, 20, 3 );

	printf("(x,y) of the first coordinate set returned by Probabilistic Hough Transform for each line\n");
	for (i = 0; i < lines->total; i++){
		// for ease of transferring variables to MATLAB
		line1 = (CvPoint*)cvGetSeqElem(lines,i);
		if (line1[0].y != 2 && line1[0].y != 148)
			printf("%d, %d; ", line1[0].x, line1[0].y);
	}
	printf("\n(x,y) of the second coordinate set returned by Probabilistic Hough Transform for each line\n");
	for (i = 0; i < lines->total; i++){
		// for ease of transferring variables to MATLAB
		line1 = (CvPoint*)cvGetSeqElem(lines,i);
		if (line1[0].y != 2 && line1[0].y != 148)
			printf("%d, %d; ", line1[1].x, line1[1].y);
	}
	printf("\n\n");


	// Test code for debugging
	//	line1 = (CvPoint*)cvGetSeqElem(lines,1);
	//	line2 = (CvPoint*)cvGetSeqElem(lines,2);
	//
	//	ss_m1 = (float)((line1[1].y -line1[0].y)/((float)(line1[1].x - line1[0].x)));
	//	ss_b1 = line1[1].y - ss_m1*line1[1].x;
	//	ss_m2 = (float)((line2[1].y -line2[0].y)/((float)(line2[1].x - line2[0].x)));
	//	ss_b2 = line2[1].y - ss_m2*line2[1].x;
	//
	//	printf("m1 %f, b1 %f; m2 %f, b2 %f \n", ss_m1, ss_b1, ss_m2, ss_b2);

	printf("Any potential fractures are located at: \n");
	for( i = 0; i < lines->total; i++ )
	{
		// see if we have two lines that intersect, such that their intersection point is
		// within 10 pixels of the ends of one of the lines.

		line1 = (CvPoint*)cvGetSeqElem(lines,i);
		if (line1[1].y == 148 || line1[1].y == 2)
		{
			continue;
		}

		for (j = i+1; j < lines->total; j++){

			line2 = (CvPoint*)cvGetSeqElem(lines,j);
			// Hough Transform points are in line1[0] and line1[1]
			// and in line2[0] and line1[1];

			if (line2[1].y == 148 || line2[1].y == 2)
				continue;

			ss_m1 = (float)((line1[1].y -line1[0].y)/((float)(line1[1].x - line1[0].x)));
			ss_b1 = line1[1].y - ss_m1*line1[1].x;
			ss_m2 = (float)((line2[1].y -line2[0].y)/((float)(line2[1].x - line2[0].x)));
			ss_b2 = line2[1].y - ss_m2*line2[1].x;


			ss_x_intersection = (ss_b2-ss_b1)/(ss_m1- ss_m2);// this formula is correct, even if it looks weird.
			ss_y_intersection = ss_m1*(ss_x_intersection) + ss_b1;

			// test code for debugging
			//printf("m1 %f, b1 %f; m2 %f, b2 %f \n", ss_m1, ss_b1, ss_m2, ss_b2);
			//printf("m1 %f, b1 %f; m2 %f, b2 %f ", (line1[1].y -line1[0].y)/(line1[1].x - line1[0].x), line1[1].y - ss_m1*line1[1].x, (line2[1].y -line2[0].y)/(line2[1].x - line2[0].x), line2[1].y - ss_m2*line2[1].x);
			//printf("line1[1].y = %d, line1[0] = %d, line1[1].x = %d, line1[0].x = %d, line2[1].y = %d, line2[0] = %d, line2[1].x = %d, line2[0].x = %d; m1 %d, b1 %d; m2 %d, b2 %d ", line1[1].y, line1[0], line1[1].x, line1[0].x, line2[1].y, line2[0], line2[1].x, line2[0].x, ss_m1, ss_b1, ss_m2, ss_b2);
			//printf("intersection: %f, %f; \n", ss_x_intersection, ss_y_intersection);

			if ( ss_x_intersection < width_cols && ss_x_intersection > -1  && ss_y_intersection < height_rows  && ss_y_intersection > 0 ) {

				if (ss_m1 < ss_m2*1.2 && ss_m1 > ss_m2*0.8){
					// do nothing - this is the case where the slopes are so similar we might have two segments of the same bone
				}
				else {
					printf("(%f, %f)\n", ss_x_intersection, ss_y_intersection);
				}

				ss_intersection_distance = (ss_x_intersection - line1[0].x)*(ss_x_intersection - line1[0].x) + (ss_y_intersection - line1[0].y)*(ss_y_intersection - line1[0].y);
				if (ss_intersection_distance < ss_max_intersection_square_dist){
					// this intersection matches!
					printf("%d, %d", ss_x_intersection, ss_y_intersection);
					continue;
				}
				ss_intersection_distance = (ss_x_intersection - line1[1].x)*(ss_x_intersection - line1[1].x) + (ss_y_intersection - line1[1].y)*(ss_y_intersection - line1[1].y);
				if (ss_intersection_distance < ss_max_intersection_square_dist){
					// this intersection matches!
					printf("%d, %d", ss_x_intersection, ss_y_intersection);
					continue;
				}
				ss_intersection_distance = (ss_x_intersection - line2[0].x)*(ss_x_intersection - line2[0].x) + (ss_y_intersection - line2[0].y)*(ss_y_intersection - line2[0].y);
				if (ss_intersection_distance < ss_max_intersection_square_dist){
					// this intersection matches!
					printf("%d, %d", ss_x_intersection, ss_y_intersection);
					continue;
				}
				ss_intersection_distance = (ss_x_intersection - line2[1].x)*(ss_x_intersection - line2[1].x) + (ss_y_intersection - line2[1].y)*(ss_y_intersection - line2[1].y);
				if (ss_intersection_distance < ss_max_intersection_square_dist){
					// this intersection matches!
					printf("%d, %d", ss_x_intersection, ss_y_intersection);
					continue;
				}

			} // end if ( ss_x_intersection < width_cols && ss_x_intersection > -1  && ss_y_intersection < height_rows  && ss_y_intersection > 0 )

		} // end for (j = i+1; j < lines->total; j++)
	} // for( i = 0; i < lines->total; i++ )

	printf("\n \n");
	printf("End Hough transform-analysis\n");
	//// END HOUGH TRANSFORM ANALYSIS




	////////////////////////////////////
	/////// OUTPUT SAVING SEQUENCE
	////////////////////////////////////
	printf("Beginning file save \n");
	// Save output file for analysis
	ptr_file = fopen("C:\\Program Files (x86)\\Texas Instruments\\c6sdk_02_00_00_00\\demos\\Fracture_detect\\Documents\\output_file2.txt", "w");
	if (ptr_file == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	// rewrite ss_image to greyscale_image just to check things.
	// this is because if our saved output file is correct, then this means
	// that we must have saved in ss_image correctly.
	row_counter = 0;
	col_counter = 0;
	for (counter = 0; counter < width_cols*height_rows; counter++) {
		if (col_counter  < width_cols){
			greyscale_image[row_counter][col_counter] = ss_image->imageData[counter];
		}

		col_counter++;

		if (col_counter == ss_image->widthStep)
		{
			col_counter = 0;
			row_counter++;
		}
	}

	for (i = 0; i < height_rows; i++)
	{
		for (j = 0; j <width_cols; j++) {
			fprintf(ptr_file, "%d\n", greyscale_image[i][j]); // you should generally output greyscale_image, not gradients
		}
	}
	fclose(ptr_file);
	printf("File saved and closed\n");
//	//// END OUTPUT SAVE SEQUENCE

	printf("OpenCV image released!\n");
	cvReleaseImage(&ss_image);//release the memory occupied by the image

	//printf("Enter something to close program: ");
	//scanf("%s", str);

	while(1);

}


/***************************** End Of File ************************************/
