# Wind-assisted dispersal in <i>Drosophila</i>

This repository contains the software, processed data, and all behavioral simulations associated with the paper "The long distance flight behavior of <i>Drosophila</i> suggests a general model for wind-assisted dispersal in insects."
For a pre-print, see: https://www.biorxiv.org/content/10.1101/2020.06.10.145169v1

The original data (~ 0.5 TB) will be available upon reasonable request upon formal publication.

Processed data is available through this repository, along with instructions for use, under the folders "./processed_field_data"

This readme assumes working knowledge of Ubuntu and python. This code is not actively maintained. It was last updated on 29 June 2020, using up-to-date versions of the required software below.

Code and data are licensed under a [Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0) License](https://creativecommons.org/licenses/by-nc-sa/4.0/ "CC BY-NC-SA 4.0").


## What you need to run our analysis
* Ubuntu (we used Ubuntu 16)
* Python (2.7)
* OpenCV (3.3.1)
* apt-get repositories: 
* Manual downloads: 




1) DATA FORMATS
This repository contains a mixture of different data formats, all offering reasonable human-readability, and all able to be read into a python environment using only open-source tools. 
* most processed data and experimental metadata are saved as .json files 
* some field data, such as anemometer recordings, are saved as .csv files
* 

## The data

If you would like to re-run our analyses starting from our raw camera-trap image stream, please contact the authors (leitchka@gmail.com). Camera data from each experiment comprise 30 to 120 GB of data; we will make available, upon reasonable request, these data from individual experiments or from the entire dataset (~0.5 TB).

Once you have the camera-trap images, you will need to follow the following instructions for making these data accessible to the analysis (below).
 
## Making the data accessible to analysis
We have provided a directory structure that will facilitate your analysis of any raw camera-trap data you've requested and received. This directory is called "./field_data_and_analysis_scripts." Subfolders with experimental date names (e.g. "2017_10_26") have been prepared for your manual deposition of the raw camera data we provide. This will require two steps:
1) First, copy the data into the folder: "./field_data_and_analysis_scripts/2017_10_26/trapcam_timelapse_TO_BE_MANUALLY_POPULATED"
2) Then, rename the folder into which you've pasted the data, so it simply reads "/trapcam_timelapse"

You'll want the resultant directory structure to look like:
		
field_data_and_analysis_scripts		
	2017_10_26
		trapcam_timelapse
			trap_A					<--------------------------------  
				mask.jpg							|
				tl_0007_0001_20161130_125959.jpg				|
				.								|
				.								|
				.       [many more image files]					|
				.						these are the folders/files 
				.							provided on request
				tl_0007_3644_20161130_160239.jpg				|
			trap_B									|
				mask.jpg							|
				tl_0009_0001_20161130_103056.jpg				|
				.								|
				.								|
				.       [et cetera]		<--------------------------------				



## Running the 





