# <i>Drosophila</i> wind-assisted dispersal

This repository contains the software, processed data, and all behavioral simulations associated with the paper "The long distance flight behavior of <i>Drosophila</i> supports an agent-based model for wind-assisted dispersal in insects."
For a pre-print, see: https://www.biorxiv.org/content/10.1101/2020.06.10.145169v1

The original data (~ 0.5 TB) will be made available upon reasonable request upon formal publication.

Processed data is available through this repository, along with instructions for use, under the folder "./processed_field_data"

This readme assumes working knowledge of Ubuntu and python. This code is not actively maintained. It was last updated on 09 April 2021.

Code and data are licensed under a [Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0) License](https://creativecommons.org/licenses/by-nc-sa/4.0/ "CC BY-NC-SA 4.0").

Data in this repository are stored in a variety of data formats, all offering reasonable human-readability, and all readable with open-source tools.
* most processed data and experimental metadata are saved as .json files
* some field data, such as anemometer recordings, are saved as .csv files

## Agent-based behavioral simulations
The folder "./agent_based_simulations_models_I_through_IV" houses:
* scripts used to simulate the groundspeed-windspeed relationships expected from fixed-azimuth models (I through IV):
  * simulations with free parameters set to values estimated from the literature (as in Fig. 6), in "/simulations"
  * sensitivity analyses to individually optimize each model (as in Supp. Figs 2 and 3) on the basis of their fit to field data, in "/simulations_sensitivity_analysis"
  * finally, comparing these individually-optimized models on the basis of relative fit to field data (as in Supp. Fig 4), in "/bayes_comparing_optimized_simulations"
* scripts used to simulate groundspeed-windspeed relationships from model I, above, with a suite of stochastic variations introduced ("/model_I_with_stochastic_heading_changes")
* an outdated folder ("/models_V_and_VI_NOT_PRESENTED_IN_PNAS_PAPER") containing preliminary attempts at stochastic simulations that we removed from the paper during the PNAS review process (not removed from repository due to low competence with git).

## Advection-diffusion behavioral simulations
All scripts used to run advection-diffusion simulations were written by Will Dickson, and are packaged in a separate GitHub repository at: https://github.com/willdickson

The folder "./advection_diffusion_simulations" houses:
* final datasets produced from the above advection-diffusion simulations, in "/advection_diffusion_datasets"
* scripts and intermediate figure sequences used to assess the above datasets with respect to their fit to field data (as in Fig. 7), in "/advection_diffusion_analyses"

## The raw camera-trap data
If you would like to re-run our analyses starting from our raw camera-trap image stream, please contact the authors (leitchka@gmail.com). Camera data from each experiment comprise 30 to 120 GB of data; we will make available, upon reasonable request, these data from individual experiments or from the entire dataset (~0.5 TB).

Once you have the camera-trap images, you will need to follow the following instructions for making these data accessible to the analysis (below).

## Key software needed to run our analyses
* Ubuntu (we used Ubuntu 16)
* Python (2.7)
* OpenCV (3.3.1)

## Making raw camera-trap data accessible to our analysis pipeline
We have provided a directory structure that will facilitate your analysis of any raw camera-trap data you've requested and received. This directory is called "./field_data_and_analysis_scripts." Subfolders with experimental date names (e.g. "2017_10_26") have been prepared for your manual deposition of the raw camera data we provide. This will require two steps:

1) First, copy the data into the folder: "./field_data_and_analysis_scripts/2017_10_26/trapcam_timelapse_TO_BE_MANUALLY_POPULATED"
2) Then, rename the folder into which you've pasted the data, so it simply reads "/trapcam_timelapse"

You'll want the resultant directory structure to look like:

		field_data_and_analysis_scripts		
			2017_10_26
				trapcam_timelapse
					trap_A					<----------------  
						mask.jpg					|
						tl_0007_0001_20161130_125959.jpg		|
						.						|
						.						|
						.       [many more image files]			|
						.					these are the folders/files
						.					provided on request
						tl_0007_3644_20161130_160239.jpg		|
					trap_B							|
						mask.jpg					|
						tl_0009_0001_20161130_103056.jpg		|
						.						|
						.						|
						.       [et cetera]		<----------------				



## Running analyses of raw camera-trap data
At this point, using "/field_data_and_analysis_scripts" as your working directory, you'll want to run the script "/run_trapcam_analysis.py," which will prompt your specification of:
* the experiments you'd like to analyze (e.g. 2017_10_26),
* the traps you'd like to analyze (e.g. trap_G),
and
* whether you'd like to tinker around with analysis parameters, or if you'd like to use default analysis parameters to perform the final analysis.
  * if you choose the former, the analysis will come up with an "in-trap/on-trap threshold" for each trap analyzed, on the basis of any dips found in a  multimodal histogram of detected-contour contrast metrics. This threshold will be saved in the parameter file, (e.g. "/2017_10_26/gaussian_analysis_parameters.json").
  * if you choose the latter, the analysis will use the saved "in-trap/on-trap threshold" to generate time-series estimations of the number of flies in, and on, the trap surface; these results will be saved in e.g. "/2017_10_26/all_traps_final_analysis_output.json"
