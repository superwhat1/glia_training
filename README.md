Collection of scripts to process fluorescence imaging data.
The main scripts for fluorescent transient detection and analysis are multi_roi_auto_data_processing, for automated transient detection on recordings with multiple regions of interest (ROIs), multi_roi_stim_data_processing, for processing multiple ROIs from recordings with timed stimuli, and matching scripts for single ROI recordings, single_roi_stim_data_processing and single_roi_auto_data_processing.
Simply provide the folder path where the data is located and an event timing file as a csv, for the stimuli pipelines, in the appropriate script.
Figures can be generated using the figure_generator script from the pipeline outputs.

The components of these pipelines are included as seperate files (deltafs, glia_training and csv cumulator). In addition there are utility scripts for running batch processing with suite2p, deleting ROIs from suite2p outputs, preprocessing of two-photon recordings and multi-recording alignment.
