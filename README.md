Collection of scripts to process fluorescence imaging data.
The main scripts for fluorescent transient detection and analysis are glia_data_processing and cellbody_data_processing, with the former being an automated event detection solution and the later requiring the input of the event timings.
Simply provide the folder path where the data is located and an event timing file as a csv, for the directed pipeline, to the pipelines.
Figures can be generated using the figure_generator script from the pipeline outputs.

The components of these pipelines are included as seperate files (deltafs, glia_training and csv cumulator). In addition there are utility scripts for running batch processing with suite2p, deleting ROIs from suite2p outputs, preprocessing of two-photon recordings and multi-recording alignment.
