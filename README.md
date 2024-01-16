# Analysis pipeline for TAH project
Semi-automated analysis of large brightfield IHC images (**IHC_pipeline.m**), with manual curation of ROI selection & artefact removal. Additional IF pipeline (**IF_pipeline.m**) + behaviour results extraction (**BEHAVIOUR_results.m**) included.

# Software requirements
* Matlab (tested with version 2021a)
* Matlab Toolboxes: Image Processing, Statistics and Machine Learning
* Fiji
* [MIJ](https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab) (add to Libraries folder)
* [CopyToImagePlus](https://github.com/kouichi-c-nakamura/copytoImagePlus) package (add to Libraries folder)

# Setup instructions
1) FIJI-Matlab communication  
a) Install Fiji if you don't have it already, e.g. on Desktop ('C:\Users\Pat\Desktop\Fiji.app')  
b) Copy mij.jar file from MIJ folder in Libraries folder into your Matlab java folder, e.g. 'C:\Program Files\MATLAB\R2021a\java\jar\mij.jar'  
c) In utils/setup_miji(): change directory folders to match your configuration. This ensures Matlab can use java to control Fiji
   e.g. Matlab folder from step b: 'C:\Program Files\MATLAB\R2021a\java\jar\mij.jar'
   Fiji installation from step a: 'C:\Users\Pat\Desktop\Fiji.app\jars\ij-1.53q.jar';'C:\Users\Pat\Desktop\Fiji.app\scripts'  

2) Directory  
In IHC_simple: change directory of analysis_folder variable to match location of the analysis folder on your computer
  e.g. analysis_folder = 'C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC';

3) Test whether the scripts work by running IHC_simple.m file after completing the setup steps. Example image has been provided inside Images folder


# Analysis steps
1) Deconvolution of H and DAB channels, done on the full image (Fiji) -> FUNCTION: funcs_pipeline/deconvolve_full
2) Selection of ROIs for each image (right hipp, left hipp, right cortex, left cortex, control area) -> FUNCTION: funcs_pipeline/select_roi
3) Quantification of % area covered by antibody -> FUNCTION: funcs_pipeline/analyse_data


To define the threshold the antibody labelling has to pass to be counted, go to utils/get_antibody_threshold and adjust the values


FILE & FOLDER NAMES ARE IMPORTANT!

Folder structure
* Data is stored under 'xxx\IHC\Images\Subject_name'
* Deconvolved images are saved inside 'xxx\IHC\Images\Subject_name\Processed'
* Selected ROIs are saved inside 'xxx\IHC\ROIs\Subject_name'
* Results are saved inside 'xxx\IHC\Results\Subject_name'

Expected file names: subject_name-antibody_name, e.g. AD-Hipp28_moc23.svs
