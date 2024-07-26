# **mercure-radiomics**
<br>

Mercure module to run [IBSI](https://theibsi.github.io/) compliant radiomics feature extraction.

Current radiomics processing packages included in the module are [Pyradiomics](https://github.com/AIM-Harvard/pyradiomics) and [MIRP](https://github.com/oncoray/mirp).This module runs as a docker container in mercure, it can be added to an existing mercure installation using docker tag : *mercureimaging/mercure-radiomics*. The module will perform radiomic feature extraction for DICOM image series sent to mercure with accompanying ROIs in DICOM RTSTRUCT or SEG format.

<br>

## **Installation**

### Add module to existing mercure installation
Follow instructions on [mercure website](https://mercure-imaging.org) on how to add a new module. Use the docker tag mercureimaging/mercure-radiomics.

<br>

### Build module for local testing, modification and development
1. Clone repo.
2. Build Docker container locally by running make (modify makefile with new docker tag as needed).
3. Test container :\
`docker run -it -v /input_data:/input -v /output_data:/output --env MERCURE_IN_DIR=/input  --env MERCURE_OUT_DIR=/output mercureimaging/mercure-radiomics`

<br>

## **Configuration**

The mercure-radiomics module requires little configuration in mercure. A rule should be configured to select a single image series and accompanying RTSTRUCT or SEG series. More information on mercure rule configuration can be found [here.](https://mercure-imaging.org/docs/usage.html)

By default, all ROIs in the segmentation file will be processed using the pyradiomics feature extractor. Settings can be configured in the mercure user interface to select fearture processor, ROIs to be processed, and processor settings. Some examples are provided below: 
```
"gaussian_filter": "false" (TBC) 
```

<br>
<br>

## **Acknowledgments**

<br>
Pyradiomics: https://github.com/AIM-Harvard/pyradiomics

Medical Image Radiomics Processor (MIRP): https://github.com/oncoray/mirp

<br>



