# **mercure-radiomics**
<br>

Mercure module to run [IBSI](https://theibsi.github.io/) compliant radiomics feature extraction.

Current radiomics processing packages included in the module are [Pyradiomics](https://github.com/AIM-Harvard/pyradiomics) and [MIRP](https://github.com/oncoray/mirp). This module runs as a docker container in mercure, it can be added to an existing mercure installation using docker tag : *mercureimaging/mercure-radiomics*. The module will perform radiomic feature extraction for DICOM image series sent to mercure with accompanying ROIs in DICOM RTSTRUCT or SEG format.

Radiomic feature measurments are output by the module in .csv format, .json results that can be retrieved from the mercure database, and DICOM Structured Reports for onward transfer and analysis.
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

By default, all ROIs in the segmentation file will be processed using the pyradiomics feature extractor. Settings can be configured in the mercure user interface to select feature processor, ROIs to be processed, and processor settings. These are shown with descriptions below: 
<pre>
"rois": ["Tumor"]                   <i>default = "ALL" to process all ROIs found.</i>
"processor": "pyradiomics"          <i>default = "pyradiomics", can also be set to "mirp".</i>
"processor_settings":"parameters"   <i>default = "default" to set to default processor settings, if 'parameters' selected 'processing_parameters' can be provided.</i>
"processing_parameters":{}          <i>default is empty, see examples below to set processor parameters.</i>
</pre>
<br>
Example of MIRP processor selection with 'TV' ROI and specified processing parameters:
<pre>
{
    "rois":["TV"],
    "processor": "mirp",
    "processor_settings":"parameters",
    "processing_parameters":{
                            "base_discretisation_method":"fixed_bin_number",
                            "base_discretisation_n_bins":32
                            }
}
</pre>
<br>
Example of Pyradiomics processor selection with 'TV' ROI and specified processing parameters:
<pre>
{
"rois":["TV"],
"processor": "pyradiomics",
"processor_settings":"parameters",
"processing_parameters":{
                        "setting": {
                            "label": 1,
                            "interpolator": "sitkBSpline",
                            "resampledPixelSpacing": null
                        },
                        "featureClass": {
                            "shape": [
                                "VoxelVolume",
                                "MeshVolume",
                                "SurfaceArea",
                                "SurfaceVolumeRatio",
                                "Compactness1",
                                "Compactness2",
                                "Sphericity",
                                "SphericalDisproportion",
                                "Maximum3DDiameter",
                                "Maximum2DDiameterSlice",
                                "Maximum2DDiameterColumn",
                                "Maximum2DDiameterRow",
                                "MajorAxisLength",
                                "MinorAxisLength",
                                "LeastAxisLength",
                                "Elongation",
                                "Flatness"
                            ]
                        }
                        }
}
</pre>
<br>
<br>

## **Acknowledgments**

<br>
Pyradiomics: https://github.com/AIM-Harvard/pyradiomics

Medical Image Radiomics Processor (MIRP): https://github.com/oncoray/mirp

<br>



