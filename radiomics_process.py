"""
radiomics_process.py
=============
Module to run IBSI compliant radiomics feature extraction in mercure.

"""

# Standard Python includes
import os
import sys
import json
import shutil
import subprocess
from pathlib import Path
from rt_utils import RTStructBuilder
import SimpleITK as sitk
import radiomics
from radiomics import featureextractor
import numpy as np
import six
from mirp import extract_features,extract_mask_labels


# Imports for loading DICOMs
import pydicom


def main(args=sys.argv[1:]):
    """
    Main entry function of the . 
    The module is called with two arguments from the function docker-entrypoint.sh:
    'testmodule [input-folder] [output-folder]'. The exact paths of the input-folder 
    and output-folder are provided by mercure via environment variables
    """
    # Print some output, so that it can be seen in the logfile that the module was executed
    print(f"Starting mercure-radiomics")

    # Check if the input and output folders are provided as arguments
    if len(sys.argv) < 3:
        print("Error: Missing arguments!")
        print("Usage: testmodule [input-folder] [output-folder]")
        sys.exit(1)

    # Check if the input and output folders actually exist
    in_folder = sys.argv[1]
    out_folder = sys.argv[2]
    if not Path(in_folder).exists() or not Path(out_folder).exists():
        print("IN/OUT paths do not exist")
        sys.exit(1)

    # # Load the task.json file, which contains the settings for the processing module
    # try:
    #     with open(Path(in_folder) / "task.json", "r") as json_file:
    #         task = json.load(json_file)
    # except Exception:
    #     print("Error: Task file task.json not found")
    #     sys.exit(1)

    # # Create default values for all module settings
    # settings = {"trained_model": "S2_osem_b10_fdg_pe2i", "series_desc": "Bowsher_" , "series_suffix":"DEFAULT"}

    # # Overwrite default values with settings from the task file (if present)
    # if task.get("process", ""):
    #     settings.update(task["process"].get("settings", {}))

    
    # filter image and mask
    current_dir = os.getcwd()
    mask_path = os.path.join(current_dir, 'mask')
    if not os.path.exists(mask_path):
        os.makedirs(mask_path)

    image_path = os.path.join(current_dir, 'image')
    if not os.path.exists(image_path):
        os.makedirs(image_path)


    #read modality and move to relevant input directories
    series = []
    for entry in os.scandir(in_folder):
        if entry.name.endswith(".dcm") and not entry.is_dir():
            
            target_path = ''
            
            dcm_file_in = Path(in_folder) /  entry.name
            ds = pydicom.dcmread(dcm_file_in)
            
            series_number = ds.SeriesInstanceUID
            if series_number not in series:
                series.append(series_number)
                if len(series) >2:
                    print("Error: More than two series in input directory.")
                    sys.exit(1)

            modality = ds.Modality
            if modality=='MR':
                target_path =image_path
                
            if modality=='RTSTRUCT':
                target_path=mask_path
                
            if (target_path):
                shutil.move(os.path.join(in_folder, entry.name), target_path)
            else:
                print("Error: Error copying files.")
                sys.exit(1)
    
    # selected_model=settings["trained_model"]
    # series_desc=settings["series_desc"]
    # series_suffix=settings["series_suffix"]
    # Load existing RT Struct. Requires the series path and existing RT Struct path
    
    mask_file_in = os.path.join(mask_path, os.listdir(mask_path)[0])
    print(mask_file_in)
    rtstruct = RTStructBuilder.create_from(
        dicom_series_path=image_path, 
        rt_struct_path=mask_file_in
        )

    # View all of the ROI names from within the image
    roi_list = rtstruct.get_roi_names()
    print(roi_list)
    selected_roi = roi_list[0]

    # Loading the 3D Mask from within the RT Struct ## need to sort out datatype!!!!!
    mask_volume = np.array(rtstruct.get_roi_mask_by_name(selected_roi))
    mask_array = sitk.GetImageFromArray(mask_volume.astype(np.uint8))
    
    
    dcm_images = [image.pixel_array for image in rtstruct.series_data ]
    image_volume = np.stack(dcm_images, axis=2)

    image_array = sitk.GetImageFromArray(image_volume.astype(np.float32))

    extractor = featureextractor.RadiomicsFeatureExtractor()
    print('Extraction parameters:\n\t', extractor.settings)
    print('Enabled filters:\n\t', extractor.enabledImagetypes)
    print('Enabled features:\n\t', extractor.enabledFeatures)	
    
    result = extractor.execute(image_array, mask_array)
    print('Result type:', type(result))  # result is returned in a Python ordered dictionary)
    print('')
    print('Calculated features')
    for key, value in six.iteritems(result):
        print('\t', key, ':', value)


    mirp_feature_data = extract_features(
        image=image_path,
        mask=mask_path,
        roi_name=[selected_roi],
        base_discretisation_method="fixed_bin_number",
        base_discretisation_n_bins=32
    )

    print(mirp_feature_data)

if __name__ == "__main__":
    main()