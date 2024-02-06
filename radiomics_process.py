"""
radiomics_process.py
=============
Module to run IBSI compliant radiomics feature extraction in mercure.

"""

# Standard Python includes
import os
import sys
import json
import stat
import shutil
import glob
import csv
import subprocess
from pathlib import Path
from rt_utils import RTStructBuilder
import SimpleITK as sitk
import radiomics
from radiomics import featureextractor
import numpy as np
import pandas as pd
import six
from mirp import extract_features,extract_mask_labels
import highdicom as hd
from pydicom.sr.codedict import codes
from rt_utils import ds_helper, image_helper

from pydicom.filereader import dcmread
from pydicom.sr.codedict import codes
from pydicom.uid import generate_uid
from highdicom.sr.content import FindingSite
from highdicom.sr.templates import Measurement, TrackingIdentifier


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
    try:
        with open(Path(in_folder) / "task.json", "r") as json_file:
            task = json.load(json_file)
    except Exception:
        print("Error: Task file task.json not found")
        sys.exit(1)

    # # Create default values for all module settings
    settings = {"rois": [1], "processor": "pyradiomics" , "processor_settings":"default"}
    
    # # Overwrite default values with settings from the task file (if present)
    if task.get("process", ""):
         settings.update(task["process"].get("settings", {}))

    
    # filter image and mask
    current_dir = os.getcwd()
    mask_path = os.path.join(current_dir, 'mask')
    if not os.path.exists(mask_path):
        os.makedirs(mask_path)

    image_path = os.path.join(current_dir, 'image')
    if not os.path.exists(image_path):
        os.makedirs(image_path)

    rt_struct_output_path = os.path.join(current_dir, 'rt_struct_output')
    if not os.path.exists(rt_struct_output_path):
        os.makedirs(rt_struct_output_path)

    settings_path = os.path.join(current_dir, 'settings')
    if not os.path.exists(settings_path):
        os.makedirs(settings_path) 

    #read modality and move to relevant input directories
    modality_list =[]
    series = []
    first_image_path=''
    image_ds_list=[]
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
            if modality=='MR' or  modality=='CT':
                target_path =image_path
                if modality not in modality_list: modality_list.append(modality)
                if first_image_path=='':first_image_path=dcm_file_in
                image_ds_list.append(ds)
                
            if modality=='RTSTRUCT':
                target_path=mask_path
                if modality not in modality_list: modality_list.append(modality)

            if modality=='SEG':
                target_path=mask_path
                if modality not in modality_list: modality_list.append(modality)
                #image_ds_list.append(ds)
            if (target_path):
                shutil.copy(os.path.join(in_folder, entry.name), target_path)
            else:
                print("Error: Error copying files.")
                sys.exit(1)
    
    selected_processor=settings["processor"]
    processor_settings=settings["processor_settings"]
    roi_settings=settings["rois"]
    parameter_json=settings["processing_parameters"]
    # series_suffix=settings["series_suffix"]
    # Load existing RT Struct. Requires the series path and existing RT Struct path
    
    mask_file_in = os.path.join(mask_path, os.listdir(mask_path)[0])
    print(mask_file_in)

    if 'RTSTRUCT' in modality_list: 
        rtstruct = RTStructBuilder.create_from(
            dicom_series_path=image_path, 
            rt_struct_path=mask_file_in
            )

        # View all of the ROI names from within the image
        roi_list = rtstruct.get_roi_names()
        print(roi_list)
        if roi_settings[0] in roi_list:
            selected_roi = roi_settings[0]
            print('ROI found, extracting features for ', selected_roi)
        else:
            print('ROI not found:', roi_settings[0])

        
        # Loading the 3D Mask from within the RT Struct ## need to sort out datatype!!!!!
        mask_volume = np.array(rtstruct.get_roi_mask_by_name(selected_roi))
        mask_array = sitk.GetImageFromArray(mask_volume.astype(np.uint8))

        dcm_images = [image.pixel_array for image in rtstruct.series_data ]
        image_volume = np.stack(dcm_images, axis=2)
    elif 'SEG' in modality_list:
        seg = hd.seg.segread(mask_file_in)
        # segment_numbers = seg.get_segment_numbers(
        #       segmented_property_type=codes.SCT.Organ
        #    )
        segment_numbers=roi_settings
        print(segment_numbers)
        print(seg.number_of_segments )
        if  seg.number_of_segments > 0 :
            source_image_uids = []
            for study_uid, series_uid, sop_uid in seg.get_source_image_uids():
                print(study_uid, series_uid, sop_uid)
                source_image_uids.append(sop_uid)
            # Retrieve a binary segmentation mask for these images for the bone segment
            mask_volume = seg.get_pixels_by_source_instance(
                    source_sop_instance_uids=source_image_uids,
                    segment_numbers=segment_numbers,
            )
            mask_volume =np.squeeze(mask_volume)
            mask_array = sitk.GetImageFromArray(mask_volume.astype(np.uint8))
            image_series_data=[]
            for root, _, files in os.walk(image_path):
                for file in files:
                    try:
                        ds = pydicom.dcmread(os.path.join(root, file))
                        if hasattr(ds, "pixel_array"):
                            image_series_data.append(ds)
                    except Exception:
                        # Not a valid DICOM file
                        continue
            dcm_images = [image.pixel_array for image in image_series_data ]
            print (len(dcm_images))
            image_volume = np.stack(dcm_images, axis=0)
        else:
            print('error')

    image_array = sitk.GetImageFromArray(image_volume.astype(np.float32))

    if selected_processor=='pyradiomics':
        if processor_settings=='default':
            extractor = featureextractor.RadiomicsFeatureExtractor()
        elif processor_settings=='parameters':
            if Path(settings_path).exists():
                settings_file = os.path.join(settings_path,"settings.json")
            with open(settings_file, "w") as write_file:
                json.dump(parameter_json, write_file, indent=4)
            p = Path(settings_file)
            p.chmod(p.stat().st_mode | stat.S_IROTH | stat.S_IXOTH | stat.S_IWOTH)
            
            extractor = featureextractor.RadiomicsFeatureExtractor()
            extractor.loadParams(settings_file)

        
        print('Extraction parameters:\n\t', extractor.settings)
        print('Enabled filters:\n\t', extractor.enabledImagetypes)
        print('Enabled features:\n\t', extractor.enabledFeatures)	
        
        result = extractor.execute(image_array, mask_array)
        print('Result type:', type(result))  # result is returned in a Python ordered dictionary)
        print('')
        print('Calculated features')
        
        #display results and convert arrays to results for json serialization
        r_dict = {}
        for key, value in six.iteritems(result):
            print('\t', key, ':', value)
            if isinstance(value, np.ndarray):
                value = value.tolist()
            r_dict[key] =value
        
        json_results= write_result_json(out_folder, r_dict)

    elif selected_processor=='mirp':
        #MIRP currently working with RTSTRUCT - could convert SEG to RTSTRUCT?
        mirp_mask_path=mask_path
        if 'SEG' in modality_list:
            rt_struct_output_filename = 'output-rt-struct_vols.dcm'
            #create new RT Struct - requires original DICOM
            rtstruct = RTStructBuilder.create_new(dicom_series_path=image_path)
            rtstruct.add_roi(
                    mask=np.moveaxis((mask_volume>0),0,2),
                    name='selected_region'
                )
            rtstruct.save(os.path.join(rt_struct_output_path, rt_struct_output_filename))
            mirp_mask_path=rt_struct_output_path
            selected_roi='selected_region'
        if processor_settings=='default':
            mirp_feature_data = extract_features(
                image=image_path,
                mask=mirp_mask_path,
                roi_name=[selected_roi],
                base_discretisation_method="fixed_bin_number",
                base_discretisation_n_bins=32
            )
        elif processor_settings=='parameters':
            args = {
                    "image":image_path,
                    "mask":mirp_mask_path,
                    "roi_name":[selected_roi]
                    }
            args.update(parameter_json)
            print(args)
            mirp_feature_data = extract_features(**args)

        print(type(mirp_feature_data[0]))   
        #convert to dictionary and print
        mirp_dict = mirp_feature_data[0].loc[0].to_dict()
        #print(mirp_dict)
        for key, value in six.iteritems(mirp_dict):
            print('\t', key, ':', value)
        json_results = write_result_json(out_folder, mirp_dict)

    
    # generate SR - reshape array if RT struct   
    if 'RTSTRUCT' in modality_list:
        mask_volume=np.flip(np.moveaxis(mask_volume,2,0),0)
    seg_sr_writer(mask_volume, image_path, out_folder, json_results)
    
    #clean up
    if os.path.exists(mask_path):
        shutil.rmtree(mask_path)
        print('mask path deleted.')
    
    if os.path.exists(image_path):
        shutil.rmtree(image_path)
        print('image path deleted.')
    
    if os.path.exists(rt_struct_output_path):
        shutil.rmtree(rt_struct_output_path)
        print('rtstruct path deleted.')
    
    if os.path.exists(settings_path):
        shutil.rmtree(settings_path)
        print('settings path deleted.')


def write_result_json(output_dir, results_dict):
    if Path(output_dir).exists():
        results_file = os.path.join(output_dir,"results.json")
        with open(results_file, "w") as write_file:
            json_string = json.dumps(results_dict, indent=4)
            write_file.write(json_string)
        p = Path(results_file)
        p.chmod(p.stat().st_mode | stat.S_IROTH | stat.S_IXOTH | stat.S_IWOTH)

        # Write the dictionary to a CSV file
        results_csv = os.path.join(output_dir,'radiomics_results.csv')
        with open(results_csv, 'w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=results_dict.keys())
            writer.writeheader()
            writer.writerow(results_dict)
    
    return json_string

def seg_sr_writer(mask, series_dir, write_dir, json_results):
        series_dir = Path(series_dir)

        image_files = sorted(series_dir.glob("*.dcm"), reverse=True)

        # Read CT Image data sets from PS3.10 files on disk
        image_datasets = [pydicom.dcmread(str(f)) for f in image_files]
        print('got past dcm_read')
        # Read CT Image data sets from PS3.10 files on disk
        
        
        # Describe the algorithm that created the segmentation
        algorithm_identification = hd.AlgorithmIdentificationSequence(
            name='mercure-radiomics',
            version='v1.0',
            family=codes.DCM.ArtificialIntelligence
        )
        
        seg_name='radiomics'
        seg_series_desc='mercure-radiomics_ROI_seg'
        sct_type = codes.SCT.Organ
        # Describe the segment
        description_segment_1 = hd.seg.SegmentDescription(
            segment_number=1,
            segment_label=seg_name,
            segmented_property_category=codes.SCT.Organ,
            segmented_property_type=sct_type,
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=algorithm_identification,
            tracking_uid=hd.UID(),
            tracking_id='RadiomicsROI'
        )

        # Create the Segmentation instance
        
        print(np.shape(mask), len(image_datasets))
        seg_dataset = hd.seg.Segmentation(
            source_images=image_datasets,
            pixel_array=mask,
            segmentation_type=hd.seg.SegmentationTypeValues.BINARY,
            segment_descriptions=[description_segment_1],
            series_instance_uid=hd.UID(),
            series_number=1001,
            sop_instance_uid=hd.UID(),
            instance_number=1,
            manufacturer='CBI',
            manufacturer_model_name='Mercure',
            software_versions='v1',
            device_serial_number='Mercure',
            series_description=seg_series_desc,
        )
        
        #print(seg_dataset)
        seg_file_path = os.path.join(write_dir, "seg_"+seg_name+".dcm")
        seg_dataset.save_as(seg_file_path)


        #write structured report:

        #TODO write IBSI measurment definitions and loop:
        # A measurement using an IBSI code (not in pydicom)
        # histogram_intensity_code = hd.sr.CodedConcept(
        #     value="X6K6",
        #     meaning="Intensity Histogram Mean",
        #     scheme_designator="IBSI",
        # )
        # hist_measurement = hd.sr.Measurement(
        #     name=histogram_intensity_code,
        #     value=-119.0738525390625,
        #     unit=codes.UCUM.HounsfieldUnit,
        # )

        # A segmentation dataset, assumed to contain a segmentation of the source
        # image above
        seg = dcmread(seg_file_path)

        # Information about the observer
        observer_person_context = hd.sr.ObserverContext(
            observer_type=codes.DCM.Person,
            observer_identifying_attributes=hd.sr.PersonObserverIdentifyingAttributes(
                name='Doe^John'
            )
        )
        observer_device_context = hd.sr.ObserverContext(
            observer_type=codes.DCM.Device,
            observer_identifying_attributes=hd.sr.DeviceObserverIdentifyingAttributes(
                uid=hd.UID(),
            )
        )


        observation_context = hd.sr.ObservationContext(
            observer_person_context=observer_person_context,
            observer_device_context=observer_device_context,
        )

        
        
        # A tracking identifier for this measurement group
        tracking_id = hd.sr.TrackingIdentifier(
        identifier='Region3D0001',
        uid=hd.UID(),
        )

        # # Define the image region using a specific segment from the segmentation
        ref_segment = hd.sr.ReferencedSegment.from_segmentation(
        segmentation=seg,
        segment_number=1,
        )

        # # Construct the measurement group
        group = hd.sr.VolumetricROIMeasurementsAndQualitativeEvaluations(
        referenced_segment=ref_segment,
        tracking_identifier=tracking_id,
        #measurements=[...],
        #qualitative_evaluations=[...],
        )
        
        radiomics_text_concept = hd.sr.CodedConcept(
             value="X99X",
             meaning="Radiomics results output",
             scheme_designator="IBSI",
        )

        radiomics_results_text = hd.sr.TextContentItem(
            name=radiomics_text_concept,
            value=json_results,
            relationship_type=hd.sr.RelationshipTypeValues.HAS_OBS_CONTEXT,
        )
        
        observation_context.append(radiomics_results_text)

        measurement_report = hd.sr.MeasurementReport(
          observation_context=observation_context,  # from above
          procedure_reported=codes.LN.CTUnspecifiedBodyRegion,
          imaging_measurements=[group],
          title=codes.DCM.ImagingMeasurementReport,
        )

        sr_series_desc='mercure-radiomics_ROI_SR'
        image_datasets.append(seg)
        # Create the Structured Report instance
        sr_dataset = hd.sr.Comprehensive3DSR(
            evidence=image_datasets,  # all datasets referenced in the report
            content=measurement_report,
            series_number=2001,
            series_instance_uid=hd.UID(),
            sop_instance_uid=hd.UID(),
            instance_number=1,
            series_description=sr_series_desc,
            manufacturer='Manufacturer'
        )
        

        #add radiomics from result.json
        sr_file_path = os.path.join(write_dir, "sr_"+seg_name+".dcm")
        sr_dataset.save_as(sr_file_path)


if __name__ == "__main__":
    main()