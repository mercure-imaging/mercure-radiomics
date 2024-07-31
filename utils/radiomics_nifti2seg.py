import os
import sys
import numpy as np
import nibabel as nib
import pydicom
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import generate_uid, ImplicitVRLittleEndian
from datetime import datetime
import argparse
import glob
from rt_utils import RTStructBuilder
import random
import string
from pathlib import Path

import glob
import highdicom as hd
from pydicom.sr.codedict import codes

def create_dicom_file(pixel_array, output_path, study_instance_uid, series_instance_uid, instance_number,ref_uid):
    # Create a new FileDataset instance
    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4'  # MR Image Storage
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
    
    
    
    ds = FileDataset(output_path, {}, file_meta=file_meta, preamble=b"\0" * 128)
    
    file_name = str(os.path.basename(output_path))
    patient_ID = file_name.split('_', 1)
    output_path_lower = output_path.lower()

    if 'flair' in output_path_lower:
        s_desc = "flair"
        s_num = 1
    elif 't2' in output_path_lower:
        s_desc = "t2"
        s_num = 2
    elif '_t1_' in output_path_lower:
        s_desc = "t1"
        s_num = 3
    elif '_t1ce_' in output_path_lower:
        s_desc= "t1ce"
        s_num = 4


    # Add a fake study date (current date in this example)
    fake_date = datetime.now().strftime("%Y%m%d")
    ds.StudyDate = '20001220'
    ds.PatientBirthDate='20001220'
    ds.PatientSex='M'
    ds.AccessionNumber='000000'
    # Optionally, you can also add a fake study time
    fake_time = datetime.now().strftime("%H%M%S")
    ds.StudyTime = '121212'

    # Generate a random Study ID (e.g., 8 characters long)
    fake_study_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))

    ds.SeriesDescription = s_desc
    # Add the fake Study ID to the Dataset
    ds.StudyID = '10000000'
    # Add the data elements
    ds.PatientName = 'ANONYMOUS'
    ds.PatientID = patient_ID #+ datetime.now().strftime('%Y%m%d')
    ds.StudyInstanceUID = study_instance_uid
    ds.SeriesInstanceUID = series_instance_uid
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.SOPClassUID = file_meta.MediaStorageSOPClassUID
    ds.Modality = "MR"
    ds.SeriesNumber = s_num
    #ds.InstanceNumber = instance_number
    ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    # tag which slice in the volume we're dealing with
    ds.ImagePositionPatient = [0, 0, instance_number]
    ds.FrameOfReferenceUID = ref_uid

    ds.BitsStored = 16
    ds.BitsAllocated = 16
    ds.SamplesPerPixel = 1
    ds.HighBit = 15

    ds.Rows, ds.Columns = pixel_array.shape
    ds.PixelSpacing = [1.0, 1.0]
    ds.SliceThickness = 1.0

    ds.PixelRepresentation = 0
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.PixelData = pixel_array.tobytes()

    # Separate directory and filename
    directory, filename = os.path.split(output_path)

    # Add subdirectory '1'
    new_directory = os.path.join(directory, str(s_num))
    if not os.path.exists(new_directory):
        os.makedirs(new_directory)
    # Put the full path back together
    new_output_path = os.path.join(new_directory, filename)
     
    ds.save_as(new_output_path)
    return output_path

def nifti_to_dicom(input_dir, output_dir, segmentation_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Find all .nii and .nii.gz files in the input directory
    nifti_files = glob.glob(os.path.join(input_dir, '*.nii')) + glob.glob(os.path.join(input_dir, '*.nii.gz'))

    if not nifti_files:
        print(f"No NIfTI files found in {input_dir}")
        return

    # Generate study information
    study_instance_uid = generate_uid()

    # Dictionary to store segmentation data
    segmentations = {}

    # Load segmentation files
    segmentation_files = glob.glob(os.path.join(segmentation_dir, '*.nii')) + glob.glob(os.path.join(segmentation_dir, '*.nii.gz'))
    for seg_file in segmentation_files:
        seg_img = nib.load(seg_file)
        seg_data = seg_img.get_fdata()
        seg_data = np.rot90(seg_data, 1, (0, 1))
        print('seg shape=', seg_data.shape)
        segmentations[os.path.basename(seg_file)] = seg_data

    dicom_series = []  # To store all DICOM files for the series

    for idx, nifti_file in enumerate(nifti_files):
        # Extract filename without extension
        base_filename = os.path.splitext(os.path.basename(nifti_file))[0]
        if base_filename.endswith('.nii'):  # Handle .nii.gz case
            base_filename = os.path.splitext(base_filename)[0]

        # Load NIfTI file
        nifti_img = nib.load(nifti_file)
        nifti_data = nifti_img.get_fdata()
        # rotate nii to match DICOM orientation
        nifti_data = np.rot90(nifti_data, 1, (0, 1))  # rotate segmentation in-plane
        # Generate series information
        series_instance_uid = generate_uid()
        f_uid = generate_uid()
        print('nifti shape=', nifti_data.shape)
        # Iterate through slices
        for slice_idx in range(nifti_data.shape[2]):
            # Extract slice data
            slice_data = nifti_data[:, :, slice_idx]

            # Normalize and convert to uint16
            slice_data = ((slice_data - slice_data.min()) / (slice_data.max() - slice_data.min())) * 65535
            slice_data = slice_data.astype(np.uint16)

            # Create and save DICOM file
            output_filename = f'{base_filename}_slice{slice_idx:03d}.dcm'
            output_path = os.path.join(output_dir, output_filename)
            dicom_path = create_dicom_file(slice_data, output_path, study_instance_uid, series_instance_uid, slice_idx + 1, f_uid)
            dicom_series.append(dicom_path)

        print(f'Converted {nifti_file} to DICOM series')

    print(f"Converted {len(nifti_files)} NIfTI file(s) to DICOM format.")
    print(dicom_series[0])
    # Create RTSTRUCT file
    if dicom_series:
        
        # Get all items in the directory
        items = os.listdir(output_dir)
        # Filter for only directories
        subdirectories = [item for item in items if os.path.isdir(os.path.join(output_dir, item))]
        # Initialize RTStructBuilder with the first DICOM file in the series
        #path_without_filename = os.path.dirname(full_path)
        for sub_dir in subdirectories:
            out_series=os.path.join(output_dir, sub_dir)
            #rtstruct_filename = os.path.join(out_series, f'{base_filename}_RTSTRUCT.dcm')
            #rtstruct = RTStructBuilder.create_new(out_series)

            # Add segments to RTSTRUCT
            # for seg_name, seg_data in segmentations.items():
            #     for region_value in np.unique(seg_data):
            #         if region_value == 0:  # Skip background
            #             continue
            #         mask = seg_data == region_value
            #         rtstruct.add_roi(
            #             mask=mask,
            #             name=f"Segmentation_{int(region_value)}"
            #         )
            
            # # Save RTSTRUCT file
            # rtstruct.save(rtstruct_filename)
            # print(f"Generated RTSTRUCT file: {rtstruct_filename}")

            # Describe the algorithm that created the segmentation
            algorithm_identification = hd.AlgorithmIdentificationSequence(
                name='nifti_to_seg',
                version='v1.0',
                family=codes.DCM.ArtificialIntelligence
            )
            print('got to loop')
            out_series = Path(out_series)
            image_files = sorted(out_series.glob("*.dcm"), reverse=True)
            #print(image_files)
            source_images = [pydicom.dcmread(str(f)) for f in image_files]
            print('got past dcm_read')

            
            series_number_out = source_images[0].SeriesNumber
            series_desc_out = source_images[0].SeriesDescription
            tumor_masks = []
            for seg_name, seg_data in segmentations.items():
                for region_value in np.unique(seg_data):
                    if region_value == 0:  # Skip background
                        continue
                    #ask=np.rollaxis(np.array(seg_data == region_value),2)
                    mask = seg_data == region_value
                    mask=np.flip(mask, axis=2)
                    tumor_masks.append(mask)
                    print('region value found:', region_value)
                    
            # Create segment descriptions for each tumor region
            segment_descriptions = [
                create_tumor_segment_description(i+1, f"Tumor {i+1}")
                for i in range(len(tumor_masks))
            ]
            seg_series = 1000 + series_number_out
            seg_series_description =  'SEG_' + series_desc_out
            # Stack all tumor masks into a single 4D array
            pixel_array = np.stack(tumor_masks, axis=3)
            pixel_array=np.moveaxis(pixel_array,2,0)
            print(pixel_array.shape)
            # Create the Segmentation object
            seg = hd.seg.Segmentation(
                source_images=source_images,
                pixel_array=pixel_array,
                #segmentation_type=hd.seg.DimensionOrganizationTypeValues.BY_SEGMENT,
                segmentation_type=hd.seg.SegmentationTypeValues.BINARY,
                segment_descriptions=segment_descriptions,
                series_instance_uid=generate_uid(),
                series_number=seg_series,
                series_description=seg_series_description,
                sop_instance_uid=generate_uid(),
                instance_number=1,
                manufacturer="Your Manufacturer",
                device_serial_number='Mercure',
                manufacturer_model_name="Your Model",
                software_versions="1.0"
            )

            # Save the DICOM SEG file
            output_path = os.path.join(out_series, "seg_.dcm")
            seg.save_as(output_path)

            print(f"DICOM SEG file with {len(tumor_masks)} tumor regions saved to {output_path}")
    


def create_tumor_segment_description(segment_number, segment_label):
    algorithm_identification = hd.AlgorithmIdentificationSequence(
            name='mercure-radiomics',
            version='v1.0',
            family=codes.DCM.ArtificialIntelligence
        )
    return hd.seg.SegmentDescription(
            segment_number=segment_number,
            segment_label=segment_label,
            segmented_property_category=codes.SCT.Organ,
            segmented_property_type=codes.SCT.Organ,
            # segmented_property_category_code_sequence=[{
            #     "CodeValue": "M-80003",
            #     "CodingSchemeDesignator": "SCT",
            #     "CodeMeaning": "Neoplasm, Primary"
            # }],
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.MANUAL,
            algorithm_identification=algorithm_identification,
            tracking_uid=hd.UID(),
            tracking_id='RadiomicsROI'
            #algorithm_name="Manual Tumor Segmentation"
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert NIfTI files to DICOM format and generate SEG files.')
    parser.add_argument('input_dir', help='Input directory containing NIfTI files')
    parser.add_argument('output_dir', help='Output directory for DICOM files')
    parser.add_argument('segmentation_dir', help='Directory containing NIfTI segmentation files')
    args = parser.parse_args()

    nifti_to_dicom(args.input_dir, args.output_dir, args.segmentation_dir)