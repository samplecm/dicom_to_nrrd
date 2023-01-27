import pydicom 
import radiomics
import os
import pyplastimatch
import six
import glob
import numpy as np
import dicom_structure_finding
from get_contour_masks import get_all_structure_masks
import nrrd
import argparse
import re

class Img_Series:
    def __init__(self, data_dir, modality):
        #object will store the imgs as one 4d array, where 4th dimension is for storing cartesian corrdinates of every pixel. 
        metadata_list = get_img_series_metadata(data_dir)
        iop = metadata_list[0].ImageOrientationPatient
        self.name = str(metadata_list[0].PatientName)
        self.acquisition_date = str(metadata_list[0].AcquisitionDate)
        self.modality = modality
      # we will refer to patient cartesian coordinates as x,y,z and the image coordinates as a,b,c
        self.cos_ax = iop[0]
        self.cos_ay = iop[1]
        self.cos_az = iop[2]
        self.cos_bx = iop[3]
        self.cos_by = iop[4]
        self.cos_bz = iop[5]

        #get cosines with c by taking cross product. 
        cross_prod = np.cross(np.array(iop[0:3]), np.array(iop[3:]))
        self.cos_cx = cross_prod[0]
        self.cos_cy = cross_prod[1]
        self.cos_cz = cross_prod[2]

        self.slice_thickness = metadata_list[0].SliceThickness
        self.pixel_spacing = metadata_list[0].PixelSpacing
        self.origin = metadata_list[0].ImagePositionPatient

        x_orig = self.origin[0]
        y_orig = self.origin[1]
        z_orig = self.origin[2]
        img_list = []
        ipp_list = []
        for data in metadata_list:
            img_list.append(data.pixel_array)
            ipp_list.append(data.ImagePositionPatient)
        row_count, col_count = img_list[0].shape    
        #will return a shape [3 , (#z) , (#y), (#x)] array where first index has coordinate specifications.
        
        #1 --> x value of each pixel
        #2 --> y value of each pixel 
        #3 --> z value of each pixel 
        #array is sorted from smallest to largest z. 
        img_array = np.zeros([len(img_list), row_count, col_count])
        coords_array = np.zeros([3,len(img_list), row_count, col_count])
        coords_array_img = np.zeros([3,len(img_list), row_count, col_count])
        for i,img in enumerate(img_list):
            img_array[i,:,:] = img

        for z_idx in range(len(img_list)):
            corner_x, corner_y, corner_z = ipp_list[z_idx]
            
            for y_idx in range(row_count):
                b= y_idx * self.pixel_spacing[0]
                for x_idx in range(col_count):
                    a = x_idx * self.pixel_spacing[1] 
                    c = z_idx * self.slice_thickness
                    z = corner_z + (b * self.cos_bz) + (a * self.cos_az)
                    x = corner_x + (b * self.cos_bx) + (a * self.cos_ax)
                    y = corner_y + (b * self.cos_by) + (a * self.cos_ay)
                    # test_a = (x-x_orig) * self.cos_ax + (y-y_orig) * self.cos_ay+ (z-z_orig) * self.cos_az
                    # test_b = (x-x_orig) * self.cos_bx + (y-y_orig) * self.cos_by+ (z-z_orig) * self.cos_bz
                    # test_c = (x-x_orig) * self.cos_cx + (y-y_orig) * self.cos_cy+ (z-z_orig) * self.cos_cz
                    coords_array[0,z_idx,y_idx, x_idx] = x
                    coords_array[1,z_idx,y_idx, x_idx] = y
                    coords_array[2,z_idx,y_idx, x_idx] = z

                    coords_array_img[0,z_idx,y_idx, x_idx] = a
                    coords_array_img[1,z_idx,y_idx, x_idx] = b
                    coords_array_img[2,z_idx,y_idx, x_idx] = c
        self.coords_array = coords_array
        self.coords_array_img = coords_array_img    
        self.image_array = img_array          
        print("")


def test():
    img_dir = "/media/sf_U_DRIVE/Profile/Desktop/Programs/IVIM_Code/Analysis_Data/SIVIM02/pre/radiomics_data/images/ivim/dicom"
    struct_dir = "/media/sf_U_DRIVE/Profile/Desktop/Programs/IVIM_Code/Analysis_Data/SIVIM02/pre/radiomics_data/structs/ivim/dicom"
    save_img_dir = "/media/sf_U_DRIVE/Profile/Desktop/Programs/IVIM_Code/Analysis_Data/SIVIM02/pre/radiomics_data/images/ivim/stk"
    save_struct_dir = "/media/sf_U_DRIVE/Profile/Desktop/Programs/IVIM_Code/Analysis_Data/SIVIM02/pre/radiomics_data/structs/ivim/stk"
    # img_dir = r"U:\Profile\Desktop\Programs\IVIM_Code\Analysis_Data\SIVIM02\pre\radiomics_data\images\ivim\dicom"
    # struct_dir = r"U:\Profile\Desktop\Programs\IVIM_Code\Analysis_Data\SIVIM02\pre\radiomics_data\structs\ivim\dicom"
    # save_img_dir = r"U:\Profile\Desktop\Programs\IVIM_Code\Analysis_Data\SIVIM02\pre\radiomics_data\images\ivim\stk"
    # save_struct_dir = r"U:\Profile\Desktop\Programs\IVIM_Code\Analysis_Data\SIVIM02\pre\radiomics_data\structs\ivim\stk"
    # img_series = Img_Series(img_dir)
    # structures_dict = dicom_structure_finding.get_contours(struct_dir)
    # structures_dict_img_coords = convert_contours_to_img_coords(img_series, structures_dict)
    
    # structures_dict_img_coords_on_plane = get_contours_on_img_planes(img_series, structures_dict_img_coords)

    # structures_masks_dict = get_all_structure_masks(img_series, structures_dict_img_coords_on_plane)
    # save_all_img_and_mask_as_nrrd(img_series, structures_masks_dict, save_paths=[save_img_dir, save_struct_dir])
    convert_all_dicoms_to_nrrd(img_dir, struct_dir, "CT", save_paths=[save_img_dir,save_struct_dir])        

def get_img_series_metadata(img_dir, modality="ct"):
    #will return a list which has the metadata for each image in slice order. 
    #list is sorted from smallest to largest z. 

    file_paths = glob.glob(os.path.join(img_dir, "*"))
    img_list = []
    for path in file_paths:
        try:
            metadata = pydicom.dcmread(path)
        except:
            print(f"Could not load {path}. Continuing...")  
            continue  
        if  modality.lower() not in metadata.Modality.lower():
            continue
        img_list.append(metadata)
    img_list.sort(key=lambda x: x.ImagePositionPatient[2])
    return img_list

def get_img_series_array(img_dir):
    
    #    
    meta_list = get_img_series_metadata(img_dir)
    
def convert_contours_to_img_coords(img_series, structures_dict):
    origin = img_series.origin
    x_orig = origin[0]
    y_orig = origin[1]
    z_orig = origin[2]
    cos_ax = img_series.cos_ax 
    cos_ay = img_series.cos_ay 
    cos_az = img_series.cos_az 
    cos_bx = img_series.cos_bx 
    cos_by = img_series.cos_by
    cos_bz = img_series.cos_bz
    cos_cx = img_series.cos_cx
    cos_cy = img_series.cos_cy
    cos_cz = img_series.cos_cz
    new_dict = {}
    for structure_name in structures_dict:
        new_contours = []
        contours = structures_dict[structure_name]
        for contour in contours:
            new_contour = []
            for point in contour:
                #first get distance of x,y,z from image origin
                x_rel = point[0] - x_orig
                y_rel = point[1] - y_orig
                z_rel = point[2] - z_orig

                #now get converted points (x,y,z) --> (a,b,c)
                a = (x_rel * cos_ax) + (y_rel * cos_ay) + (z_rel * cos_az)
                b = (x_rel * cos_bx) + (y_rel * cos_by) + (z_rel * cos_bz)
                c = (x_rel * cos_cx) + (y_rel * cos_cy) + (z_rel * cos_cz)

                new_contour.append([a,b,c])
            new_contours.append(new_contour)   
        new_dict[structure_name] = new_contours 
    return new_dict  

   

def get_contours_on_img_planes(img_series, structures_dict):
    #this function will convert the contours lists such that all contours align with image planes in img_series. 

    #first get the img coord "c" value of each img.
    c_vals = []
    coords_array_img = img_series.coords_array_img
    slice_thickness = float(img_series.slice_thickness)
    for c in coords_array_img[2,:,0,0]:
        c_vals.append(c)
    c_vals = np.array(c_vals)
    new_dict = {}

    for structure_name in structures_dict:
        new_contours = []    
        contours = structures_dict[structure_name]
        for contour in contours:
            if len(contour) == 0:
                continue
            c = contour[0][2]
            closest_indices = np.argsort(np.abs(c_vals-c))
            closest_img_c = c_vals[closest_indices[0]]     #find two closest image slices c vals
            second_closest_img_c = c_vals[closest_indices[1]]

            if np.abs(closest_img_c-c) > slice_thickness:
                print(f"Couldn't find image slice in range for contour at c = {c}")
                continue
            elif np.abs(closest_img_c-c) <= 0.5:
                new_contour = []
                for point in contour:
                    new_point = [point[0], point[1], closest_img_c] 
                    new_contour.append(new_point)
                new_contours.append(new_contour)    
            elif np.abs(closest_img_c-c) < slice_thickness:
                print(f"Warning: contour for {structure_name} at c = {c} assigned to image slice more than 0.5mm away, but less than slice thickness.")
                new_contour = []
                for point in contour:
                    new_point = [point[0], point[1], closest_img_c] 
                    new_contour.append(new_point)
                new_contours.append(new_contour)  
                    
        new_dict[structure_name] = new_contours    
    return new_dict

def save_all_img_and_mask_as_nrrd(img_series, structures_masks_dict, save_paths=os.getcwd()):
    for structure_name in structures_masks_dict:
        save_img_and_mask_as_nrrd(img_series, structure_name, structures_masks_dict[structure_name], save_paths)
    return
def save_img_and_mask_as_nrrd(img_series, structure_name, structure_masks, save_paths):
    img_file_name = img_series.name + "_" + img_series.acquisition_date + ".nrrd"
    struct_file_name = img_series.name + "_" + img_series.acquisition_date + "_" + structure_name + ".nrrd"
    if type(save_paths) is list and len(save_paths) > 1:
        img_save_path = os.path.join(save_paths[0], img_file_name)
        struct_save_path = os.path.join(save_paths[1], struct_file_name)
    elif type(save_paths) is list and len(save_paths) == 1:
        img_save_path = os.path.join(save_paths[0], img_file_name)
        struct_save_path = os.path.join(save_paths[0], struct_file_name)
    elif type(save_paths) == str:
        img_save_path = os.path.join(save_paths, img_file_name)
        struct_save_path = os.path.join(save_paths, struct_file_name)    
    if save_paths is None:  
        img_save_path = os.path.join(os.getcwd(), img_file_name) 
        struct_save_path = os.path.join(os.getcwd(), struct_file_name)   
    
    
    #first need to swap the rows and columns of image and masks because nrrd wants first dimension to be width, second to be height . also sort from largest z to smallest instead of smallest to largest
    img = np.swapaxes(img_series.image_array, 0,2)
    mask = np.swapaxes(structure_masks,0,2)
    #make the header
    header = {'kinds': ['domain', 'domain', 'domain'], 'units': ['mm','mm', 'mm'], 'spacings': [float(img_series.pixel_spacing[1]), float(img_series.pixel_spacing[0]), float(img_series.slice_thickness)]} #'space directions': np.array([[1,0,0], [0,1,0],[0,0,1]])
    try:
        nrrd.write(img_save_path, img, header)    
        print(f"Wrote nrrd file {img_file_name} to {img_save_path}")
    except:
        print(f"Failed to write nrrd file {img_file_name} to {img_save_path}")
         
    try:
        nrrd.write(struct_save_path, mask, header)
        print(f"Wrote nrrd file {struct_file_name} to {struct_save_path}")
    except:    
        print(f"Failed to write nrrd file {struct_file_name} to {struct_save_path}")       
    return

def convert_all_dicoms_to_nrrd(img_dir, struct_dir, modality, save_paths=None):   
    img_series = Img_Series(img_dir, modality)
    structures_dict = dicom_structure_finding.get_contours(struct_dir)
    structures_dict_img_coords = convert_contours_to_img_coords(img_series, structures_dict)
    
    structures_dict_img_coords_on_plane = get_contours_on_img_planes(img_series, structures_dict_img_coords)

    structures_masks_dict = get_all_structure_masks(img_series, structures_dict_img_coords_on_plane)
    save_all_img_and_mask_as_nrrd(img_series, structures_masks_dict, save_paths=save_paths)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="converter for dicom images and structures to make corresponding nrrd files for use with pyradiomics.")
    parser.add_argument("--img_dir", help="specify the folder containing a single patient's image series.", default=None, type=str)
    parser.add_argument("--structure_dir", help="specify the folder containing a single patient's RTSTRUCT file.", default=None, type=str)
    parser.add_argument("--modality", help="specify image modality you wish to convert (Default is CT. Other options are PET or MR)", default="None", type=str)
    parser.add_argument("--save_dir", help="Specify where you wish to save your nrrd files. if You wish to place converted structures in a separate folder, also specify the argument \"save_dir_struct\"", default=None, type=str)
    parser.add_argument("--save_dir_struct", help="Specify a separate folder you'd like to save the struct nrrd files to.", default=None, type=str)
    args = parser.parse_args()
    v = vars(args)
    n_args = sum([1 for a in v.values() if a])
    print("Starting DICOM to NRRD converter.")
    print(f"Supplied image directory: {args.img_dir}")
    print(f"Supplied structure directory: {args.structure_dir}")
    print(f"Supplied output directory: {args.save_dir}")
    print(f"Supplied modality: {args.modality}")

    img_dir = None
    structure_dir = None
    modality = "CT"
    save_dir = None
    if os.path.exists(args.img_dir):
        img_dir = args.img_dir
    if os.path.exists(args.structure_dir):
        structure_dir = args.structure_dir
    if os.path.exists(args.save_dir):
        save_dir = args.save_dir
    if args.img_dir == None:
        while True:
            try:
                input_val = input("\n Please specify the directory containing dicom data for a single patient image series. If structure(s) are in a different folder, please separate directories with a comma.\n>> ")
                input_val = list(input_val.split(","))

                if len(input_val) == 1:    #same directory for images/structures
                    if os.path.exists(input_val[0].strip()):
                        print(f"Loading all data from {input_val}")
                        img_dir =input_val[0].strip()
                        structure_dir = input_val[0].strip()
                        break
                if len(input_val) == 2:    #separate dir for images/structures
                    input_val[0] = input_val[0].strip()    #remove any white space in path
                    input_val[1] = input_val[1].strip()    #remove any white space in path
                    if os.path.exists(input_val[0]) and os.path.exists(input_val[1]):
                        print(f"Loading image data from {input_val[0]}")
                        print(f"Loading structure data from {input_val[1]}")
                        img_dir = input_val[0]
                        structure_dir = input_val[1]
                        break
                if len(input_val) >= 3:
                   print("More than two paths specified. Please include up to 2 separate paths (one for images, one for structures)") 

            except KeyboardInterrupt:
                quit()
            except: pass    
    elif not os.path.exists(args.img_dir):
        while True:
            try:
                input_val = input(f"\n Could not find {args.img_dir}. Please specify a new directory containing image data. If structure(s) are in a different folder, please separate directories with a comma.\n>> ")
                input_val = list(input_val.split(","))

                if len(input_val) == 1:    #same directory for images/structures
                    if os.path.exists(input_val[0].strip()):
                        print(f"Loading all data from {input_val}")
                        img_dir =input_val[0].strip()
                        structure_dir = input_val[0].strip()
                        break
                if len(input_val) == 2:    #separate dir for images/structures
                    input_val[0] = input_val[0].strip()    #remove any white space in path
                    input_val[1] = input_val[1].strip()    #remove any white space in path
                    if os.path.exists(input_val[0]) and os.path.exists(input_val[1]):
                        print(f"Loading image data from {input_val[0]}")
                        print(f"Loading structure data from {input_val[1]}")
                        img_dir = input_val[0]
                        structure_dir = input_val[1]
                        break
                if len(input_val) >= 3 or len(input_val) == 0:
                   print("0 or more than two paths specified. Please include up to 2 separate paths (one for images, one for structures)\n>> ") 

            except KeyboardInterrupt:
                quit()
            except: pass  

    if args.structure_dir == None:
        structure_dir = img_dir
        
    if args.save_dir == None:
        while True:
            try:
                input_val = input("\n Please specify the directory to save nrrd images to. If structure(s) are to be saved in a different folder, please separate directories with a comma.\n>> ")
                input_val = list(input_val.split(","))
                if len(input_val) == 1:
                    if os.path.exists(input_val[0].strip()):
                        save_dir = input_val[0].strip()
                        print(f"saving all nrrd data into {input_val}")
                        break
                if len(input_val) == 2:
                    input_val[0] = input_val[0].strip()    #remove any white space in path
                    input_val[1] = input_val[1].strip()    #remove any white space in path
                    if os.path.exists(input_val[0]) and os.path.exists(input_val[1]):
                        print(f"saving image nrrd data into {input_val[0]}")
                        print(f"saving structure nrrd data into {input_val[1]}")
                        save_dir = input_val
                        break
                if len(input_val) >= 3:
                   print("More than two paths specified. Please include up to 2 separate paths (one for images, one for structures)\n>> ") 

            except KeyboardInterrupt:
                quit()
            except: pass   
    elif not os.path.exists(args.save_dir):
        while True:
            try:
                input_val = input(f"\n Could not find {args.save_dir}. Please specify the directory to save nrrd images to. If structure(s) are to be saved in a different folder, please separate directories with a comma.\n>> ")
                input_val = list(input_val.split(","))
                if len(input_val) == 1:
                    if os.path.exists(input_val[0].strip()):
                        save_dir = input_val[0].strip()
                        print(f"saving all nrrd data into {input_val}")
                        break
                if len(input_val) == 2:
                    input_val[0] = input_val[0].strip()    #remove any white space in path
                    input_val[1] = input_val[1].strip()    #remove any white space in path
                    if os.path.exists(input_val[0]) and os.path.exists(input_val[1]):
                        print(f"saving image nrrd data into {input_val[0]}")
                        print(f"saving structure nrrd data into {input_val[1]}")
                        save_dir = input_val
                        break
                if len(input_val) >= 3:
                   print("More than two paths specified. Please include up to 2 separate paths (one for images, one for structures)\n>> ") 

            except KeyboardInterrupt:
                quit()
            except: pass       

    if args.save_dir_struct is not None and not os.path.exists(args.save_dir_struct):
        while True:
            try:
                input_val = input(f"\n Could not find {args.save_dir_struct} for structure nrrd saving. Please re-specify the directory to save nrrd images to for structures. If structure(s) are to be saved in a different folder, please separate directories with a comma.\n>> ")
                input_val = list(input_val.split(","))
                if len(input_val) == 1:
                    if os.path.exists(input_val[0].strip()):
                        save_dir = input_val[0].strip()
                        print(f"saving all nrrd data into {input_val}")
                        break
                if len(input_val) == 2:
                    input_val[0] = input_val[0].strip()    #remove any white space in path
                    input_val[1] = input_val[1].strip()    #remove any white space in path
                    if os.path.exists(input_val[0]) and os.path.exists(input_val[1]):
                        print(f"saving image nrrd data into {input_val[0]}")
                        print(f"saving structure nrrd data into {input_val[1]}")
                        save_dir = input_val
                        break
                if len(input_val) >= 3:
                   print("More than two paths specified. Please include up to 2 separate paths (one for images, one for structures)") 

            except KeyboardInterrupt:
                quit()
            except: pass      
    elif args.save_dir_struct is not None and os.path.exists(args.save_dir_struct) and len(save_dir) == 1:
        save_dir = [save_dir[0], args.save_dir_struct]          
    if args.modality == None:
        print("no modality specified, using CT as default.")
    
    convert_all_dicoms_to_nrrd(img_dir, structure_dir, modality, save_dir)    
