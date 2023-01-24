from PIL import Image, ImageDraw
from shapely.geometry.polygon import Polygon
import numpy as np
import copy
import matplotlib.pyplot as plt
import mask_interpolation
def get_all_structure_masks(img_series, structures_dict):
    new_dict = {}
    coords_array = img_series.coords_array_img
    for structure_name in structures_dict:
        contours = structures_dict[structure_name]
        mask_array = get_contour_masks(contours, coords_array)
        new_dict[structure_name] = mask_array
    return new_dict
def get_contour_masks(contours, array):
    num_slices, len1, len2 = array.shape[1:]
    contour_masks = np.zeros((num_slices,len1,len2)) 
    contours = cartesian_to_pixel_coordinates(clone_list(contours), array)

    for idx in range(num_slices):#loop through all slices creating a mask for the contours
        contour_mask_filled = Image.new('L', (len2, len1), 0)
        slice_found = False
        
        for i, contour in enumerate(contours):
            if len(contour) < 3:
                continue
            if abs(int(round(contour[0][2], 2)*100) - int(round(array[2,idx,0,0], 2)*100)) < 2: #if contour is on the current slice and only want to make new image when starting with first island
                slice_found = True
                contourPoints = []
                for point in contour:
                    contourPoints.append((int(point[0]), int(point[1]))) #changed
                ImageDraw.Draw(contour_mask_filled).polygon(contourPoints, outline= 1, fill = 1)         

                # plt.imshow(contour_mask_filled)
                # plt.show()
                # unshown=False
                # print("")
            if slice_found == True:
                break        
 
        contour_masks[idx, :,:] = np.array(contour_mask_filled)       
    #lastly, want to interpolate any missing slices. 
    for idx in range(1,contour_masks.shape[0]-1):
        if np.amax(contour_masks[idx,:,:])==0 and np.amax(contour_masks[idx-1,:,:])==1 and np.amax(contour_masks[idx+1,:,:]) == 1:
            interp_img = mask_interpolation.interp_shape(contour_masks[idx-1,:,:], contour_masks[idx+1,:,:], 0.5)
            contour_masks[idx,:,:] = interp_img
            # fig,ax = plt.subplots(3)
            # ax[0].imshow(contour_masks[idx-1,:,:])
            # ax[1].imshow(contour_masks[idx,:,:])
            # ax[2].imshow(contour_masks[idx+1,:,:])
            # plt.show()
            # print("")
    return contour_masks                

def clone_list(list):
    listCopy = copy.deepcopy(list)
    return listCopy  

def cartesian_to_pixel_coordinates(contours, array):
    #convert x and y values for a contour into the pixel indices where they are on the pet array
    xVals = array[0,0,0,:]
    yVals = array[1,0,:,0]
    for contour in contours: 
        if len(contour) == 0: continue
        for point in contour:
            point[0] = min(range(len(xVals)), key=lambda i: abs(xVals[i]-point[0]))
            point[1] = min(range(len(yVals)), key=lambda i: abs(yVals[i]-point[1]))
    return contours  
