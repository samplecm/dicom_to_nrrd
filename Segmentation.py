import numpy as np 
import pydicom
from pydicom import dcmread
import pickle
import os
import glob
import pathlib
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from PIL import Image, ImageDraw
from operator import itemgetter
import matplotlib.pyplot as plt
from fastDamerauLevenshtein import damerauLevenshtein
from Contours import Contours
import Chopper


def SubsegmentPointFiller(contours):
    for i, subsegment in enumerate(contours): 
        for j, slice in enumerate(subsegment): 
            contours[i][j] = PointFiller(slice)
    return contours        

def PointFiller(slice):
    #Add interpolated points between points far apart in a contour 
    if len(slice) == 0:
        return slice
    noLargeSpaces = True
    i = 0
    while i < len(slice)-1:
        if abs(slice[i][0] - slice[i+1][0]) > 2:
            x = 0.5 * (slice[i][0] + slice[i+1][0])
            y = 0.5 * (slice[i][1] + slice[i+1][1])
            z = slice[i][2]
            newPoint = [x,y,z]
            slice.insert(i+1, newPoint)
            i = max(i-2,0)
        else:
            i += 1
    i=0        
    while i < len(slice)-1:
        if abs(slice[i][1] - slice[i+1][1]) > 2:
            x = 0.5 * (slice[i][0] + slice[i+1][0])
            y = 0.5 * (slice[i][1] + slice[i+1][1])
            z = slice[i][2]
            newPoint = [x,y,z]
            slice.insert(i+1, newPoint)
            i = max(i-2, 0)
        else:
            i += 1  
          
    return slice            





def ProcessData(organ : str, numCuts):     
    print("In ProcessData")

    filesFolder = os.path.join(pathlib.Path(__file__).parent.absolute(), "Patients")
    filesFolderList = sorted(os.listdir(os.path.join(pathlib.Path(__file__).parent.absolute(), "Patients")))
    patientFolders = []

    #Go through and get all patient file names in a list, and create a new folder for storing data in for each: 
    for folder in filesFolderList:
        patientFolders.append(folder)
        if not os.path.isdir(os.path.join(pathlib.Path(__file__).parent.absolute(), "Processed_Data/") + folder):
            os.mkdir(os.path.join(pathlib.Path(__file__).parent.absolute(), "Processed_Data/") + folder)


    for p in range(len(patientFolders)):

        patientFiles = sorted(glob.glob(os.path.join(filesFolder, patientFolders[p], "*")))
        #get the RTSTRUCT dicom file and get patient 's CT scans: 
        for fileName in patientFiles:
            if "STRUCT" in fileName:
                structFile = fileName   
        # for elem in dcmread(structFile).iterall():
        #     print(elem)
        structsMeta = dcmread(structFile).data_element("ROIContourSequence")
        structure, structureROINum= FindStructure(dcmread(structFile).data_element("StructureSetROISequence"), organ)  
        

        if structureROINum != 1111:
            contourList = []
            for contourInfo in structsMeta:
                if contourInfo.get("ReferencedROINumber") == structureROINum: #get the correct contour for the given organ
                    for contoursequence in contourInfo.ContourSequence: 
                        contourList.append(contoursequence.ContourData)
                        #But this is easier to work with if we convert from a 1d to a 2d list for contours ( [ [x1,y1,z1], [x2,y2,z2] ... ] )
                        numContourPoints = 0 
                        tempContour = []
                        i = 0
                        while i < len(contourList[-1]):
                            x = float(contourList[-1][i])
                            y = float(contourList[-1][i + 1])
                            z = float(contourList[-1][i + 2])
                            tempContour.append([x, y, z ])
                            i += 3
                            if (numContourPoints > 0):
                                pointListy = np.vstack((pointListy, np.array([x,y,z])))
                            else:
                                pointListy = np.array([x,y,z])   
                            numContourPoints+=1         
                        contourList[-1] = tempContour 
                            

            #Now save the contours as objects along with the organ name, dicom structure name
            contours = Contours(organ, structure, contourList)

            Chopper.OrganChopper(contours, numCuts, organ)
            Visuals.plotStructure(contours.segmentedContours18)
            saveContoursPath = os.path.join(pathlib.Path(__file__).parent.absolute(), "Processed_Data/") + patientFolders[p] + "/contours_" + structure + ".txt"
            with open(saveContoursPath, "wb") as fp:
                pickle.dump(contours, fp)  




def FindStructure(metadata, organ, invalidStructures = []):
    #Here we take the string for the desired structure (organ) and find the matching structure for each patient. 
    #The algorithm is to first make sure that the organ has a substring matching the ROI with at least 3 characters,
    #then out of the remaining possiblities, find top 3 closest fits with damerau levenshtein algorithm, then check to make sure that they are allowed to match according to rules defined in AllowedToMatch(). There should then ideally
    # be only one possible match, but if there are two, take the first in the list.   
    #Get a list of all structures in structure set: 
    unfilteredStructures = []
    for element in metadata:
        if element.get("ROIName").lower() not in invalidStructures:
            unfilteredStructures.append(element.get("ROIName").lower())           
    #Now find which is the best fit.
    #First filter out any structures without at least a 3 character substring in common
    structures = []
    for structure in unfilteredStructures:
        valid = True
        if LongestSubstring(structure, organ) < 3:
            valid = False
        if not AllowedToMatch(organ, structure): 
            valid = False
        #Add to structures if valid
        if valid:
            structures.append(structure)
    #Now test string closeness to find
    closestStrings = [["No Match",100],["No Match",100],["No Match",100]] #has to be in the top 3 closest strings to check next conditions
    for structure in structures:
        closeness = StringDistance(structure, organ)
        closestStrings.sort(key=itemgetter(1)) #Sort by the closeness value, and not the structure names
        for value in range(len(closestStrings)):
            if closeness < closestStrings[value][1]: #If closer than a value already in the top 3
                closestStrings[value] = [structure, closeness]
                break
    
    if closestStrings[0][1] == 100:
        return "No Match", 1111    
    #Now return the organ that is remaining and has closest string
    for element in metadata:
        if element.get("ROIName").lower() == closestStrings[0][0]:
            roiNumber = element.get("ROINumber")
    try:
        return closestStrings[0][0], roiNumber 
    except:
        return "No Match", 1111 #error code for unfound match.    

def AllowedToMatch(s1, s2):
    s1 = s1.lower()
    s2 = s2.lower()
    allowed = True
    keywords = []
    #You can't have only one organ with one of these keywords...
    keywords.append("prv")
    keywords.append("brain")
    keywords.append("ptv")
    keywords.append("stem")
    keywords.append("node")
    keywords.append("cord")
    keywords.append("chi")
    keywords.append("opt")
    keywords.append("oral")
    keywords.append("nerv")
    keywords.append("par")
    keywords.append("globe")
    keywords.append("lip")
    keywords.append("cav")
    keywords.append("sub")
    keywords.append("test")
    keywords.append("fact")
    #keywords can't be in only one of two string names: 
    for keyword in keywords:
        num = 0
        if keyword in s1:
            num += 1
        if keyword in s2:
            num += 1
        if num == 1:
            allowed = False        

    #Cant have left and no l in other, or rightt and no r
    if "left" in s1:
        if "l" not in s2:
            allowed = False      
    if "left" in s2:
        if "l" not in s1:
            allowed = False    
    #its tricky matching up left and right organs sometimes with all the conventions used... this makes sure that both are left or both are right
    if ("_l_" in s1) or (" l " in s1) or  (" l-" in s1) or ("-l-" in s1) or (" l_" in s1) or ("_l " in s1) or ("-l " in s1) or ("left" in s1) or ("l " == s1[0:2]) or ("_lt_" in s1) or (" lt " in s1) or  (" lt-" in s1) or ("-lt-" in s1) or (" lt_" in s1) or ("_lt " in s1) or ("-lt " in s1) or ("lt " == s1[0:3]):
        if not (("lpar" in s2) or ("lsub" in s2) or ("_l_" in s2) or (" l " in s2) or  (" l-" in s2) or ("-l-" in s2) or (" l_" in s2) or ("_l " in s2) or ("-l " in s2) or ("left" in s2) or ("l " == s2[0:2])or ("_lt_" in s2) or (" lt " in s2) or  (" lt-" in s2) or ("-lt-" in s2) or (" lt_" in s2) or ("_lt " in s2) or ("-lt " in s2) or ("lt " == s2[0:3])):   
            allowed = False  
    if (("_l_" in s2) or (" l " in s2) or  (" l-" in s2) or ("-l-" in s2) or (" l_" in s2) or ("_l " in s2) or ("-l " in s2) or ("left" in s2) or ("l " == s2[0:2])or ("_lt_" in s2) or (" lt " in s2) or  (" lt-" in s2) or ("-lt-" in s2) or (" lt_" in s2) or ("_lt " in s2) or ("-lt " in s2)or ("lt " == s2[0:3])):  
        if not (("lpar" in s1) or ("lsub" in s1) or ("_l_" in s1) or (" l " in s1) or  (" l-" in s1) or ("-l-" in s1) or (" l_" in s1) or ("_l " in s1) or ("-l " in s1) or ("left" in s1) or ("l " == s1[0:2]) or ("_lt_" in s1) or (" lt " in s1) or  (" lt-" in s1) or ("-lt-" in s1) or (" lt_" in s1) or ("_lt " in s1) or ("-lt " in s1) or ("lt " == s1[0:3])):
            allowed = False        
    
    if ("_r_" in s1) or (" r " in s1) or  (" r-" in s1) or ("-r-" in s1) or (" r_" in s1) or ("_r " in s1) or ("-r " in s1) or ("right" in s1) or ("r " == s1[0:2])or ("_rt_" in s1) or (" rt " in s1) or  (" rt-" in s1) or ("-rt-" in s1) or (" rt_" in s1) or ("_rt " in s1) or ("-rt " in s1)or ("right" in s1):
        if not (("rpar" in s2) or ("rsub" in s2) or ("_r_" in s2) or (" r " in s2) or  (" r-" in s2) or ("-r-" in s2) or (" r_" in s2) or ("_r " in s2) or ("-r " in s2) or ("right" in s2) or ("r " == s2[0:2]) or ("_rt_" in s2) or (" rt " in s2) or  (" rt-" in s2) or ("-rt-" in s2) or (" rt_" in s2) or ("_rt " in s2) or ("-rt" in s2) ):   
            allowed = False
    if (("_r_" in s2) or (" r " in s2) or  (" r-" in s2) or ("-r-" in s2) or (" r_" in s2) or ("_r " in s2) or ("-r " in s2) or ("right" in s2) or ("r " == s2[0:2]) or ("_rt_" in s2) or (" rt " in s2) or  (" rt-" in s2) or ("-rt-" in s2) or (" rt_" in s2) or ("_rt " in s2) or ("-rt" in s2) ): 
        if not (("rpar" in s1) or ("rsub" in s1) or ("_r_" in s1) or (" r " in s1) or  (" r-" in s1) or ("-r-" in s1) or (" r_" in s1) or ("_r " in s1) or ("-r " in s1) or ("right" in s1) or ("r " == s1[0:2])or ("_rt_" in s1) or (" rt " in s1) or  (" rt-" in s1) or ("-rt-" in s1) or (" rt_" in s1) or ("_rt " in s1) or ("-rt " in s1)):
            allowed = False
    return allowed


def StringDistance(s1, s2):
    return damerauLevenshtein(s1,s2,similarity=False)

    return d[lenstr1-1,lenstr2-1]
def LongestSubstring(s1,s2):
    m = len(s1)
    n = len(s2)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if s1[i] == s2[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(s1[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(s1[i-c+1:i+1])
    return longest  



if __name__ == "__main__":
    numCuts = [2, 1, 2]
    ProcessData("Left Parotid", numCuts)

    #Now loop through contours and plot
