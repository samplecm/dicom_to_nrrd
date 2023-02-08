import math
from typing import final
from Contours import Contours
import copy
import Segmentation
import time


def clone_list(list):
    listCopy = copy.deepcopy(list)
    return listCopy  

def mandible_chopper(contours):
    print("Sub-segmenting Mandible... ")
    finalContours = []
    contoursX = []

    contoursX = x_chop(contours.whole_roi_img, 5)
    # contoursX = list(filter(lambda layer: layer != [[]], contours))
    # contours.sort(key=Chopper.get_z_val)  
    for j in range(len(contoursX)):
        temp = z_chop(contoursX[j], 1)
        for j in range(len(temp)):
            finalContours.append(copy.deepcopy(temp[j]))

     
    
    #Now I want to separate any connected islands into separate contours, so they can be properly rendered 
    #so the structure of the list will now include one additional nested list. each subsegment slice is now a list of all slice islands first (most will be length 1)
    finalContours = separate_islands(clone_list(finalContours))
    finalContours = Segmentation.SubsegmentPointFiller(finalContours)
    contours.segmentedContours = finalContours

def cord_chopper(contours: Contours):
    #this  function chops the spinal cord axially into chunks of 2cm vertical length. This starts from the top of the cord, and any remaining <2cm chunk left at the bottom is ignored. 
    print("Subsegmenting spinal cord ... ")
    max_z = -10000
    min_z = 10000
    orig_contours = copy.deepcopy(contours.whole_roi_img)
    for slice in orig_contours:
        z = slice[0][0][2]
        if z > max_z:
            max_z = z
        if z < min_z: 
            min_z = z    
    z = max_z - 20 #in mm 
    new_contours = []
    while z > min_z:
        #first check if this z value already has a contour
        slice_exists = False
        for i in range(len(orig_contours)):
            if orig_contours[i][0][0][2] == z:
                slice_exists = True
                break

        if not slice_exists:    

            contoursZ = closest_contour_z(z, orig_contours) #Get the closest contour above and below the cut
            newContour = []
            for j in range(len(orig_contours[contoursZ[0]][0])): #Go over all points in closest contour
                point1 = orig_contours[contoursZ[0]][0][j]      
                point2 = closest_point(point1, orig_contours[contoursZ[1]])    #Now get the closest point in second closest contour
                newPoint = InterpolateXY(point1, point2, z)    #now interpolate between the two

                #add new point to new contour
                newContour.append(newPoint)
            #add this new contour to newContoursList
            orig_contours.append([newContour])    
            
        z -= 20 #in mm     

    orig_contours.sort(key=get_z_val)  
    z = max_z - 20
    while z > min_z:
        new_contours.append([])
        for j in range(len(orig_contours)): #Go over all points in closest contour    
            if orig_contours[j][0][0][2] >= z and orig_contours[j][0][0][2] <= z + 20:
                new_contours[-1].append(copy.deepcopy(orig_contours[j]))
        z -= 20        
    contours.segmentedContours = new_contours
    print("")


def organ_chopper(contours : Contours, numCuts):
    print(f"Subsegmenting {contours.roiName} ...")
    #chop into 18ths
    numCutsX = numCuts[0]
    numCutsY = numCuts[1]
    numCutsZ = numCuts[2]


    z_chop(contours, numCutsZ) #Make axial cuts and store in contours.segmentedContours3

    contoursY = []
    if numCutsY > 0:
        for i in range(len(contours.segmentedContours3)):
            temp = y_chop(contours.segmentedContours3[i], numCutsY)
            for j in range(len(temp)):
                contoursY.append(copy.deepcopy(temp[j]))
    else:
        contoursY = copy.deepcopy(contours.segmentedContours3)
    contours.segmentedContours6 = contoursY
    finalContours = []
    if numCutsX > 0:
        for i in range(len(contoursY)):
            temp = x_chop(contoursY[i], numCutsX)
            for j in range(len(temp)):
                finalContours.append(copy.deepcopy(temp[j]))
    else:
        finalContours = contoursY

    
    #Now I want to separate any connected islands into separate contours, so they can be properly rendered 
    #so the structure of the list will now include one additional nested list. each subsegment slice is now a list of all slice islands first (most will be length 1)
    finalContours = separate_islands(clone_list(finalContours))
    finalContours = Segmentation.SubsegmentPointFiller(finalContours)
    finalContours = ReOrderParotids(finalContours, contours.roiName, numCuts)
    contours.segmentedContours = finalContours
 
def separate_islands(contours):
    #First loop through slices and check for cases of lines in x or y direction that traverse both ways in the contour.
    newContours = []
    for s, subsegment in enumerate(contours):
        for i, slice in enumerate(subsegment):  
            if len(slice) == 0:
                continue
            bisect_indices = [] #this stores potential breakup indices. 
            xLines = []
            yLines = []
            point_idx = 0
            while point_idx < len(slice)-1:
                #look for lines in x/y direction that traverse both ways. 
                lineEnd_idx = point_idx #default value for iterating

                if slice[point_idx][0] == slice[point_idx+1][0]:
                    lineEnd_idx = point_idx+1
                    #first see if line stretches for more than the next point:
                    for j in range(point_idx+2, len(slice)):
                        if slice[point_idx][0] == slice[j][0]:
                            lineEnd_idx = j
                        else:
                            break    
                        
                    found = False
                    for item in xLines: #check if already a line at this location, and if there is append this index to that list
                        if abs(item[0] - slice[point_idx][0]) < 0.1: #if within 0.1mm of each other
                            item.append([[point_idx, lineEnd_idx], slice[lineEnd_idx][1], slice[point_idx][1]])
                            found=True
                    if found==False:
                        xLines.append([slice[point_idx][0], [[point_idx, lineEnd_idx], slice[lineEnd_idx][1],slice[point_idx][1]]])  #values are xLocation, displacement of y values


                if slice[point_idx][1] == slice[point_idx+1][1]:
                    lineEnd_idx = point_idx+1
                    #first see if line stretches for more than the next point:
                    for j in range(point_idx+2, len(slice)):
                        if slice[point_idx][1] == slice[j][1]:
                            lineEnd_idx = j
                        else:
                            break
                    found = False
                    for item in yLines: #check if already a line at this location, and if there is append this index to that list
                        if abs(item[0] - slice[point_idx][1]) < 0.1:
                            item.append([[point_idx, lineEnd_idx], slice[lineEnd_idx][0],slice[point_idx][0]])
                            found=True
                    if found==False: 
                        yLines.append([slice[point_idx][1], [[point_idx, lineEnd_idx], slice[lineEnd_idx][0],slice[point_idx][0]]])
                point_idx += 1 + (lineEnd_idx - point_idx)        
            #now check for bisecting lines at same x or y:
            # for point_idx in range(len(slice)-1):
            #     #look for lines in x/y direction that traverse both ways. 
            #     if slice[point_idx][0] == slice[point_idx+1][0]:
            #         found = False
            #         for item in xLines: #check if already a line at this location, and if there is append this index to that list
            #             if abs(item[0] - slice[point_idx][0]) < 0.1:
            #                 item.append([point_idx, slice[point_idx+1][1], slice[point_idx][1]])
            #                 found=True
            #         if found==False:
            #             xLines.append([slice[point_idx][0], [point_idx, slice[point_idx+1][1],slice[point_idx][1]]])  #values are xLocation, displacement of y values


            #     if slice[point_idx][1] == slice[point_idx+1][1]:
            #         found = False
            #         for item in yLines: #check if already a line at this location, and if there is append this index to that list
            #             if abs(item[0] - slice[point_idx][1]) < 0.1:
            #                 item.append([point_idx, slice[point_idx+1][0],slice[point_idx][0]])
            #                 found=True
            #         if found==False: 
            #             yLines.append([slice[point_idx][1], [point_idx, slice[point_idx+1][0],slice[point_idx][0]]])
            # #now check for bisecting lines at same x or y:
            lines = []
            for xVal in xLines:
                if len(xVal) <= 2:
                    continue
                for idx in range(1, len(xVal)):
                    lineSegment = xVal[idx]
                    y11 = lineSegment[1]
                    y12 = lineSegment[2]
                    min_y1 = min([y11, y12])
                    max_y1 = max([y11, y12])
                    #Now check if this oversects with any other lines along this x point

                    for line in lines: 
                        y21 = line[2]
                        y22 = line[1]
                        min_y2 = min([y21, y22])
                        max_y2 = max([y21, y22])

                        minY = min([min_y1, min_y2])
                        if minY == min_y1 and max_y1 > min_y2:
                            #they do intersect! Add to the potential breakup list
                            bisect_indices.append([line[0], lineSegment[0]])
                            break
                        elif minY == min_y2 and max_y2 > min_y1:
                            #they do intersect! Add to the potential breakup list
                            bisect_indices.append([line[0], lineSegment[0]])
                            break
                    lines.append(lineSegment)    
                 
            for yVal in yLines:
                if len(yVal) < 2:
                    continue
                lines = []   
                for idx in range(1, len(yVal)):
                    lineSegment = yVal[idx]
                    x11 = lineSegment[1]
                    x12 = lineSegment[2]
                    min_x1 = min([x11, x12])
                    max_x1 = max([x11, x12])
                    #Now check if this oversects with any other lines along this x point

                    for line in lines: 
                        x21 = line[2]
                        x22 = line[1]
                        min_x2 = min([x21, x22])
                        max_x2 = max([x21, x22])

                        minX = min([min_x1, min_x2])
                        if minX == min_x1 and max_x1 > min_x2:
                            #they do intersect! Add to the potential breakup list
                            bisect_indices.append([line[0], lineSegment[0]])
                            break
                        elif minX == min_x2 and max_x2 > min_x1:
                            #they do intersect! Add to the potential breakup list
                            bisect_indices.append([line[0], lineSegment[0]])
                            break
                    lines.append(lineSegment)       

            #Now go through bisect_indices and make first appropriate division
            restructured = False
            for indices in bisect_indices:
                indices.sort(key=lambda x: x[0])
                # #First make sure they are more than 10 points apart:
                # spacedApart = indices[1] - indices[0]
                # if spacedApart < 10:
                #     continue
                #Now take the points between these values and divide them into their own new list 
                newList1 = []
                newList2 = []
                for point_idx in range(len(slice)):
                    if point_idx <= indices[0][0]:
                        newList1.append(slice[point_idx])
                    elif point_idx <= indices[1][0] and point_idx >= indices[0][1]:
                        newList2.append(slice[point_idx])
                    elif point_idx >= indices[1][1]:
                        newList1.append(slice[point_idx])     
                restructured = True
                break          
            #for testing:

            #now the new format:
            if restructured:
                if len(newList1) > len(newList2):
                    # slice.append(newList1)
                    # slice.append(newList2)
                    newList1 = Segmentation.PointFiller(newList1)
                    contours[s][i] = closed_looper(newList1)
                else:  
                    # slice.append(newList1)
                    # slice.append(newList2)  
                    newList2 = Segmentation.PointFiller(newList2)
                    contours[s][i] = closed_looper(newList2)
    return contours
                    
                    




def z_chop(contours, numCutsZ):
    
    newContoursList = []
    if type(contours) != list:
        contours_array = contours.whole_roi_img
        zCuts = best_cut_z(contours.whole_roi_img, numCutsZ) #Get the list of z values where structure is to be chopped
    else:
        contours_array = contours    
        zCuts = best_cut_z(contours, numCutsZ) #Get the list of z values where structure is to be chopped

    #commenting out part for adding in interpolated contours (biases dose calc)
    # for i in range(len(zCuts)):
    #     contoursZ = closest_contour_z(zCuts[i], contours) #Get the closest contour above and below the cut
    #     newContour = []
    #     for j in range(len(contours_array[contoursZ[0]][0])): #Go over all points in closest contour

    #         point1 = contours_array[contoursZ[0]][0][j]
    #         #Now get the closest point in second closest contour
    #         point2 = closest_point(point1, contours_array[contoursZ[1]])
    #         #now interpolate between the two
    #         newPoint = InterpolateXY(point1, point2, zCuts[i])

    #         #add new point to new contour
    #         newContour.append(newPoint)
    #     #add this new contour to newContoursList
    #     newContoursList.append([newContour])    
    #Now add all the original contours to the newcontourslist and sort it by z value
    for i in range(len(contours_array)):
        newContoursList.append(contours_array[i])
    newContoursList.sort(key=get_z_val)

    #Now create a list of contour lists for each axial slice
    axialContours = []
    for i in range(numCutsZ + 1):
        axialContours.append([])

    for i in range(len(newContoursList)):
        for j in range(numCutsZ):
            if newContoursList[i][0][0][2] <= zCuts[j]: #z value of first point in contour
                if j == 0:
                    axialContours[j].append(copy.deepcopy(newContoursList[i]))
                elif newContoursList[i][0][0][2] >= zCuts[j-1]: #Needs to be larger than the cut below to be in a middle or top region
                    axialContours[j].append(copy.deepcopy(newContoursList[i]))

        if newContoursList[i][0][0][2] >= zCuts[-1]: 
            axialContours[-1].append(copy.deepcopy(newContoursList[i]))        

            

    if type(contours) != list:        
        contours.segmentedContours3 =  axialContours  

    return axialContours           
            

def get_z_val(e): 
    #function for returning z value of a contour
    if len(e[0]) == 0:
        raise ValueError("Contour of 0 points has no z value")
    return e[0][0][2]

def GetContourArea(contour):
    area = 0
    if len(contour) == 0:
        return 0

    for island in contour[0:1]: #taking first island only here
        for i in range(len(island)-1):
            area += island[i][0] * island[i+1][1] - island[i+1][0] * island[i][1]
        area += island[-1][0] * island[0][1] - island[0][0] * island[-1][1]    
    return 0.5 * abs(area)

def best_cut_z(contours, numCuts):
    zCuts = [] 
    numContours = len(contours)
    totalVolume = 0
    deltaZ = abs(contours[0][0][0][2]- contours[1][0][0][2]) #axial distance between slices 
    contourAreas = []
    for i in range(numContours-1):
        area=0
        area += GetContourArea(contours[i])
        contourAreas.append(area)
        if i != numContours - 1:  #Dont include the final area in the volume calculation, since it would be over-extending
            totalVolume += area * deltaZ

    #Now find the correct cut spots which result in equal volume subsegments  
    for i in range(1, numCuts + 1):
        volumeGoal = i / (numCuts + 1)    #This is the ratio of volume that each subsegment should consist of (python integer division returns a float)
        contIndex = 0
        subVolume = 0
        while (subVolume < volumeGoal * totalVolume): #This loop will result in contIndex which tells us that the new slice should be somewhere between contourIndex and contourIndex - 1
            subVolume += contourAreas[contIndex] * deltaZ
            contIndex += 1
        volumeBelow = 0
        for j in range(contIndex - 1):
            avgArea = 0.5 * (contourAreas[j] + contourAreas[j+1])#first get the average area between two contours, used to approximate volume between the two (makes no difference while using equal slice diffs)
            volumeBelow += avgArea * deltaZ
        #now get the avg area for the slicing region: 
        avgArea = 0.5 * (contourAreas[contIndex - 1] + contourAreas[contIndex])

        zSlice = contours[contIndex - 1][0][0][2] + (volumeGoal * totalVolume - volumeBelow) / avgArea
        zCuts.append(zSlice)
    return zCuts    


def closest_contour_z(z, contours):
    if type(contours) != list:
        contours = contours.whole_roi_img
        
    temp = 1000
    closestContours = [0,0]

    #First find closest contour
    for i in range(len(contours)):
        if len(contours[i][0]) == 0:
            continue
        contourDistance = abs(contours[i][0][0][2] - z)    
        if contourDistance < temp:
            closestContours[0] = i
            temp = contourDistance


    #Now get second closest contour
    temp = 1000
    for i in range(len(contours)):
        if len(contours[i][0]) == 0:
            continue
        if i == closestContours[0]: #Can't count the closest contour again
            continue        
        contourDistance = abs(contours[i][0][0][2] - z)

        if contourDistance < temp:
            closestContours[1] = i
            temp = contourDistance

    return closestContours        


def closest_point(point, contour, island_idx=None):
    if island_idx == None:
        m = 1000
        closest_pointIdx1 = 1000
        closest_pointIdx2 = 1000
        for p, polygon in enumerate(contour):
            for i in range(len(polygon)):
                diff = math.sqrt((point[0] - polygon[i][0])**2 + (point[1] - polygon[i][1])**2)
                if diff < m:
                    closest_pointIdx1 = p
                    closest_pointIdx2 = i
                    m = diff
        if closest_pointIdx1 == 1000:
            raise Exception("Could not find closest point.")
        return contour[closest_pointIdx1][closest_pointIdx2]    
    else:   
        m = 1000
        closest_pointIdx = 1000
        try:
            for i in range(len(contour)):
                diff = math.sqrt((point[0] - contour[i][0])**2 + (point[1] - contour[i][1])**2)
                if diff < m:
                    closest_pointIdx = i
                    m = diff
            if closest_pointIdx == 1000:
                raise Exception("Could not find closest point.")
            return contour[closest_pointIdx]        
        except IndexError:
            return None        
        

def InterpolateXY(point1, point2, z):
    if point1[2] == z:
        return point1
    elif point2[2] == z: 
        return point2

    xSlope = (point2[0] - point1[0]) / (point2[2] - point1[2])
    ySlope = (point2[1] - point1[1]) / (point2[2] - point1[2])   

    newX = point1[0] + xSlope * (z - point1[2])
    newY = point1[1] + ySlope * (z - point1[2])

    newPoint = [newX, newY, z]

    return newPoint


def y_chop(contours, numCutsY):  
    yCuts = best_cut_y(contours, numCutsY)

    #add intersection points to contours
    for i in range(len(contours)):
        for j in range(len(yCuts)):
            contours[i] = add_intersections_y(contours[i], yCuts[j])
            #contours[i] = closed_looper(contours[i])

    #Now divide into separate parts
    finalContours = []
    divisions = [] #list for each y division for current contour

    #make the list the correct size so there is an item for each y div
    for div in range(len(yCuts)+ 1):
        finalContours.append([])

    for i in range(len(contours)):
        divisions.clear()
        for y in range(len(yCuts)+1):
            divisions.append([])  
        for y in range(len(yCuts)+1): #a section for each cut, +1
            for j in range(len(contours[i][0])):
                if y == 0:
                    if contours[i][0][j][1] <= yCuts[y]:
                        divisions[y].append(contours[i][0][j])

                elif y == len(yCuts):
                    if contours[i][0][j][1] >= yCuts[y-1]:
                        divisions[y].append(contours[i][0][j])

                else:
                    if contours[i][0][j][1] >= yCuts[y-1] and contours[i][0][j][1] <= yCuts[y]:
                        divisions[y].append(contours[i][0][j])
   
        for div in range(len(divisions)):
            #temp = closed_looper(divisions[div])
            if len(divisions[div]) == 0:
                continue
            finalContours[div].append(copy.deepcopy([divisions[div]])) #extra brackets to separate different contours
    return finalContours        

def best_cut_y(contours, numCutsY):
    errorTolerance = 0.001
    area = 0
    maxY = -1000
    minY = 1000
    yCuts= []

    #get total area of contours
    for contour in contours:
        area += GetContourArea(contour)

    for cut in range(numCutsY):
        areaGoal = (cut + 1) / (numCutsY + 1)

        for j in range(len(contours)):
            if len(contours[j]) > 0:
                for row in range(len(contours[j][0])):
                    if contours[j][0][row][1] > maxY:
                        maxY = contours[j][0][row][1]
                    if contours[j][0][row][1] < minY:
                        minY = contours[j][0][row][1]

        tempContours = []
        cutContours = [] 

        error = 1000
        time1 = time.time()
        while error > errorTolerance:     
            tempContours.clear()
            yCut = (minY + maxY) / 2
            newArea = 0

            for i in range(len(contours)):
                cutContours.clear()
                temp = add_intersections_y(contours[i], yCut)
                #temp = closed_looper(temp)
                tempContours.append(temp)
                for j in range(len(tempContours[-1][0])):
                    if tempContours[-1][0][j][1] <= yCut:
                        cutContours.append(tempContours[-1][0][j].copy())
                if len(cutContours) > 0:
                    #cutContours = closed_looper(cutContours)
                    newArea += GetContourArea([cutContours])
            error = abs((newArea / area) - areaGoal)   
            if (newArea / area < areaGoal):
                minY = yCut
            elif (newArea / area > areaGoal):
                maxY = yCut
            if time.time() - time1 > 5:
                break #stop after 5 seconds    

             
        yCuts.append(yCut) 
    return yCuts
def ReOrderSMs(contours, name : str):
    reOrderedList = []
    if "left" in name.lower():
        reOrderedList.append(clone_list(contours[0]))
        reOrderedList.append(clone_list(contours[4]))
        reOrderedList.append(clone_list(contours[2]))
        reOrderedList.append(clone_list(contours[6]))
        reOrderedList.append(clone_list(contours[1]))
        reOrderedList.append(clone_list(contours[5]))
        reOrderedList.append(clone_list(contours[3]))
        reOrderedList.append(clone_list(contours[7]))
    if "right" in name.lower():
        reOrderedList.append(clone_list(contours[1]))
        reOrderedList.append(clone_list(contours[5]))
        reOrderedList.append(clone_list(contours[3]))
        reOrderedList.append(clone_list(contours[7]))
        reOrderedList.append(clone_list(contours[0]))
        reOrderedList.append(clone_list(contours[4]))
        reOrderedList.append(clone_list(contours[2]))
        reOrderedList.append(clone_list(contours[6]))  
    return reOrderedList      
        



def ReOrderParotids(contours, organName, numCuts):
    #Reorder from inferior --> superior, medial --> lateral, anterior --> posterior
    j = 0
    final_contours = []

    if 'l' in organName.lower():
        final_contours.append(contours[0])
        final_contours.append(contours[6])
        final_contours.append(contours[12])
        final_contours.append(contours[1])
        final_contours.append(contours[7])
        final_contours.append(contours[13])
        final_contours.append(contours[2])
        final_contours.append(contours[8])
        final_contours.append(contours[14])
        final_contours.append(contours[3])
        final_contours.append(contours[9])
        final_contours.append(contours[15])
        final_contours.append(contours[4])
        final_contours.append(contours[10])
        final_contours.append(contours[16])
        final_contours.append(contours[5])
        final_contours.append(contours[11])
        final_contours.append(contours[17])

    else: 
        final_contours.append(contours[2])
        final_contours.append(contours[8])
        final_contours.append(contours[14])
        final_contours.append(contours[1])
        final_contours.append(contours[7])
        final_contours.append(contours[13])
        final_contours.append(contours[0])
        final_contours.append(contours[6])
        final_contours.append(contours[12])
        final_contours.append(contours[5])
        final_contours.append(contours[11])
        final_contours.append(contours[17])
        final_contours.append(contours[4])
        final_contours.append(contours[10])
        final_contours.append(contours[16])
        final_contours.append(contours[3])
        final_contours.append(contours[9])
        final_contours.append(contours[15])
        
        
    return final_contours

def add_intersections_y(contour, yCut):
    finalContour = copy.deepcopy(contour)
    numAdded =1 #start at one, increment after adding each point, to keep track of where to add additional point (add to index)
    z = contour[0][0][2]
    #index 0 outside of loop:
    if contour[0][0][1]  > yCut and contour[0][-1][1] < yCut:
        if contour[0][0][0] == contour[-1][0][0]: #if xs are same don't need to interpolate
            xNew = contour[0][0][0]
        else:        
            m = (contour[0][0][1]- contour[0][-1][1]) / (contour[0][0][0] - contour[0][-1][0])
            xNew = (yCut - contour[0][-1][1]) / m + contour[0][-1][0]
        finalContour = AddPoint(finalContour, 0, [xNew, yCut, z])
        numAdded = numAdded + 1
    if contour[0][0][1] < yCut and contour[0][-1][1] > yCut:
        if contour[0][0][0] == contour[0][-1][0]: #if xs are same don't need to interpolate
            xNew = contour[0][0][0]
        else:  
            m = (contour[0][-1][1] - contour[0][0][1]) / (contour[0][-1][0] - contour[0][0][0])
            xNew = (yCut - contour[0][0][1])/m + contour[0][0][0]
        finalContour = AddPoint(finalContour, 0, [xNew, yCut, z])    
        numAdded = numAdded + 1
    for i in range(0, len(contour[0])-1):
        if contour[0][i][1] < yCut and contour[0][i+1][1] > yCut:
            if contour[0][i][0] == contour[0][i+1][0]: #if xs are same don't need to interpolate
                xNew = contour[0][i][0]
            else:  
                m = (contour[0][i+1][1] - contour[0][i][1]) / (contour[0][i+1][0] - contour[0][i][0])
                xNew = (yCut - contour[0][i][1]) / m + contour[0][i][0]
            finalContour = AddPoint(finalContour, i + numAdded, [xNew, yCut, z])
            numAdded = numAdded + 1
        elif contour[0][i][1] > yCut and contour[0][i+1][1] < yCut:
            if contour[0][i][0] == contour[0][i+1][0]: #if xs are same don't need to interpolate
                xNew = contour[0][i][0]
            else:    
                m = (contour[0][i+1][1] - contour[0][i][1])/(contour[0][i+1][0] - contour[0][i][0])    
                xNew = (yCut - contour[0][i][1])/m + contour[0][i][0]
            finalContour = AddPoint(finalContour, i + numAdded, [xNew, yCut, z])
            numAdded = numAdded + 1
    return finalContour        

def AddPoint(contour, index, point):
    b = []
    for idx in range(len(contour[0:1])): 
        b.append([])
        if index > 0:
            for j in range(index):
                b[idx].append(contour[idx][j].copy())
        b[idx].append(point)
        if index == len(contour[idx]):
            continue
        for j in range(index, len(contour[idx])):
            b[idx].append(contour[idx][j].copy())
    return b 

def closed_looper(contour):
    if len(contour) == 0:
        return []
    numPoints = len(contour)
    x1 = contour[0][0]
    x2 = contour[numPoints-1][0]
    y1 = contour[0][1]
    y2 = contour[numPoints-1][1]

    if x1 != x2 or y1 != y2:
        contour.append(contour[0].copy())
    return contour    


def x_chop(contours, numCuts):
    xCuts = best_cut_x(contours, numCuts)
 
    for i in range(len(contours)): #First add intersection points
        for j in range(len(xCuts)):
            contours[i] = add_intersections_x(contours[i], xCuts[j])
        #contours[i] = closed_looper(contours[i])

        #now divide into separate parts
    finalContours = []
    divisions = [] #list for each y division for current contour

    #make the list the correct size so there is an item for each x div
    for div in range(len(xCuts)+ 1):
        finalContours.append([])

    for i in range(len(contours)):
        if len(contours[i]) == 0:
            continue
        divisions.clear()
        for x in range(len(xCuts)+1):
            divisions.append([])       
        for x in range(len(xCuts)+1): #a section for each cut, +1

            for j in range(len(contours[i][0])):
                if x == 0:
                    if contours[i][0][j][0] <= xCuts[0]:
                        divisions[x].append(contours[i][0][j].copy())
                elif x == len(xCuts):
                    if contours[i][0][j][0] >= xCuts[x-1]:
                        divisions[x].append(contours[i][0][j].copy())
                else:
                    if contours[i][0][j][0] >= xCuts[x-1] and contours[i][0][j][0] <= xCuts[x]:
                        divisions[x].append(contours[i][0][j].copy())
        
        for div in range(len(divisions)):
            if len(divisions[div]) == 0:
                continue
            temp = closed_looper(divisions[div])
  
            finalContours[div].append(copy.deepcopy([temp]))#[divisions[div]]))
    return finalContours   


def best_cut_x(contours, numCuts):
    errorTolerance = 0.001
    volume = 0
    xCuts= []

    #get total area of contours
    for i in range(len(contours)-1):
        if len(contours[i]) > 0 and len(contours[i+1]) > 0:
            try:
                deltaZ = abs(contours[i+1][0][0][2]-contours[i][0][0][2])
                volume = volume + (GetContourArea(contours[i]) + GetContourArea(contours[i+1])) * deltaZ
            except: pass    

    for cut in range(numCuts):
        maxX = -1000
        minX = 1000
        areaGoal = (cut + 1) / (numCuts + 1)

        for j in range(len(contours)):
            if len(contours[j][0]) > 0:
                for row in range(len(contours[j][0])):
                    if contours[j][0][row][0] > maxX:
                        maxX = contours[j][0][row][0]
                    if contours[j][0][row][0] < minX:
                        minX = contours[j][0][row][0]

        tempContours = []
        cutContours = [] 
        numIters = 0
        error = 1000
        time1 = time.time()
        while error > errorTolerance and numIters < 1000:     
            
            tempContours.clear()
            xCut = (minX + maxX) / 2
            newAreas = []
            newVolume = 0

            for i in range(len(contours)):
                cutContours.clear()
                temp = add_intersections_x(contours[i], xCut)
                #temp = closed_looper(temp)
                if temp == []: 
                    continue
                tempContours.append(temp)
                for j in range(len(tempContours[-1][0])):
                    if tempContours[-1][0][j][0] <= xCut:
                        cutContours.append(copy.deepcopy(tempContours[-1][0][j]))
                if len(cutContours) > 0:
                    cutContours = [cutContours]
                    newAreas.append([GetContourArea(cutContours), cutContours[0][0][2]]) #keep the z value as well for computing volume
                else: 
                    newAreas.append([])    
            for area_idx in range(len(newAreas)-1):    
                if len(newAreas[area_idx]) > 0 and len(newAreas[area_idx+1]) > 0:
                    deltaZ = abs(newAreas[area_idx][1]-newAreas[area_idx+1][1])
                    newVolume = newVolume + (newAreas[area_idx][0] + newAreas[area_idx+1][0]) * deltaZ    
            error = abs((newVolume / volume) - areaGoal)   
            if (newVolume / volume < areaGoal): 
                minX = xCut            
                
            elif (newVolume / volume > areaGoal):
                maxX = xCut    
            #Re-adjust the max/min if stuck
               
            numIters = numIters + 1
            if time.time() - time1 > 5:
                break #stop after 5 seconds    
  
        xCuts.append(xCut) 

    return xCuts

def add_intersections_x(contour, xCut):
    if len(contour) == 0:
        return []
    if len(contour[0]) == 0:
        return []
    finalContour = copy.deepcopy(contour)
    numAdded =1 #start at one, increment after adding each point, to keep track of where to add additional point (add to index)
    for slice in contour:
        try:
            z = slice[0][2]
        except: continue    
    #index 0 outside of loop:
    if contour[0][0][0]  > xCut and contour[0][-1][0] < xCut:
        if contour[0][-1][1] == contour[0][0][1]: #if ys are same don't need to interpolate
            yNew = contour[0][0][1]
        else:        
            m = (contour[0][0][0]- contour[0][-1][0]) / (contour[0][0][1] - contour[0][-1][1])
            yNew = (xCut - contour[0][-1][0]) / m + contour[0][-1][1]
        finalContour = AddPoint(finalContour, 0, [xCut, yNew, z])
        numAdded = numAdded + 1
    if contour[0][0][0] < xCut and contour[0][-1][0] > xCut:
        if contour[0][0][1] == contour[0][-1][1]: #if ys are same don't need to interpolate
            yNew = contour[0][0][1]
        else:        
            m = (contour[0][-1][0] - contour[0][0][0]) / (contour[0][-1][1] - contour[0][0][1])
            yNew = (xCut - contour[0][0][0])/m + contour[0][0][1]
        finalContour = AddPoint(finalContour, 0, [xCut, yNew, z]) 
        numAdded += 1   

    for i in range(0, len(contour[0])-1): #using first island
        if contour[0][i][0] < xCut and contour[0][i+1][0] > xCut:
            if contour[0][i+1][1] == contour[0][i][1]: #if ys are same don't need to interpolate
                yNew = contour[0][i][1]
            else:    
                m = (contour[0][i+1][0] - contour[0][i][0])/(contour[0][i+1][1] - contour[0][i][1])    
                yNew = (xCut - contour[0][i][0])/m + contour[0][i][1]
            finalContour = AddPoint(finalContour, i + numAdded, [xCut, yNew, z])
            numAdded = numAdded + 1
        elif contour[0][i][0] > xCut and contour[0][i+1][0] < xCut:
            if contour[0][i+1][1] == contour[0][i][1]: #if ys are same don't need to interpolate
                yNew = contour[0][i][1]
            else:    
                m = (contour[0][i+1][0] - contour[0][i][0])/(contour[0][i+1][1] - contour[0][i][1])    
                yNew = (xCut - contour[0][i][0])/m + contour[0][i][1]
            finalContour = AddPoint(finalContour, i + numAdded, [xCut, yNew, z])
            numAdded = numAdded + 1
    return finalContour        

def ReOrderParotids(contours, organName, numCuts):
    #Reorder from inferior --> superior, medial --> lateral, anterior --> posterior
    j = 0
    finalContours = []

    if 'l' in organName.lower():
        for i in range(len(contours)):
            if i % (numCuts[2] + 1) == 0 and i > 0:
                j = j + 1
            index = j % ((numCuts[0] + 1) * (numCuts[1] + 1) * (numCuts[2] + 1))     
            finalContours.append(contours[index])
            j = j + (numCuts[0] + 1) * (numCuts[1] + 1)
    else: 
        j = numCuts[0]
        for i in range(len(contours)):
            if j >= (numCuts[0] + 1) * (numCuts[1] + 1) * (numCuts[2] + 1):
                j = j - ( (numCuts[0] + 1) * (numCuts[1] + 1) * (numCuts[2] + 1) +1)
                if j == -1:
                    j = j + (numCuts[0] + 1) * (numCuts[1] + 1)
            finalContours.append(contours[j])
            j = j + (numCuts[0]+1) * (numCuts[1] + 1)
    return finalContours    