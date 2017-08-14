# -*- coding: utf-8 -*-
'''
Created on Fri Jul 14 13:57:16 2017

@version: 1.0

@author: Ryan M. Imrie
@email: ryan.m.imrie@gmail.com

EndosomeToolkit.py is a collection of functions that antomate several steps in 
the processing and analysis of fluorescent microscopy images of endosomes.
In addition, these functions may be broadly applicable to any image where the 
subject of analysis is spherical.

For a detailed description please read "EndosomeToolkitDocumentation.pdf"

'''

# Copyright (C) 2017, Ryan M. Imrie
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# The GNU General Public License can be found at: https://www.gnu.org/licenses/

def fileCleanup(inputDirectory):
    
    '''
    [Usage]
    fileCleanup(inputdirectory)
    
    [Description]
    Tidies up files output from zen batch image converter.
    Will look for z## and c# for z-slice and channel data and retain only these
    in the file name.
    '''

    import os

    for file in os.listdir(inputDirectory):

        if '.tif' in file:
            currentFile = str(file)
            currentFilePath = os.path.join(inputDirectory, file)
            currentFileZ = ''
            currentFileC = ''
            for i in range (0, len(str(currentFile))):
                if currentFile[i] == 'z' and currentFile[i-1] == '_':
                    int1 = ''
                    int2 = ''
                    try:
                        int1 = int(currentFile[i+1])
                        int2 = int(currentFile[i+2])
                    except Exception:
                        pass
                    if type(int1) == int and type(int2) == int:
                        currentFileZ = currentFile[i] + str(int1) + str(int2)
                if currentFile[i] == 'c' and currentFile[i+2] == '_':
                    int3 = ''
                    try:
                        int3 = int(currentFile[i+1])
                    except Exception:
                        pass
                    if type(int3) == int:
                        currentFileC = currentFile[i] + str(int3)
            os.rename(currentFilePath, os.path.join(inputDirectory, currentFileC + currentFileZ + '.tif'))                        

def getMetadata(inputDirectory, inputFile):

    '''
    [Usage]
    getMetadata(inputDirectory, inputFile)
    
    [Description]
    Retrieves the xml metadata saved within .czi files from Zeiss microscopes
    and saves this in a .txt file in the input directory given.    
    '''

    import os

	
	# Print start-up message, combine inputDirectory and inputFile to single path
    print('-----------\nGetMetadata\n-----------\n')
    cziPath = os.path.join(inputDirectory, inputFile)
    print('Getting metadata from:\n' + cziPath + '\n')

	# If file isn't .czi file, print error and exit
    if not '.czi' in inputFile:
        print('Error: File must be a .czi file')
        return

	# Read all lines from czi file. If line isn't text, ignore line. 
    with open(cziPath, 'r', errors='ignore') as f:
        cziFile = f.readlines()

	# Variables for loop
    newFile = []
    metadataFound = False
    currentLine = 1

	# Store all lines between the line containing '<Metadata>' and the line containing '</ImageDocument>'
    for line in cziFile:
		# Print progress update for user every 10000 lines processed
        if currentLine % 10000 == 0 or currentLine == len(cziFile):
            print('\r\tScanning line {} of {}'.format(currentLine, len(cziFile)),
                  end="", flush=True)
        currentLine += 1

        if "<Metadata>" in line:
            metadataFound = True
        elif "</ImageDocument>" in line:
            metadataFound = False
        if metadataFound == True:
            newFile.append(line)
	# Write stored lines to file
    with open (os.path.join(inputDirectory, inputFile[:-4] + 'metadata.txt'), 'w') as f:
        for line in newFile:
            f.write(line)

    print('\n\nComplete\n')
    return;

        
def cellMask(inputDirectory, channels = 2, maskChannel = 1, autoAdjust = False,
             threshold = 100, masterMask = True, fillHoles = True):
    
    '''
    [Usage]
    cellMask(inputDirectory, option1 = , option2 = ,...)
    
    [Description]
    Generates a cell mask from a given channel and applies that mask to every
    channel in the z-stack. Assuming the given channel represents a cell,
    everything outside of that cell will be removed.
    
    [Options]    
    channels = integer (No. of different channels in the z-stack)
    maskChannel = integer (Channel no. to use to generate cell masks)
    autoAdjust = True/False (Set thresholding value to average grey value of z-stack)
    threshold = integer (If autoAdjust = False, acts as threshold value, if autoAdjust = True, is multiplied to average grey value to give threshold value)
    masterMask = True/False (Generate masterMask as max-projection of raw cell masks)
    fillHoles = True/False (Perform morphological closing on cell masks and fill holes)
    '''

    import os
    import glob
    import numpy as np
    from PIL import Image
    from scipy.misc import toimage
    from scipy.ndimage import morphology
    from skimage.color import rgb2grey

	
	# Print start-up message
    print('--------\nCellMask\n--------\n')
    print('Processing files from directory:\n' + inputDirectory + '\n')

	# If inputDirectory doesn't exist, print error and exit
    if not os.path.exists(os.path.join(inputDirectory)):
        print('\nError: Could not find directory')
        return

	# Generate paths for all output directories
    cellMasksDir = os.path.join(inputDirectory, 'cellMasks')
    maskedCellsDir = os.path.join(inputDirectory, 'maskedCells')
    cellMasksRawDir = os.path.join(inputDirectory, 'cellMasksRaw')
    
	# If output directories don't exist, make them	
    for d in [cellMasksDir, maskedCellsDir, cellMasksRawDir]:
        if not os.path.exists(d):
            os.mkdir(d)

	# Calculate number of z-slices by comparing total no. of images to no. of channels
    stackSize = len(glob.glob(os.path.join(inputDirectory, '*.tif')))
    stackSize //= channels

	# If autoAdjust set to True, open and calculate average greyscale value for each image in the mask channel
    if autoAdjust == True:
        averages = []
        for i in range(stackSize):
            print('\r\tCalculating adaptive threshold: {:02d} of {:02d}'
                  .format(i+1, stackSize), end="", flush=True)

            file = 'c{}z{:02d}.tif'.format(maskChannel, i+1)
            image = Image.open(os.path.join(inputDirectory, file))
            averages.append(rgb2grey(np.array(image)))
		# Set the threshold value to threshold * average of average greyscale values
        threshold *= np.average(averages)

    with open(os.path.join(inputDirectory, 'cellMaskLog.txt'), 'w') as f:
        f.write('cellMask Options')
        f.write('maskChannel\t{}\n'.format(maskChannel))
        f.write('autoAdjust\t{}\n'.format(autoAdjust))
        f.write('threshold\t{}\n'.format(threshold))
        f.write('fillHoles\t{}\n'.format(fillHoles))
        
    print('\n\tThreshold Value: {:.3f}'.format(threshold))
    
	# For every z-slice
    for i in range(stackSize):
    
        print('\r\tApplying Threshold: {:02d} of {:02d}'.format(i+1,stackSize), 
              end="", flush=True)
		# Open mask channel image for that z-slice as a greyscale image
        file = 'c{}z{:02d}.tif'.format(maskChannel, i+1)
        image = Image.open(os.path.join(inputDirectory, file))
        imageMask = rgb2grey(np.array(image))
                
        # Find pixels where imageMask is lower/higher than threshold
        lwr = np.where(imageMask <= threshold)
        hgr = np.where(imageMask > threshold)
        # Set the lower/higher pixels to False/True (i.e. 0/1) respectively
        imageMask[lwr] = False
        imageMask[hgr] = True
        imageMask = imageMask.astype(bool)
		
		# Save this 'raw' cell mask as tif
        toimage(imageMask).save(os.path.join(cellMasksRawDir, 'z{:02d}.tif'.format(i+1)))

		# If fillHoles set to true, apply binary closing and fillHoles to raw cell mask
        if fillHoles == True:
            imageMask = morphology.binary_closing(imageMask, structure=np.ones((2,2)))
            imageMask = morphology.binary_fill_holes(imageMask)
			# Save this cell mask
            toimage(imageMask).save(os.path.join(cellMasksDir, 'z{:02d}.tif'.format(i+1)))

		# For each channel in the z-slice
        for n in range(1, channels + 1):
            file = 'c{}z{:02d}.tif'.format(n, i+1)
            # Open the original image
            image = Image.open(os.path.join(inputDirectory, file))
            # Multiply the original image to the cell mask
            new_image = rgb2grey(np.array(image)) * imageMask
            # Save this masked image
            toimage(new_image, cmin=0.0, cmax=65535).save(os.path.join(maskedCellsDir, file))

	# If masterMask set to true, generate maximum z-projection of raw cell masks
    if masterMask == True:
        maxImage = np.array(Image.open(os.path.join(cellMasksRawDir, 'z01.tif')))
        maskDir = os.path.join(inputDirectory, 'masterMasks')
        if not os.path.exists(maskDir):
            os.mkdir(maskDir)  
        print('\r\n', end="", flush=True)
        for i in range(stackSize):
            print('\r\tGenerating Master Mask: {:02d} of {:02d}'.format(i+1,stackSize),
                  end="", flush=True)
            currentImage = np.array(Image.open(os.path.join(
                    cellMasksRawDir, 'z{:02d}.tif'.format(i+1))))
            maxImage = np.maximum(maxImage, currentImage) 
		# Save max projection
        toimage(maxImage).save(os.path.join(inputDirectory, 'masterMasks', 'masterMask.tif'))  

    print('\n\nComplete\n')

    return;          

def masterMask(inputDirectory, maskFile, fillHoles = True):

    '''
    [Usage]
    masterMask(inputDirectory, maskFile, option1)
       
    [Description]
    Applies a master mask of the given filename to the original .tifs within the
    given directory and saves the output in a new subdirectory. Designed as a 
    follow-up to cellMask, and requires the directory and file layout output by
    that function.
    
    [Options]    
    fillHoles = True/False (Perform morphological closing on cell masks and fill holes)
    '''

    import os
    import glob
    import numpy as np
    from PIL import Image
    from scipy.misc import toimage    

	# Print start-up message
    print('----------\nMasterMask\n----------\n')
    print('Processing files from directory:\n' + inputDirectory + '\n')
    print('MasterMask:\n' + maskFile + '\n')

	# If inputDirectory doesn't exist, print error and exit	
    if not os.path.exists(os.path.join(inputDirectory)):
        print('\nError: Could not find directory')
        return
		
	# Generate paths for all directories created by cellMask
    maskDir = os.path.join(inputDirectory, 'maskedCells')    
    cellMaskDir = os.path.join(inputDirectory, 'cellMasks')
    cellMaskRawDir = os.path.join(inputDirectory, 'cellMasksRaw')

	# If subdirectories for master mask output don't exist, make them	
    for d in [maskDir, cellMaskDir,cellMaskRawDir]:
        if not os.path.exists(os.path.join(d, maskFile[:-4])):
            os.mkdir(os.path.join(d, maskFile[:-4]))   			

	# Calculate number of z-slices by comparing total no. of mask images to total no. of original images
    stackSize = len(glob.glob(os.path.join(cellMaskRawDir, '*.tif')))        
    tifs = len(glob.glob(os.path.join(inputDirectory, '*.tif')))
    channels = int(tifs / stackSize)

	# Open master mask image as a binary image
    masterMask = np.array(Image.open(os.path.join(inputDirectory, 'masterMasks', maskFile)))
    masterMask = masterMask.astype(bool)		
	
	# For every z-slice	
    for i in range(stackSize):
        print('\r\tApplying masterMask: {:02d} of {:02d}'.format(i+1,stackSize),
              end="", flush=True)
        currentFile = 'z{:02d}.tif'.format(i+1)
        # Generate subpath for output directory and current file
        outDir = os.path.join(maskFile[:-4], currentFile)

		# Open mask for current z-slice as a binary image
        currentMask = np.array(Image.open(os.path.join(cellMaskRawDir, currentFile)))
        currentMask = currentMask.astype(bool)
		# Multiply current z-slice mask to master mask
        newMask = currentMask * masterMask

		# Save this edited 'raw' mask
        toimage(newMask).save(os.path.join(cellMaskRawDir, outDir))              

		# If fillHoles set to True
        if fillHoles == True:
			# Open processed mask for current z-slice
            currentMask = np.array(Image.open(os.path.join(cellMaskDir, currentFile)))
            currentMask = currentMask.astype(bool)
			# Multiply current z-slice processed mask to master mask
            newMask = currentMask * masterMask

			# Save this edited mask
            toimage(newMask).save(os.path.join(cellMaskDir, outDir))              

		# For each channel in the z-slice
        for n in range(channels):
            currentFile = 'c{}z{:02d}.tif'.format(n+1, i+1)
            outDir = os.path.join(maskFile[:-4], currentFile)
            # Open the original image
            image = np.array(Image.open(os.path.join(inputDirectory, currentFile)))
            # Multiply the original image to the new cell mask
            newImage = image * newMask
            # Save this masked image
            toimage(newImage).save(os.path.join(maskDir, outDir))

    print('\n\nComplete\n')


    return;

def isolateEndosomes(inputDirectory, imagePath, channel = 1, threshold = 35, autoAdjust = True,
                     scaling = 4, minDistance = 15, minSize = 1.5, maxSize = 10, useMetadata = True): 

    import os
    import glob
    import math
    import numpy as np
    from PIL import Image
    from scipy import ndimage
    from scipy.misc import toimage, imresize    
    from skimage.feature import peak_local_max
    from skimage.morphology import watershed, remove_small_objects

	# Print start-up message, combine inputDirectory and imagePath to single path
    print('----------------\nIsolateEndosomes\n----------------\n')
    fileDir = os.path.join(inputDirectory, imagePath)
    print('Processing files from directory:\n' + fileDir + '\n')
    print('Processing channel:\nc{}\n'.format(channel))

	# Calculate number of z-slices by counting number of images in a single channel
    stackSize = 0
    for file in os.listdir(fileDir):
        if 'c{}'.format(channel) in file and '.tif' in file:
            stackSize += 1

	# Generate path for output directory
    outDir = os.path.join(fileDir, 'isolatedEndosomes')
	# If output directory doesn't exist, make it
    if not os.path.exists(outDir):
        os.mkdir(outDir)

	# if useMetadata set to True
    if useMetadata == True:
        dataDir = fileDir

        dataDir = os.path.normpath(fileDir)
        dataDir = dataDir.split(os.sep)
        for i in range(0, len(dataDir)):
            if 'maskedCells' in dataDir[i]:
                metaDataName = dataDir[i-1]
        
		# Open the metadata file generated by getMetadata and collect 2D scaling value
        with open(os.path.join(inputDirectory, metaDataName + 'metadata.txt'), 'r') as f:
            for line in f:
                if 'ScalingX' in line:
                    microscopeScaling = float(''.join(line.split())[10:-11])
					# Convert from metre to micrometre
                    microscopeScaling *= 1000000
                if 'ScalingZ' in line:
                    microscopeScalingZ = float(''.join(line.split())[10:-11])
					# Convert from metre to micrometre
                    microscopeScalingZ *= 1000000                 
        voxScaling = microscopeScalingZ / microscopeScaling

    else:
		# Essentially, ignore microscope scaling
        microscopeScaling = 1
        voxScaling = 1

	# Create log file for isolateEndosome run
    with open(os.path.join(outDir, 'isolateEndosomesLog.txt'), 'w') as f:
        f.write('Option Values:\n')
        f.write('autoAdjust\t{}\n'.format(autoAdjust))
        f.write('threshold\t{}\n'.format(threshold))
        f.write('scaling\t{}\n'.format(scaling))
        f.write('minSize\t{}\n'.format(minSize))
        f.write('maxSize\t{}\n'.format(maxSize))
        f.write('minDistance\t{}\n'.format(minDistance))

	# Apply microscopeScaling and scaling to minSize, maxSize, and minDistance
	# This converts the micrometre value of minSize to pixel units
    
    minSize = ((minSize**0.5) / microscopeScaling)**2
    minSize = int(minSize * (scaling**2))
    maxSize = ((maxSize**0.5) / microscopeScaling)**2
    maxSize = int(maxSize * (scaling**2))
    minDistance = minDistance / microscopeScaling
    minDistance = int(minDistance * scaling)
 
    # Use area and volume formulas, along with voxel scaling to calculate total min and max voxels
    # To be used for imageJ 3D Object Counter
    minRad = int((minSize / math.pi)**0.5)
    minVol = int(((4/3) * (math.pi * (minRad**3))) / (voxScaling))
    maxRad = int((maxSize / math.pi)**0.5)
    maxVol = int(((4/3) * (math.pi * (maxRad**3))) / (voxScaling))

    # Appends size conversions and total voxel calculations to log file    
    with open(os.path.join(outDir, 'isolateEndosomesLog.txt'), 'a') as f:
        f.write('\nSize conversions:\n')
        f.write('minSize (Pixels)\t{}\n'.format(minSize))
        f.write('maxSize (Pixels)\t{}\n'.format(maxSize))
        f.write('minDistance (Pixels)\t{}\n'.format(minDistance))
        f.write('Microscope Scaling (micrometre/pixel)\t{}\n'.format(microscopeScaling))
        f.write('Voxel Scaling (XY:Z)\t{}\n'.format(voxScaling))
        
        f.write('\nFor imageJ 3D Object Counter:\n')        
        f.write('Min Volume (Voxels)\t{}\n'.format(minVol))
        f.write('Max Volume (Voxels)\t{}\n'.format(maxVol))
        
	# For every z-slice	
    for i in range(stackSize):
        currentThreshold = threshold
        currentFile = 'c{}z{:02d}.tif'.format(channel, i+1)

        # Open image        
        image = np.array(Image.open(os.path.join(fileDir, currentFile)))

        # If autoAdjust set to True, add average grey value of image to threshold
        if autoAdjust == True:
            average = np.mean(image)
            currentThreshold += average

        # Set all values in image below threshold to 0
        image[image < currentThreshold] = 0

        # Resize image by value of scaling
        image = imresize(image, scaling * 100)

        print('\r\tApplying watershed transform: {:02d} of {:02d} (Threshold: {:.3f})'.format(i+1,
              stackSize, currentThreshold),end="", flush=True)

        # If the image is blank, skip to next image        
        if np.mean(image) == 0:
            newImage = image

        else:               
            # Distance transform
            distance = ndimage.distance_transform_edt(image)

            # Detect peak local maxima 
            localMaxi = peak_local_max(distance, indices=False, 
                                               min_distance=minDistance)

            # Enlarge local maxima
            localMaxi = ndimage.filters.maximum_filter(localMaxi, size = (2,2))
            markers, objectCount = ndimage.label(localMaxi)

            # For each enlarged peak, calculate the centre of the peak        
            objectRange = []
            for i in range(0,objectCount):
                objectRange.append(i)   
            localMaxiCoords = ndimage.measurements.center_of_mass(localMaxi, labels = markers, index = objectRange)       
            localMaxiCoords = [list(tuples) for tuples in localMaxiCoords]

            # Split coordinates for each peak centre into x and y
            localMaxiCoordsx = []
            localMaxiCoordsy = []
            for i in localMaxiCoords:
                try:
                    localMaxiCoordsx.append(int(i[0]))
                    localMaxiCoordsy.append(int(i[1]))
                except:
                    pass                

            # Clear original image of local maxima        
            localMaxi = np.zeros((np.shape(localMaxi)))

            # Place a white pixel at the coordinates for each peak centre
            localMaxi[localMaxiCoordsx, localMaxiCoordsy] = 1

            # Apply a unique label to each peak
            markers, objectCount = ndimage.label(localMaxi)

            # Watershed transform        
            labels = watershed(-distance, markers, mask=image)

            # Detect edges between differently labelled objects
            sx = ndimage.sobel(labels, axis=0, mode='constant')
            sy = ndimage.sobel(labels, axis=1, mode='constant')
            sob = np.hypot(sx, sy)
            # Make every edge white
            sob[sob > 0] = 1
            #Invert the edge image
            sob = ~sob.astype(bool)

            # Set all values in image that remain to 1
            image[image > 0] = 1
            image= image.astype(bool)
            # Multiply edge image to image
            newImage = sob * image

            # Remove objects of total pixel count < minSize        
            newImage = remove_small_objects(newImage, minSize)

        # Save segmented image
        toimage(newImage).save(os.path.join(outDir, currentFile))

    return;

'''
    
def euclidDistance(inputDirectory, objectFile, centreFile):
    
    import os
    import numpy as np
    
    print('--------------\n3dLocalisation\n--------------\n')
    print('Processing files from directory:\n' + inputDirectory + '\n')

    
    with open (os.path.join(directory, centreFile), 'r') as f:
        centreVolume = 0.0
        for line in f:
            if float(line.rsplit('\t')[]) > centreVolume:
                centre = np.array((line.rsplit('\t')[],
                                   line.rsplit('\t')[],
                                   line.rsplit('\t')[]))
               
        
    with open (os.path.join(directory, objectFile), 'r') as f:
        localisation = []
        for line in f:
            currentObject = np.array((line.rsplit('\t')[],
                                      line.rsplit('\t')[],
                                      line.rsplit('\t')[]))
            
            localisation.append[np.linalg.norm(currentObject-centre)]
            
    with open(os.path.join(directory, '3dLocalisation.txt'), 'w') as f:
        for i in range (0, len(localisation)):
            f.write('{}\t{}'.format(i+1, localisation[i]))
            
    
    return;

'''