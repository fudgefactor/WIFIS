import pyfits
import numpy as np
import os
import matplotlib.pyplot as plt

saturation_cutoff = 30000

class hxrg_frame:
    def __init__(self,filename):
        self.filename = filename

        hdulist = pyfits.open(filename)

        # Detector/ASIC properties
        # MUX Type, 1->H1RG, 2->H2RG (implemented), 3->H4RG
        self.mux_type = int(hdulist[0].header['MUXTYPE'])
        if(self.mux_type == 2):
            self.detstr = 'H2RG'
            self.nx = 2048
            self.ny = 2048
        #ASIC ID
        self.asic_id = hdulist[0].header['ASIC_NUM']   
        #SCA ID
        self.sca_id = hdulist[0].header['SCA_ID']
        # Number of readout channels
        self.noutputs = int(hdulist[0].header['NOUTPUTS'])
        # Gain setting of ASIC
        self.asic_gain = int(hdulist[0].header['ASICGAIN'])

        # Observation properties
        # UTC Julian time of exposure
        self.jacqtime = float(hdulist[0].header['ACQTIME'])
        # Frame time for exposure
        self.frametime = float(hdulist[0].header['FRMTIME'])
        # Exposure time in ramp
        self.seqnum_m = float(hdulist[0].header['SEQNUM_M'])
        if self.seqnum_m == 0:
            self.exptime = float(0)
        else:
            self.exptime = float(hdulist[0].header['INTTIME'])
        # Image data
        self.imgdata = hdulist[0].data                           
        self.imgdata_arr = np.array(self.imgdata)
        hdulist.close()

        # Has the image been reference corrected
        self.corrected = False

    # Rudimentary reference pixel correction using the top and bottom rows done on a per channel basis
    def ref_correct_frame(self):
        if (self.noutputs == 32):
            
            for i in range(32):
                avgdrift = 0.5*(np.average(self.imgdata[0:4,i*64:i*64+64])+np.average(self.imgdata[2044:2048,i*64:i*64+64]))
                self.imgdata[:,i*64:i*64+64] -= avgdrift

            self.corrected = True

    #returns the values of reference pixels
    def get_ref_pixels(self):
        self.h_reference_pixels_1 = self.imgdata[:, 0:3]
        self.h_reference_pixels_2 = self.imgdata[:, self.nx-4:self.nx-1]
        self.v_reference_pixels_1 = self.imgdata[0:3, :]
        self.v_reference_pixels_2 = self.imgdata[self.ny-4:self.ny-1, :]
        

class h2rg_ramp:
    def __init__(self, folder, ramp_parameters):
        self.folder = folder
        #a list containg ramp parameters, user input
        # idx = 0: nResets
        # idx = 1: nGroups
        # idx = 2: nFrames
        # idx = 3: nDrops
        # idx = 4: nRamps
        #will use the above parameters to generate filenames of the fits files
        self.ramp_parameters = ramp_parameters
        self.currentDir = os.getcwd()
        os.chdir(folder)
        
    def __exit__(self):
        os.chdir(self.currentDir)
    
    #reads the bad pixel mask fits. badPixelMask_FN should also include the path to the file
    def take_badPixMask(self, badPixelMask_FN):
        self.badPixMask = hxrg_frame(badPixelMask_FN)
    
    #Generage list of filenames in the ramp based on the ramp parameters
    def generateFileList(self):
        
        filenames = []
        nResets = self.ramp_parameters[0]
        nGroups = self.ramp_parameters[1]
        nFrames = self.ramp_parameters[2]
        nRamps = self.ramp_parameters[4]
        for ramp in range(0, nRamps):
            rampStr = "%02d" % (ramp +1)
            for reset in range(0, nResets):
                resetStr = "%02d" % (reset+1)
                resetFN = 'H2RG_R' + rampStr + '_M00_N' + resetStr + '.fits'
                filenames.append(resetFN)
            for group in range(0,  nGroups):
                groupStr = "%02d" % (group +1)
                for frame in range(0,  nFrames):
                    frameStr = "%02d" % (frame+1)
                    frameFN = 'H2RG_R' + rampStr + '_M' + groupStr + '_N' + frameStr + '.fits'
                    filenames.append(frameFN)
        
        return filenames
    
    #read the ramp
    def read_ramp(self):
        fileList = self.generateFileList()
        
        frameList = []
        
        for file in fileList:
            frameList.append(hxrg_frame(file))
            
        self.frameList = frameList
        
    #Fit slopes for each pixel
    def fit_slopes(self):
        
        pixelValue = []
        expTime = []
        
        print 'Reading fits files ...'
        
        for image in self.frameList:
            if image.seqnum_m == 0:
                continue
            pixelValue.append(np.array(image.imgdata))
            expTime.append(image.exptime)
            print image.exptime
        expTime = np.array(expTime)
        
        print 'done'
        
        #validPixels = self.badPixMask  !!!Will add later
        validPixels = image.imgdata*0
        validPixels[5:500, 5:500] = 1
        
        validPixIdx_temp = np.where(validPixels == 1)
        
        validPixX = validPixIdx_temp[0]
        validPixY = validPixIdx_temp[1]

        fittingData = []
        fittedParameters = []
        
        validPixIdx = zip(validPixX, validPixY)
        
        print 'Generating data cube ...'
        pixelArray = np.array(pixelValue) #VERY RESOURCE HEAVY!!!!!!!!!!!!!!!!!
        print 'done'
        
        print 'Running slope fit ...'
        
        for pixel in validPixIdx:
            idxX = pixel[0]
            idxY = pixel[1]
            columnData = pixelArray[:, idxX, idxY]
            validColumnDataIdx = np.where(columnData <= saturation_cutoff)

            validColumnData = list(columnData[i] for i in validColumnDataIdx)
            validExpTime = list(expTime[i] for i in validColumnDataIdx)

            fittedParameters.append(np.polyfit(validExpTime[0],  validColumnData[0], 1))
        print 'done'
        
