import pyfits
import numpy as np

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
        self.exptime = float(hdulist[0].header['EXPTIME'])
        # Image data
        self.imgdata = hdulist[0].data                           
        
        hdulist.close()

        # Has the image been reference corrected
        self.corrected = False

    # Rudimentary reference pixel correction using the top and bottom rows done on a per channel basis
    def ref_correct_frame(self):
        if (self.noutputs == 32):
            
            for i in range(32):
                avgdrift = 0.5*(np.average(self.imgdata[0:4,i*64:i*64+64])+
                                np.average(self.imgdata[2044:2048,i*64:i*64+64]))
                self.imgdata[:,i*64:i*64+64] -= avgdrift

            self.corrected = True

                
            

