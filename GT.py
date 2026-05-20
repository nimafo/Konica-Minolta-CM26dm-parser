
"""
File: GT.py
Author: Nima Forouzandeh
Date: 2024-10-01
Version: 0.0
"""
import pandas as pd
import numpy as np
import colour
import ast


class GT:
    """
        Description:
            This script provides Python functions to process ground truth data from reflectance spectrophotometer CM-26dg.
            
            Inputs:
            - input_file: a CSV file containing reflectance data from the spectrophotometer.
            - n: material ID to process.
            - delimiter: delimiter used in the input file. Default is ','. for tab-separated files, use '\t'.
            
            Functions:
            __calculate_rgb: calculates the RGB values and visible reflectance for a given material.
            
            Attributes:
            - input_file: a CSV file containing reflectance data from the spectrophotometer.
            - dataset: pandas dataframe containing the reflectance data.
            - material_count: the number of materials in the input file.
            - material: the reflectance data for the nth material and L*(D65), a*(D65), b*(D65), Target Name
            - SCI: the reflectance data for the SCI target.[%]
            - SCE: the reflectance data for the SCE target.[%]
            - spectral_SCI: the spectral reflectance for the SCI target.[%]
            - spectral_SCE: the spectral reflectance for the SCE target.[%]
            - Lab_SCI: the CIE L*a*b* values for the SCI target with D65 illuminant.
            - Lab_SCE: the CIE L*a*b* values for the SCE target with D65 illuminant.
            - RGB_SCI: the RGB values for the SCI target.[0-1]
            - RGB_SCE: the RGB values for the SCE target.[0-1]
            - reflectance_SCI: the visible reflectance for the SCI target.[0-1]
            - reflectance_SCE: the visible reflectance for the SCE target.[0-1]
            - specularity: the specularity of the material.[0-1]
            """
    def __init__(self, input_file,n=None, delimiter = ',', decimal = '.'):
        ''' Dataset Level'''
        self.dataset = pd.read_csv(input_file, index_col = False, delimiter = delimiter, header=0, decimal=decimal)
        self.material_count = int(len(self.dataset.iloc[0:,0:0])/2)
        self.wavelength_range = range(360,750,10)
        if n:
            if n > self.material_count: # n is 
                raise ValueError(f"Material {n} does not exist. There are only {self.material_count} materials.")
        '''Material Level'''
        nth_material = self.dataset.iloc[2*n-2:2*n, :]
        columns_to_select = ['Data Name']+['L*(D65)']+['a*(D65)']+['b*(D65)']+['Target Name'] + [f'{i}nm' for i in self.wavelength_range]
        self.material = nth_material.loc[:, columns_to_select]
        self.SCI = self.material[self.material['Target Name'] == 'SCI'].iloc[0]
        self.SCE = self.material[self.material['Target Name'] == 'SCE'].iloc[0]
        self.spectral_SCI = self.SCI[[f'{i}nm' for i in self.wavelength_range]]
        self.spectral_SCI.index = [int(i[:-2]) for i in self.spectral_SCI.index] # remove the nm from the column names
        self.spectral_SCE = self.SCE[[f'{i}nm' for i in self.wavelength_range]]
        self.spectral_SCE.index = [int(i[:-2]) for i in self.spectral_SCE.index] # remove the nm from the column names
        self.Lab_SCI = self.SCI[['L*(D65)', 'a*(D65)', 'b*(D65)']]
        self.Lab_SCE = self.SCE[['L*(D65)', 'a*(D65)', 'b*(D65)']]
        self.RGB_SCI = self._lab_to_rgb(self.Lab_SCI)  
        self.RGB_SCE = self._lab_to_rgb(self.Lab_SCE)
        self.reflectance_SCI = self._rgb_to_vis(self.RGB_SCI)
        self.reflectance_SCE = self._rgb_to_vis(self.RGB_SCE)
        self.specularity = self._get_specularity()
    def _lab_to_rgb(self, Lab):
        """
            Convert from CIE L*a*b* (D65) to RGB (sRGB). Illuminant D65 is assumed.
            
            Parameters:
                Lab (array-like): L*, a*, b* values in the CIE L*a*b* color space under D65 illuminant.
            
            Returns:
                RGB (array): Corresponding RGB values, clipped to the [0, 1] range.
            """
        D65_white = colour.CCS_ILLUMINANTS['CIE 1931 2 Degree Standard Observer']['D65']
        XYZ = colour.Lab_to_XYZ(Lab, illuminant=D65_white)#Lab to XYZ
        RGB = colour.XYZ_to_RGB( # XYZ to RGB (sRGB)
            XYZ,
            illuminant_XYZ=D65_white,
            illuminant_RGB=D65_white,
            matrix_XYZ_to_RGB=colour.RGB_COLOURSPACES['sRGB'].matrix_XYZ_to_RGB
        )

        if np.any(RGB > 1):
            print(rf"Warning: RGB values are out of range. values are clipped to [0, 1]. actual values: {RGB}")
            RGB = np.clip(RGB, 0, 1)
        
        return RGB
    def _rgb_to_vis(self, RGB):
        """
            Convert from RGB to visible reflectance.
            
            Parameters:
                RGB (array-like): RGB values in the sRGB color space.
            
            Returns:
                vis (float): Visible reflectance.
            """
        vis = 0.265*RGB[0] + 0.670*RGB[1] + 0.065*RGB[2]
        #print(rf"Visible reflectance: {vis} from RGB: {RGB}")
        if vis > 1:
            print("Warning: Visible reflectance is out of range. values are clipped to [0, 1].")
            vis = np.clip(vis, 0, 1)
        return vis
    def _get_specularity(self):
        """
            Calculate the specularity of the material.
            
            Returns:
                specularity (float): Specularity of the material.
            """
        if self.reflectance_SCI > self.reflectance_SCE:
            return self.reflectance_SCI - self.reflectance_SCE
        else:
            print("SCI reflectance is less than SCE reflectance. Specularity is 0.")
            return 0




class GT_fromSpectralDB:
    """
    Process a spectral material database in CSV format and extract properties for a specific material by row number.
    Compatible with SpectralDB-style structure, providing outputs similar to the original GT class.
    """
    def __init__(self, input_file, n=None, verbose=False):
        self.dataset = pd.read_csv(input_file)
        self.material_count = len(self.dataset)
        self.wavelength_range = list(range(360, 750, 10))

        if n is None or n < 1 or n > self.material_count:
            raise ValueError(f"Material index {n} is invalid. Choose between 1 and {self.material_count}.")

        row = self.dataset.iloc[n - 1]
        self.material = row
        self.SCI = self._parse_spectrum(row['SCIMeasures'])
        self.SCE = self._parse_spectrum(row['SCEMeasures'])

        self.spectral_SCI = pd.Series(self.SCI)
        self.spectral_SCE = pd.Series(self.SCE)

        self.Lab_SCI = np.array([row['L'], row['a'], row['b']])
        self.Lab_SCE = np.array([row['L'], row['a'], row['b']])  # assuming same Lab for SCE if not given separately

        self.RGB_SCI = self._lab_to_rgb(self.Lab_SCI)
        self.RGB_SCE = self._lab_to_rgb(self.Lab_SCE)

        self.reflectance_SCI = self._rgb_to_vis(self.RGB_SCI)
        self.reflectance_SCE = self._rgb_to_vis(self.RGB_SCE)

        self.specularity = self._get_specularity()
        self.verbose = verbose

    def _parse_spectrum(self, spectrum_str):
        try:
            parsed = ast.literal_eval(spectrum_str)
            return {int(k): float(v) for k, v in parsed.items()}
        except:
            return {}

    def _lab_to_rgb(self, Lab):
        D65_white = colour.CCS_ILLUMINANTS['CIE 1931 2 Degree Standard Observer']['D65']
        XYZ = colour.Lab_to_XYZ(Lab, illuminant=D65_white)
        RGB = colour.XYZ_to_RGB(
            XYZ,
            illuminant_XYZ=D65_white,
            illuminant_RGB=D65_white,
            matrix_XYZ_to_RGB=colour.RGB_COLOURSPACES['sRGB'].matrix_XYZ_to_RGB
        )
        if np.any(RGB > 1):
            print(rf"Warning: RGB values are out of range. Clipping applied. Original: {RGB}")
            RGB = np.clip(RGB, 0, 1)
        return RGB

    def _rgb_to_vis(self, RGB):
        vis = 0.265 * RGB[0] + 0.670 * RGB[1] + 0.065 * RGB[2]
        if vis > 1:
            print("Warning: Visible reflectance is out of range. Clipping applied.")
            vis = np.clip(vis, 0, 1)
        return vis

    def _get_specularity(self):
        if self.reflectance_SCI > self.reflectance_SCE:
            return self.reflectance_SCI - self.reflectance_SCE
        else:
            print("SCI reflectance is less than or equal to SCE reflectance. Specularity set to 0.")
            return 0.0
