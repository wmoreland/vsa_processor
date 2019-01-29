#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################################################
#    vsa_processor: calculates vesicle number densities from areas
#    Copyright (C) 2019  William M Moreland
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    http://www.gnu.org/licenses/gpl-3.0.html
############################################################################

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def calculate_NAi(input_areas, input_info):
    """
    Takes measured areas and related meta data and calculates the number of
    vesicles of size i per unit area (NAi). This is done separately for each
    magnification.
    """

    # Calculate reference area and vesicle and phenocryst fractions of images
    info = input_info.copy()
    info.index = info['image']
    info['ref_area'] = ((info['width'] * info['height'])
                        / info['scale_factor'] ** 2
                        - ((info['width'] * info['height'])
                           / info['scale_factor'] ** 2
                           * (info['edge_grey'] / 255)))
    info['ves_frac'] = info['vesicles_grey'] / 255
    info['phen_frac'] = info['crystals_grey'] / 255

    # The area input contains several columns per magnification which must be
    # merged.
    col_scan, col_25, col_100, col_250 = [], [], [], []
    for column in input_areas.columns:
        if re.search(r'scan', column):
            col_scan.append(re.search(r'scan', column).group())
        if re.search(r'billet+', column):
            col_scan.append(re.search(r'billet', column).group())
        if re.search(r'25[a-z]+', column):
            col_25.append(re.search(r'25[a-z]+', column).group())
        if re.search(r'100[a-z]+', column):
            col_100.append(re.search(r'100[a-z]+', column).group())
        if re.search(r'250[a-z]+', column):
            col_250.append(re.search(r'250[a-z]+', column).group())

    areas = pd.DataFrame({'scan': input_areas['scan'],
                          '25': input_areas[col_25].stack().reset_index()[0],
                          '100': input_areas[col_100].stack().reset_index()[0],
                          '250': input_areas[col_250].stack().reset_index()[0]}
                         )

    # Calculate diameters from areas
    diam = pd.DataFrame({'scan': (4/np.pi * areas['scan']) ** 0.5,
                         '25': (4/np.pi * areas['25']) ** 0.5,
                         '100': (4/np.pi * areas['100']) ** 0.5,
                         '250': (4/np.pi * areas['250']) ** 0.5})

    # The bins are calculated based on a minimum circular vesicle size which
    # can be accurately represented by square pixels. In this case the minimum
    # is 11 pixels or 7 µm at ×250 magnification. Geometric bins are used
    # rather than linear as the vesicles cover several orders of magnitude in
    # size.
    min_size_px = 11
    scale_factor_250 = 1693.33
    max_size_mm = 5
    bins = [min_size_px/scale_factor_250]
    geobin = bins[0]
    while geobin < max_size_mm:
        geobin *= 10 ** 0.1
        bins.append(geobin)

    # Create histograms from diameters and bins
    freq = pd.DataFrame({'scan': np.histogram(diam['scan'].dropna(), bins)[0],
                         '25': np.histogram(diam['25'].dropna(), bins)[0],
                         '100': np.histogram(diam['100'].dropna(), bins)[0],
                         '250': np.histogram(diam['250'].dropna(), bins)[0],
                         'bins': bins[:-1]})
    freq = freq[['bins', 'scan', '25', '100', '250']]

    # Calculate number of vesicles of size i per unit area (NAi)
    NAi = pd.DataFrame({'scan': (freq['scan']
                                 / info.loc['scan', 'ref_area'].sum()),
                        '25': (freq['25']
                               / info.loc[col_25, 'ref_area'].sum()),
                        '100': (freq['100']
                                / info.loc[col_100, 'ref_area'].sum()),
                        '250': (freq['250']
                                / info.loc[col_250, 'ref_area'].sum()),
                        'bins': (bins[:-1])})
    NAi = NAi[['bins', 'scan', '25', '100', '250']]
    return NAi, info, areas, diam, freq


def merge_NAi(NAi, info, sample, plot=False):
    """
    Merge NAi data to create smooth transition between magnifications using the
    minimized ΔNA method of Shea et al. (2010)
    """
    NAi_merged, sutures, source = [], [], []

    switch_250 = True
    switch_100 = True
    switch_25 = True
    for i in range(len(NAi)):

        # Initialise unified Na list with first value from 250 list
        if not NAi_merged:
            NAi_merged.append(NAi['250'][i])
            source.append('250')
            continue

        # Check that using the next 250 value would give a smaller delta
        if switch_250:
            if NAi_merged[-1] - NAi['250'][i] < NAi_merged[-1] - NAi['100'][i]:
                NAi_merged.append(NAi['250'][i])
                source.append('250')
                continue
            else:
                sutures.append([NAi['bins'][i], NAi['250'][i]])
                switch_250 = False

        # Check that using the next 100 value would give a smaller delta
        if switch_100:
            if NAi_merged[-1] - NAi['100'][i] < NAi_merged[-1] - NAi['25'][i]:
                NAi_merged.append(NAi['100'][i])
                source.append('100')
                continue
            else:
                sutures.append([NAi['bins'][i], NAi['100'][i]])
                switch_100 = False

        # Check that using the next 25 value would give a smaller delta
        if switch_25:
            if NAi_merged[-1] - NAi['25'][i] < NAi_merged[-1] - NAi['scan'][i]:
                NAi_merged.append(NAi['25'][i])
                source.append('25')
                continue
            else:
                sutures.append([NAi['bins'][i], NAi['25'][i]])
                switch_25 = False

        # Fill the remainder of the list with scan values
        NAi_merged.append(NAi['scan'][i])
        source.append('scan')
    NAi['NAi'] = NAi_merged
    NAi['source'] = source

    # If present, correct for phenocrysts
    if info['phen_frac'].mean() != 0:

        # Group the phenocryst fraction values by image magnification
        NAi_adj, phen_scan, phen_25, phen_100, phen_250 = [], [], [], [], []
        for i in info.index:
            if re.search(r'scan', i):
                phen_scan.append(info.loc[i, 'phen_frac'])
            if re.search(r'billet', i):
                phen_scan.append(info.loc[i, 'phen_frac'])
            if re.search(r'25[a-z]+', i):
                phen_25.append(info.loc[i, 'phen_frac'])
            if re.search(r'100[a-z]+', i):
                phen_100.append(info.loc[i, 'phen_frac'])
            if re.search(r'250[a-z]+', i):
                phen_250.append(info.loc[i, 'phen_frac'])

        # Adjust NAi using the relevant adjustment factor
        for i in NAi.index:
            if NAi.loc[i, 'source'] == 'scan':
                NAi_adj.append(
                        NAi.loc[i, 'NAi']/(1-np.array(phen_scan).mean()))
            if NAi.loc[i, 'source'] == '25':
                NAi_adj.append(NAi.loc[i, 'NAi']/(1-np.array(phen_25).mean()))
            if NAi.loc[i, 'source'] == '100':
                NAi_adj.append(NAi.loc[i, 'NAi']/(1-np.array(phen_100).mean()))
            if NAi.loc[i, 'source'] == '250':
                NAi_adj.append(NAi.loc[i, 'NAi']/(1-np.array(phen_250).mean()))
        NAi['NAi adj.'] = NAi_adj

    # Plot the various NAi data to check smoothness of merged set
    if plot:
        fig, ax = plt.subplots()
        for mag, c in zip(['250', '100', '25', 'scan'], ['r', 'g', 'b', 'y']):
            ax.plot(NAi['bins'], NAi[mag], label=mag, color=c)
            ax.plot(NAi['bins'], NAi['NAi'], color='k', linestyle='--',
                    label='Merged')
            ax.set_xscale('log')
            ax.set_title(sample)
            ax.set_xlabel('Vesicle diameter (mm)')
            ax.set_ylabel('NA$_i$ (mm$^{-1}$)')
            ax.legend()

    return NAi


def calculate_NVi(NAi, info, clast_ves):
    """
    Uses the method of Sahagian and Proussevitch (1998) to convert vesicle area
    to volume and then calculates the number of vesicles per unit volume and
    the resulting vesicle number densities.
    """
    # Remove empty bins
    for value in NAi['NAi']:
        if NAi['NAi'].iloc[-1] == 0:
            NAi = NAi.drop(NAi.index[-1])

    # Make new dataframe. Needs new bins to fill the 32 rows needed for the
    # calculated alphas below
    NVi = pd.DataFrame()
    geobins = [NAi['bins'].max()]
    for i in range(31):
        newbin = geobins[-1] / 10 ** 0.1
        geobins.append(newbin)
    NVi['bins'] = geobins[::-1]

    # Insert NAi so that the values for the largest bins are aligned to bottom
    if info['phen_frac'].mean() != 0:
        NAi_temp = list(NAi['NAi adj.'])
    else:
        NAi_temp = list(NAi['NAi'])
    for i in range(len(NVi) - len(NAi)):
        NAi_temp.insert(0, np.nan)
    NVi['NAi'] = NAi_temp
    NVi['NAi'].fillna(0, inplace=True)

    # Create volumetric bin and class columns. Set class as index
    NVi['VolBins'] = NVi['bins']**3
    NVi['Class'] = (np.arange(1, 33))[::-1]
    NVi.set_index('Class', drop=True, inplace=True)

    # Conversion coefficients (alpha, Sahagian & Proussevitch 1998)
    alpha = [1.47610582647948E-10, 2.90280443748631E-10, 5.70843214808730E-10,
             1.12257968617433E-09, 2.20758824251399E-09, 4.34130053222808E-09,
             8.53735359658935E-09, 1.67891682811450E-08, 3.30171001740199E-08,
             6.49313818862787E-08, 1.27696678119265E-07, 2.51141733322558E-07,
             4.93947059438401E-07, 9.71575138306916E-07, 1.91128675756108E-06,
             3.76060963066430E-06, 7.40149346605457E-06, 1.45740496129929E-05,
             2.87178157121063E-05, 5.66506478435492E-05, 1.11946135410306E-04,
             2.21811628924847E-04, 4.41358745081414E-04, 8.84056423472059E-04,
             1.78950062705712E-03, 3.68381080582122E-03, 7.79463714537751E-03,
             1.72711049259740E-02, 4.14945114741756E-02, 1.16190477752079E-01,
             4.56122645308450E-01, 1.64612085334339E+00]
    NVi['Alphas'] = alpha

    # Calculate Hbar values
    Hbar = []
    for i in NVi.index:
        hbar = ((NVi['VolBins'][i] + NVi['VolBins'][i]*10**0.3) / 2)**(1/3)
        Hbar.append(hbar)
    NVi['Hbar'] = Hbar

    # Calculate number of vesicles of size i per unit volume (NVi)
    NVi_col = []
    for i in NVi.index:
        row = NVi.loc[NVi.index == i]
        i = i
        h_term = (1/row.Hbar)
        common_term = NVi.loc[1, 'Alphas'] * row.NAi
        for iteration in range(i-1):
            j = iteration + 1
            iter_term = NVi.loc[j+1, 'Alphas'] * NVi.loc[i-j, 'NAi']
            common_term = common_term - iter_term
        nv = (h_term * common_term).loc[i]
        NVi_col.append(nv if nv > 0 else 0)
    NVi['NVi'] = NVi_col

    # Calculate Spherical volume
    SpherVol = []
    for i in NVi.index:
        spervol = (np.pi * (NVi['Hbar'][i] ** 3)) / 6
        SpherVol.append(spervol)
    NVi['SpherVol'] = SpherVol

    # Calculate volume fraction
    NVi['VolFrac'] = NVi['NVi'] * NVi['SpherVol']
    TotVolFrac = NVi['VolFrac'].sum()
    Correction_factor = (clast_ves / 100) / TotVolFrac
    NVi['VolFracAdj'] = NVi['VolFrac'] * Correction_factor

    # Calculate cumulative volume fraction
    CumVolFrac = []
    cumvolfrac = 0
    for i in NVi.index:
        current_row = NVi['VolFracAdj'][i]
        cumvolfrac = cumvolfrac + current_row
        CumVolFrac.append(cumvolfrac)
    NVi['CumVolFrac'] = CumVolFrac

    # Calculate total vesicle number densities
    NA = NVi['NAi'].sum()
    NV = NVi['NVi'].sum()
    NVcorr = (NV * 100) / (100 - clast_ves)

    return NVi, NA, NV, NVcorr


def output_data(sample, clast_ves, areas, info, diam, freq, NAi, NVi,
                NA, NV, NVcorr, save=False):
    # Arrange data in a dictionary
    data = {}
    data[sample] = {'Sample number': sample,
                    'Clast porosity': clast_ves / 100,
                    'Areas': areas,
                    'Info': info,
                    'Diameters': diam,
                    'Frequencies': freq,
                    'NAi': NAi,
                    'NVi': NVi,
                    'NA': NA,
                    'NV': NV,
                    'NVcorr': NVcorr}
    if save:
        np.save('sample_data.npy', data)
    return data


def vsa_processor(input_areas, input_info, sample, clast_ves, plot=False,
                  save=False):

    NAi, info, areas, diam, freq = calculate_NAi(input_areas=input_areas,
                                                 input_info=input_info)
    NAi = merge_NAi(NAi=NAi, info=info, sample=sample, plot=plot)
    NVi, NA, NV, NVcorr = calculate_NVi(NAi=NAi, info=info,
                                        clast_ves=clast_ves)
    data = output_data(sample, clast_ves, areas, info, diam, freq, NAi, NVi,
                       NA, NV, NVcorr, save=save)
    return data


areas = pd.read_csv('TestData/sample_areas.csv')
info = pd.read_csv('TestData/sample_info.csv')
data = vsa_processor(input_areas=areas, input_info=info, sample='Test',
                     clast_ves=74.6, save=True)
