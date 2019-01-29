#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################################################
#    <one line to give the program's name and a brief idea of what it does.>
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


def create_bins():
    # Set up the bins
    min_size_px = 11
    scale_factor_250 = 1693.33
    max_size_mm = 5
    bins = [min_size_px/scale_factor_250]
    geobin = bins[0]
    while geobin < max_size_mm:
        geobin *= 10 ** 0.1
        bins.append(geobin)
    return bins


def process_areas(areas):
    # Consolidate area data
    col_scan, col_25, col_100, col_250 = [], [], [], []

    for column in areas.columns:
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

    areas = pd.DataFrame({'scan': areas['scan'],
                          '25': areas[col_25].stack().reset_index()[0],
                          '100': areas[col_100].stack().reset_index()[0],
                          '250': areas[col_250].stack().reset_index()[0]})
    return areas, col_scan, col_25, col_100, col_250


def process_metadata(info):
    # Finish creating information dataframe
    info.index = info['image']
    info['ref_area'] = ((info['width'] * info['height'])
                        / info['scale_factor'] ** 2
                        - ((info['width'] * info['height'])
                           / info['scale_factor'] ** 2
                           * (info['edge_grey'] / 255)))
    info['ves_frac'] = info['vesicles_grey'] / 255
    info['phen_frac'] = info['crystals_grey'] / 255
    return info


def calc_diameters(areas):
    # Calculate diameter from area
    diam = pd.DataFrame({'scan': (4/np.pi * areas['scan']) ** 0.5,
                         '25': (4/np.pi * areas['25']) ** 0.5,
                         '100': (4/np.pi * areas['100']) ** 0.5,
                         '250': (4/np.pi * areas['250']) ** 0.5})
    return diam


def calc_histogram(diam, bins):
    # Create histograms from diameters and bins
    freq = pd.DataFrame({'scan': np.histogram(diam['scan'].dropna(), bins)[0],
                         '25': np.histogram(diam['25'].dropna(), bins)[0],
                         '100': np.histogram(diam['100'].dropna(), bins)[0],
                         '250': np.histogram(diam['250'].dropna(), bins)[0],
                         'bins': bins[:-1]})
    freq = freq[['bins', 'scan', '25', '100', '250']]
    return freq


def calc_NAi(freq, info, col_25, col_100, col_250, bins):
    # Calculate number of bubbles of size i per unit area
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
    return NAi


def merge_NAi(NAi):
    # Merge NAi data to create smooth transition between magnifications
    NAi_merged, sutures, source = [], [], []

    switch_250 = True
    switch_100 = True
    switch_25 = True
    for i in range(len(NAi)):
        if not NAi_merged:
            # Initialise unified Na list with first value from 250 list
            NAi_merged.append(NAi['250'][i])
            source.append('250')
            continue
        if switch_250:
            # Check that using the next 250 value would give a smaller delta
            if NAi_merged[-1] - NAi['250'][i] < NAi_merged[-1] - NAi['100'][i]:
                NAi_merged.append(NAi['250'][i])
                source.append('250')
                continue
            else:
                sutures.append([NAi['bins'][i], NAi['250'][i]])
                switch_250 = False
        if switch_100:
            # Check that using the next 100 value would give a smaller delta
            if NAi_merged[-1] - NAi['100'][i] < NAi_merged[-1] - NAi['25'][i]:
                NAi_merged.append(NAi['100'][i])
                source.append('100')
                continue
            else:
                sutures.append([NAi['bins'][i], NAi['100'][i]])
                switch_100 = False
        if switch_25:
            # Check that using the next 25 value would give a smaller delta
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
    return NAi


def plot_NA(NAi, sample):
    # Plot the various NAi data to check smoothness of merged set
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


def phenocryst_adjustment(info, NAi):
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
    return NAi


def remove_empty_bins(NAi):
    # Remove empty bins
    for value in NAi['NAi']:
        if NAi['NAi'].iloc[-1] == 0:
            NAi = NAi.drop(NAi.index[-1])
            return NAi


def create_NVi_dataframe(NAi):
    # Make new dataframe. Needs new bins to fill the 32 rows needed for the
    # alphas
    NVi = pd.DataFrame()
    geobins = [NAi['bins'].max()]
    for i in range(31):
        newbin = geobins[-1] / 10 ** 0.1
        geobins.append(newbin)
    NVi['bins'] = geobins[::-1]
    return NVi


def populate_NAi(info, NAi, NVi):
    # Insert NAi so that the values for the largest bins are aligned to bottom
    if info['phen_frac'].mean() != 0:
        NAi_temp = list(NAi['NAi adj.'])
    else:
        NAi_temp = list(NAi['NAi'])
    for i in range(len(NVi) - len(NAi)):
        NAi_temp.insert(0, np.nan)
    NVi['NAi'] = NAi_temp
    NVi['NAi'].fillna(0, inplace=True)  # replace NaN with 0s
    return NVi


def populate_vol_class_alpha(NVi):
    # Create volumetric bin and class columns. Set class as index
    NVi['VolBins'] = NVi['bins']**3
    NVi['Class'] = (np.arange(1, 33))[::-1]
    NVi.set_index('Class', drop=True, inplace=True)

    # Set alpha column
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
    return NVi


def calc_Hbar(NVi):
    # Calculating Hbar values
    Hbar = []
    for i in NVi.index:
        hbar = ((NVi['VolBins'][i] + NVi['VolBins'][i]*10**0.3) / 2)**(1/3)
        Hbar.append(hbar)
    NVi['Hbar'] = Hbar
    return NVi


def calc_Nv(NVi):
    # Calculating Nv
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
    return NVi


def calc_spherical_vol(NVi):
    # Calculating Spherical volume
    SpherVol = []
    for i in NVi.index:
        spervol = (np.pi * (NVi['Hbar'][i] ** 3)) / 6
        SpherVol.append(spervol)
    NVi['SpherVol'] = SpherVol
    return NVi


def calc_vol_frac(NVi, clast_ves):
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
    return NVi


def calc_number_densities(NVi, clast_ves):
    # Calculate total vesicle number densities
    NA = NVi['NAi'].sum()
    NV = NVi['NVi'].sum()
    NVcorr = (NV * 100) / (100 - clast_ves)
    return NA, NV, NVcorr


def create_dic(sample, clast_ves, areas, info, diam, freq, NAi, NVi,
               NA, NV, NVcorr):
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
    return data


def process_data(areas, info, sample, clast_ves, plot=False, save=True):

    bins = create_bins()
    areas, col_scan, col_25, col_100, col_250 = process_areas(areas)
    info = process_metadata(info)
    diam = calc_diameters(areas)
    freq = calc_histogram(diam, bins)
    NAi = calc_NAi(freq, info, col_25, col_100, col_250, bins)
    NAi = merge_NAi(NAi)
    if plot:
        plot_NA(NAi=NAi, sample=sample)
    NAi = phenocryst_adjustment(info, NAi)
    NAi = remove_empty_bins(NAi)
    NVi = create_NVi_dataframe(NAi)
    NVi = populate_NAi(info, NAi, NVi)
    NVi = populate_vol_class_alpha(NVi)
    NVi = calc_Hbar(NVi)
    NVi = calc_Nv(NVi)
    NVi = calc_spherical_vol(NVi)
    NVi = calc_vol_frac(NVi, clast_ves)
    NA, NV, NVcorr = calc_number_densities(NVi, clast_ves)
    data = create_dic(sample, clast_ves, areas, info, diam, freq, NAi, NVi,
                      NA, NV, NVcorr)
    if save:
        np.save('sample_data.npy', data)
    return data


areas = pd.read_csv('TestData/sample_areas.csv')
info = pd.read_csv('TestData/sample_info.csv')
data = process_data(areas=areas, info=info, sample='Test', clast_ves=74.6,
                    save=True)
