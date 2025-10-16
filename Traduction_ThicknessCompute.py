# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 17:27:44 2024

@author: Admin
"""

import numpy as np

lignes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
colonnes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
# lignes = ['A','B','C','D']
# lignes = ['B']

# chemin_dossier_acquisition_avec_csv = 'C:\\Users\\Admin\\Downloads\\wetransfer_cas_zc_rfl_v1_cas250327-02_250519-skyr_drinks_galaya_2025-05-27_1501\\'
chemin_dossier_acquisition_avec_csv = 'F:\\valid_ZC28\\plaques_reference_ZC28\\AxhVLA250403-20\\'

# liste_acquisitions = ['Aplat_SCExploCas241010-06_valid_28102024_ZC11_HG_01',
#                       'Aplat_SCExploCas241010-17_valid_28102024_ZC11_HG_01',
#                       'Aplat_SCExploCas240925-06_valid_28102024_ZC11_HG_01',
#                       'Aplat_SCExploCas240925-12_valid_28102024_ZC11_HG_01',
#                       'Aplat_SCExploCas241010-01_valid_29102024_ZC11_HG_01',
#                       'Aplat_SCExploCas240925-01_valid_29102024_ZC11_HG_01',
#                       'Aplat_SCExploCas241010-01_valid_29102024_ZC11_HG_02',
#                       'Aplat_SCExploCas240925-01_valid_29102024_ZC11_HG_02',
#                       'Aplat_SCExploCas240925-12_valid_29102024_ZC11_HG_02',
#                       'Aplat_SCExploCas240925-06_valid_29102024_ZC11_HG_02',
#                       'Aplat_SCExploCas241010-06_valid_29102024_ZC11_HG_02',
#                       'Aplat_SCExploCas241010-17_valid_29102024_ZC11_HG_02',
#                       'Aplat_SCExploCas241010-17_valid_29102024_ZC11_HG_03',
#                       'Aplat_SCExploCas241010-01_valid_29102024_ZC11_HG_03',
#                       'Aplat_SCExploCas240925-12_valid_29102024_ZC11_HG_03',
#                       'Aplat_SCExploCas240925-06_valid_29102024_ZC11_HG_03',
#                       'Aplat_SCExploCas240925-01_valid_29102024_ZC11_HG_03',
#                       'Aplat_SCExploCas241010-06_valid_29102024_ZC11_HG_03']

# liste_acquisitions = ['Aplat_SCExploCas240925-01_valid_ZC19_23102024_HG_01',
#                       'Aplat_SCExploCas240925-01_valid_ZC19_23102024_HG_02',
#                       'Aplat_SCExploCas240925-01_valid_ZC19_23102024_HG_03',
#                       'Aplat_SCExploCas240925-06_valid_ZC19_23102024_HG_01',
#                       'Aplat_SCExploCas240925-06_valid_ZC19_23102024_HG_02',
#                       'Aplat_SCExploCas240925-06_valid_ZC19_23102024_HG_03',
#                       'Aplat_SCExploCas240925-12_valid_ZC19_23102024_HG_01',
#                       'Aplat_SCExploCas240925-12_valid_ZC19_23102024_HG_02',
#                       'Aplat_SCExploCas240925-12_valid_ZC19_23102024_HG_03',
#                       'Aplat_SCExploCas241010-01_valid_ZC19_23102024_HG_01',
#                       'Aplat_SCExploCas241010-01_valid_ZC19_23102024_HG_02',
#                       'Aplat_SCExploCas241010-01_valid_ZC19_23102024_HG_03',
#                       'Aplat_SCExploCas241010-06_valid_ZC19_23102024_HG_01',
#                       'Aplat_SCExploCas241010-06_valid_ZC19_23102024_HG_02',
#                       'Aplat_SCExploCas241010-06_valid_ZC19_23102024_HG_03',
#                       'Aplat_SCExploCas241010-17_valid_ZC19_23102024_HG_01',
#                       'Aplat_SCExploCas241010-17_valid_ZC19_23102024_HG_02',
#                       'Aplat_SCExploCas241010-17_valid_ZC19_23102024_HG_03']

# liste_acquisitions = ['Aplat_SCExploAxHV240425-09_valid_06112024_ZC11_HG_01',
#                       'Aplat_SCExploAxHV240425-10_valid_06112024_ZC11_HG_01',
#                       'Aplat_SCExploAxHV240627-08_valid_06112024_ZC11_HG_01',
#                       'Aplat_SCExploAxHV240722-10_valid_06112024_ZC11_HG_01',
#                       'Aplat_SCExploAxHV240822-06_valid_06112024_ZC11_HG_01',
#                       'Aplat_SCExploAxHV241015-18_valid_06112024_ZC11_HG_01']

# liste_acquisitions = ['Aplat_SCExploAxHV240722-10_valid_07112024_ZC19_HG_03',
#                       'Aplat_SCExploAxHV240722-10_valid_07112024_ZC19_HG_02',
#                       'Aplat_SCExploAxHV240425-09_valid_07112024_ZC19_HG_03',
#                       'Aplat_SCExploAxHV240425-09_valid_07112024_ZC19_HG_02',
#                       'Aplat_SCExploAxHV240822-06_valid_07112024_ZC19_HG_03',
#                       'Aplat_SCExploAxHV240822-06_valid_07112024_ZC19_HG_02',
#                       'Aplat_SCExploAxHV240627-08_valid_07112024_ZC19_HG_03',
#                       'Aplat_SCExploAxHV240627-08_valid_07112024_ZC19_HG_02',
#                       'Aplat_SCExploAxHV241015-18_valid_07112024_ZC19_HG_03',
#                       'Aplat_SCExploAxHV241015-18_valid_07112024_ZC19_HG_02',
#                       'Aplat_SCExploAxHV240425-10_valid_07112024_ZC19_HG_03',
#                       'Aplat_SCExploAxHV240425-10_valid_07112024_ZC19_HG_02',
#                       'Aplat_SCExploAxHV240822-06_valid_06112024_ZC19_HG_01',
#                       'Aplat_SCExploAxHV240627-08_valid_06112024_ZC19_HG_01',
#                       'Aplat_SCExploAxHV241015-18_valid_06112024_ZC19_HG_01',
#                       'Aplat_SCExploAxHV240722-10_valid_06112024_ZC19_HG_01',
#                       'Aplat_SCExploAxHV240425-09_valid_06112024_ZC19_HG_01',
#                       'Aplat_SCExploAxHV240425-10_valid_06112024_ZC19_HG_01']

# liste_acquisitions = ['AxHV_ZC_rfl_v1_SCAxhVLA250122-18_ZC13_CHAUD_28042025','AxHV_ZC_rfl_v1_SCAxhVLA250122-18_ZC13_FROID_28042025','AxHV_ZC_rfl_v1_SCAxhVLA250122-18_ZC13_MOY_01_28042025','AxHV_ZC_rfl_v1_SCAxhVLA250122-18_ZC13_MOY_2_28042025']
# liste_acquisitions = ['Cas_ZC_rfl_v1_Cas250304-01_ZC13_CHAUD_28042025','Cas_ZC_rfl_v1_Cas250304-01_ZC13_FROID_28042025','Cas_ZC_rfl_v1_Cas250304-01_ZC13_MOY_1_28042025','Cas_ZC_rfl_v1_Cas250304-01_ZC13_MOY_2_28042025','Cas_ZC_rfl_v1_Cas250313-17_ZC13_CHAUD_28042025','Cas_ZC_rfl_v1_Cas250313-17_ZC13_FROID_28042025','Cas_ZC_rfl_v1_Cas250313-17_ZC13_MOY_1_28042025','Cas_ZC_rfl_v1_Cas250313-17_ZC13_MOY_2_28042025','AxHV_ZC_rfl_v1_AxhVLA250206-15_ZC13_CHAUD_28042025','AxHV_ZC_rfl_v1_AxhVLA250206-15_ZC13_FROID_28042025','AxHV_ZC_rfl_v1_AxhVLA250206-15_ZC13_MOY_1_28042025','AxHV_ZC_rfl_v1_AxhVLA250206-15_ZC13_MOY_2_28042025']
# liste_acquisitions = ['AxHV_ZC_rfl_v1_SCAxhVLA250122-18_ZC13_MOY_2_28042025']
# liste_acquisitions = ['Cas_ZC_rfl_v1_SCBG250423-11_ZC22_gamme_large_16052025','Cas_ZC_rfl_v1_SCBG250423-20_gamme_large_ZC22_19052025']
# liste_acquisitions = ['Cas_ZC_rfl_v1_Cas250327-02_250519-Skyr_Drinks_GALAYA']
# liste_acquisitions = ['AxHV_ZC_rfl_v1_AxhVLA250408-14_bug_ZU_neg_ZC21_01','AxHV_ZC_rfl_v1_AxhVLA250408-14_bug_ZU_neg_ZC21_02']
liste_acquisitions = ['AxHV_ZC_rfl_v1_AxhVLA250403-20_REFERENCE_ZYMOPTIQ_EXPERT']

""" Paramètres Réflecto """
minIntensityRange455 = 150
maxIntensityRange455 = 250
minIntensityRange730 = 150
maxIntensityRange730 = 250
maxIntensityRatioRejectedT = 0.97
angle = 0
maxThickness = 200
polymerIndice730 = 1.50
# polymerIndice730 = 1.51
# polymerIndice730 = 1.49
siliciumIndice730 = 3.65
polymerIndice455 = 1.51
# polymerIndice455 = 1.545
# polymerIndice455 = 1.48
siliciumIndice455 = 4.55

""" Paramètres Statistiques """
MinAreaCount = 2
MinThickness = 20
MaxThickness = 200
ThicknessBinSize = 3
MinThicknessFrequency = 2
ThicknessDefinition = 'Median'


def computeThickness(intensity455, intensity730, minRefInt455, maxRefInt455, minRefInt730, maxRefInt730):
    res = np.nan
    if intensity455 < minRefInt455:
        if (minRefInt455 - intensity455) <= minIntensityRange455:
            intensity455 = minRefInt455
        else:
            return np.nan
    elif intensity455 > maxRefInt455:
        if (intensity455 - maxRefInt455) <= maxIntensityRange455:
            intensity455 = maxRefInt455
        else:
            return np.nan
    if intensity730 < minRefInt730:
        if (minRefInt730 - intensity730) <= minIntensityRange730:
            intensity730 = minRefInt730
        else:
            return np.nan
    elif intensity730 > maxRefInt730:
        if (intensity730 - maxRefInt730) <= maxIntensityRange730:
            intensity730 = maxRefInt730
        else:
            return np.nan
    if (intensity455 > maxIntensityRatioRejectedT * maxRefInt455) and (
            intensity730 > maxIntensityRatioRejectedT * maxRefInt730):
        return np.nan
    if (intensity455 >= minRefInt455) and (intensity455 <= maxRefInt455) and (intensity730 >= minRefInt730) and (
            intensity730 <= maxRefInt730):
        normalizeIntensity455 = (intensity455 - minRefInt455) / (maxRefInt455 - minRefInt455)
        normalizeIntensity730 = (intensity730 - minRefInt730) / (maxRefInt730 - minRefInt730)

        thicknessList455 = []

        n0 = 1.0
        n1 = polymerIndice455
        n2 = siliciumIndice455
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity455 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 455
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList455.append(float(h1))

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList455.append(float(h2))

                phase1 += np.pi
                phase2 += np.pi

        thicknessList730 = []

        n0 = 1.0
        n1 = polymerIndice730
        n2 = siliciumIndice730
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity730 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 730
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList730.append(float(h1))

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList730.append(float(h2))

                phase1 += np.pi
                phase2 += np.pi

            index455 = -1
            index730 = -1
            delta = -1

            for i in range(len(thicknessList455)):
                h1 = thicknessList455[i]
                for j in range(len(thicknessList730)):
                    h2 = thicknessList730[j]
                    if ((delta == -1) or (delta > abs(h2 - h1))):
                        delta = abs(h2 - h1)
                        index455 = i
                        index730 = j
            if (delta >= 0):
                res = (thicknessList455[index455] + thicknessList730[index730]) / 2
    return res


def computeThicknessweighted(intensity455, intensity730, minRefInt455, maxRefInt455, minRefInt730, maxRefInt730,
                             ResidualWeight):
    res = np.nan
    if intensity455 < minRefInt455:
        if (minRefInt455 - intensity455) <= minIntensityRange455:
            intensity455 = minRefInt455
        else:
            return np.nan
    elif intensity455 > maxRefInt455:
        if (intensity455 - maxRefInt455) <= maxIntensityRange455:
            intensity455 = maxRefInt455
        else:
            return np.nan
    if intensity730 < minRefInt730:
        if (minRefInt730 - intensity730) <= minIntensityRange730:
            intensity730 = minRefInt730
        else:
            return np.nan
    elif intensity730 > maxRefInt730:
        if (intensity730 - maxRefInt730) <= maxIntensityRange730:
            intensity730 = maxRefInt730
        else:
            return np.nan
    if (intensity455 > maxIntensityRatioRejectedT * maxRefInt455) and (
            intensity730 > maxIntensityRatioRejectedT * maxRefInt730):
        return np.nan
    if (intensity455 >= minRefInt455) and (intensity455 <= maxRefInt455) and (intensity730 >= minRefInt730) and (
            intensity730 <= maxRefInt730):
        normalizeIntensity455 = (intensity455 - minRefInt455) / (maxRefInt455 - minRefInt455)
        normalizeIntensity730 = (intensity730 - minRefInt730) / (maxRefInt730 - minRefInt730)

        thicknessList455 = []

        n0 = 1.0
        n1 = polymerIndice455
        n2 = siliciumIndice455
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity455 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 455
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList455.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList455.append(h2)

                phase1 += np.pi
                phase2 += np.pi

        thicknessList730 = []

        n0 = 1.0
        n1 = polymerIndice730
        n2 = siliciumIndice730
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity730 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 730
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList730.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList730.append(h2)

                phase1 += np.pi
                phase2 += np.pi

            index455 = -1
            index730 = -1
            delta = -1

            for i in range(len(thicknessList455)):
                h1 = thicknessList455[i]
                for j in range(len(thicknessList730)):
                    h2 = thicknessList730[j]
                    if ((delta == -1) or (delta > abs(h2 - h1))):
                        delta = abs(h2 - h1)
                        index455 = i
                        index730 = j
            if (delta >= 0):
                weight455 = ResidualWeight + np.sin(
                    (np.pi / (maxRefInt455 - minRefInt455)) * (intensity455 - minRefInt455))
                weight730 = ResidualWeight + np.sin(
                    (np.pi / (maxRefInt730 - minRefInt730)) * (intensity730 - minRefInt730))
                totalWeight = weight455 + weight730
                res = (weight455 * thicknessList455[index455] + weight730 * thicknessList730[index730]) / totalWeight
    return res


def computeThicknessweightedrestricted(intensity455, intensity730, minRefInt455, maxRefInt455, minRefInt730,
                                       maxRefInt730, ResidualWeight, restriction_rate):
    res = np.nan
    if intensity455 < minRefInt455:
        if (minRefInt455 - intensity455) <= minIntensityRange455:
            intensity455 = minRefInt455
        else:
            return np.nan
    elif intensity455 > maxRefInt455:
        if (intensity455 - maxRefInt455) <= maxIntensityRange455:
            intensity455 = maxRefInt455
        else:
            return np.nan
    if intensity730 < minRefInt730:
        if (minRefInt730 - intensity730) <= minIntensityRange730:
            intensity730 = minRefInt730
        else:
            return np.nan
    elif intensity730 > maxRefInt730:
        if (intensity730 - maxRefInt730) <= maxIntensityRange730:
            intensity730 = maxRefInt730
        else:
            return np.nan
    if (intensity455 > maxIntensityRatioRejectedT * maxRefInt455) and (
            intensity730 > maxIntensityRatioRejectedT * maxRefInt730):
        return np.nan
    ''' NOUVEAUX TEST '''
    if intensity455 < (1 + restriction_rate) * minRefInt455 and intensity455 > minRefInt455:
        intensity455 = minRefInt455
    if intensity455 > (1 - restriction_rate) * maxRefInt455 and intensity455 < maxRefInt455:
        intensity455 = maxRefInt455
    if intensity730 < (1 + restriction_rate) * minRefInt730 and intensity730 > minRefInt730:
        intensity730 = minRefInt730
    if intensity730 > (1 - restriction_rate) * maxRefInt730 and intensity730 < maxRefInt730:
        intensity730 = maxRefInt730
    ''' '''

    if (intensity455 >= minRefInt455) and (intensity455 <= maxRefInt455) and (intensity730 >= minRefInt730) and (
            intensity730 <= maxRefInt730):
        normalizeIntensity455 = (intensity455 - minRefInt455) / (maxRefInt455 - minRefInt455)
        normalizeIntensity730 = (intensity730 - minRefInt730) / (maxRefInt730 - minRefInt730)

        thicknessList455 = []

        n0 = 1.0
        n1 = polymerIndice455
        n2 = siliciumIndice455
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity455 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 455
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList455.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList455.append(h2)

                phase1 += np.pi
                phase2 += np.pi

        thicknessList730 = []

        n0 = 1.0
        n1 = polymerIndice730
        n2 = siliciumIndice730
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity730 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 730
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList730.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList730.append(h2)

                phase1 += np.pi
                phase2 += np.pi

            # print('liste épaisseurs 455 : ',thicknessList455, '\nliste épaisseurs 730 : ',thicknessList730)

            index455 = -1
            index730 = -1
            delta = -1
            ecart = 0

            for i in range(len(thicknessList455)):
                h1 = thicknessList455[i]
                for j in range(len(thicknessList730)):
                    h2 = thicknessList730[j]
                    if ((delta == -1) or (delta > abs(h2 - h1))):
                        delta = abs(h2 - h1)
                        index455 = i
                        index730 = j
            if (delta >= 0):
                weight455 = ResidualWeight + np.sin(
                    (np.pi / (maxRefInt455 - minRefInt455)) * (intensity455 - minRefInt455))
                weight730 = ResidualWeight + np.sin(
                    (np.pi / (maxRefInt730 - minRefInt730)) * (intensity730 - minRefInt730))
                totalWeight = weight455 + weight730
                res = (weight455 * thicknessList455[index455] + weight730 * thicknessList730[index730]) / totalWeight
                ecart = abs(thicknessList455[index455] - thicknessList730[index730])

    return res, ecart


def computeThicknessweightedrestrictedabsolute(intensity455, intensity730, minRefInt455, maxRefInt455, minRefInt730,
                                               maxRefInt730, ResidualWeight, marge_restriction):
    res = np.nan
    if intensity455 < minRefInt455:
        if (minRefInt455 - intensity455) <= minIntensityRange455:
            intensity455 = minRefInt455
        else:
            return np.nan
    elif intensity455 > maxRefInt455:
        if (intensity455 - maxRefInt455) <= maxIntensityRange455:
            intensity455 = maxRefInt455
        else:
            return np.nan
    if intensity730 < minRefInt730:
        if (minRefInt730 - intensity730) <= minIntensityRange730:
            intensity730 = minRefInt730
        else:
            return np.nan
    elif intensity730 > maxRefInt730:
        if (intensity730 - maxRefInt730) <= maxIntensityRange730:
            intensity730 = maxRefInt730
        else:
            return np.nan
    if (intensity455 > maxIntensityRatioRejectedT * maxRefInt455) and (
            intensity730 > maxIntensityRatioRejectedT * maxRefInt730):
        return np.nan
    ''' NOUVEAUX TEST '''
    if intensity455 < minRefInt455 + marge_restriction and intensity455 > minRefInt455:
        intensity455 = minRefInt455
    if intensity455 > maxRefInt455 - marge_restriction and intensity455 < maxRefInt455:
        intensity455 = maxRefInt455
    if intensity730 < minRefInt730 + marge_restriction and intensity730 > minRefInt730:
        intensity730 = minRefInt730
    if intensity730 > maxRefInt730 - marge_restriction and intensity730 < maxRefInt730:
        intensity730 = maxRefInt730
    ''' '''

    if (intensity455 >= minRefInt455) and (intensity455 <= maxRefInt455) and (intensity730 >= minRefInt730) and (
            intensity730 <= maxRefInt730):
        normalizeIntensity455 = (intensity455 - minRefInt455) / (maxRefInt455 - minRefInt455)
        normalizeIntensity730 = (intensity730 - minRefInt730) / (maxRefInt730 - minRefInt730)

        thicknessList455 = []

        n0 = 1.0
        n1 = polymerIndice455
        n2 = siliciumIndice455
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity455 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 455
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList455.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList455.append(h2)

                phase1 += np.pi
                phase2 += np.pi

        thicknessList730 = []

        n0 = 1.0
        n1 = polymerIndice730
        n2 = siliciumIndice730
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity730 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 730
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList730.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList730.append(h2)

                phase1 += np.pi
                phase2 += np.pi

            index455 = -1
            index730 = -1
            delta = -1

            for i in range(len(thicknessList455)):
                h1 = thicknessList455[i]
                for j in range(len(thicknessList730)):
                    h2 = thicknessList730[j]
                    if ((delta == -1) or (delta > abs(h2 - h1))):
                        delta = abs(h2 - h1)
                        index455 = i
                        index730 = j
            if (delta >= 0):
                weight455 = ResidualWeight + np.sin(
                    (np.pi / (maxRefInt455 - minRefInt455)) * (intensity455 - minRefInt455))
                weight730 = ResidualWeight + np.sin(
                    (np.pi / (maxRefInt730 - minRefInt730)) * (intensity730 - minRefInt730))
                totalWeight = weight455 + weight730
                res = (weight455 * thicknessList455[index455] + weight730 * thicknessList730[index730]) / totalWeight
    return res


def computeThicknessweightedrestricted_02(intensity455, intensity730, minRefInt455, maxRefInt455, minRefInt730,
                                          maxRefInt730, ResidualWeight, restriction_rate):
    res = np.nan
    if intensity455 < minRefInt455:
        if (minRefInt455 - intensity455) <= minIntensityRange455:
            intensity455 = minRefInt455
        else:
            return np.nan
    elif intensity455 > maxRefInt455:
        if (intensity455 - maxRefInt455) <= maxIntensityRange455:
            intensity455 = maxRefInt455
        else:
            return np.nan
    if intensity730 < minRefInt730:
        if (minRefInt730 - intensity730) <= minIntensityRange730:
            intensity730 = minRefInt730
        else:
            return np.nan
    elif intensity730 > maxRefInt730:
        if (intensity730 - maxRefInt730) <= maxIntensityRange730:
            intensity730 = maxRefInt730
        else:
            return np.nan
    if (intensity455 > maxIntensityRatioRejectedT * maxRefInt455) and (
            intensity730 > maxIntensityRatioRejectedT * maxRefInt730):
        return np.nan
    ''' NOUVEAUX TEST '''
    if intensity455 < (1 + restriction_rate) * minRefInt455 and intensity455 > minRefInt455:
        intensity455 = minRefInt455
    if intensity455 > (1 - restriction_rate) * maxRefInt455 and intensity455 < maxRefInt455:
        intensity455 = maxRefInt455
    if intensity730 < (1 + restriction_rate) * minRefInt730 and intensity730 > minRefInt730:
        intensity730 = minRefInt730
    if intensity730 > (1 - restriction_rate) * maxRefInt730 and intensity730 < maxRefInt730:
        intensity730 = maxRefInt730
    ''' '''

    if (intensity455 >= minRefInt455) and (intensity455 <= maxRefInt455) and (intensity730 >= minRefInt730) and (
            intensity730 <= maxRefInt730):
        normalizeIntensity455 = (intensity455 - minRefInt455) / (maxRefInt455 - minRefInt455)
        normalizeIntensity730 = (intensity730 - minRefInt730) / (maxRefInt730 - minRefInt730)

        thicknessList455 = []

        n0 = 1.0
        n1 = polymerIndice455
        n2 = siliciumIndice455
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity455 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 455
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList455.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList455.append(h2)

                phase1 += np.pi
                phase2 += np.pi

        thicknessList730 = []

        n0 = 1.0
        n1 = polymerIndice730
        n2 = siliciumIndice730
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity730 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 730
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList730.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList730.append(h2)

                phase1 += np.pi
                phase2 += np.pi

            index455 = -1
            index730 = -1
            delta = -1

            for i in range(len(thicknessList455)):
                h1 = thicknessList455[i]
                for j in range(len(thicknessList730)):
                    h2 = thicknessList730[j]
                    if ((delta == -1) or (delta > abs(h2 - h1))):
                        delta = abs(h2 - h1)
                        index455 = i
                        index730 = j
            if (delta >= 0):
                weight455 = ResidualWeight + np.sin(
                    (np.pi / (maxRefInt455 * (1 - restriction_rate) - minRefInt455 * (1 + restriction_rate))) * (
                                intensity455 - minRefInt455 * (1 + restriction_rate)))
                weight730 = ResidualWeight + np.sin(
                    (np.pi / (maxRefInt730 * (1 - restriction_rate) - minRefInt730 * (1 + restriction_rate))) * (
                                intensity730 - minRefInt730 * (1 + restriction_rate)))
                totalWeight = weight455 + weight730
                res = (weight455 * thicknessList455[index455] + weight730 * thicknessList730[index730]) / totalWeight
    return res


def computeThicknessweightedrestrictedabsolute_02(intensity455, intensity730, minRefInt455, maxRefInt455, minRefInt730,
                                                  maxRefInt730, ResidualWeight, marge_restriction):
    res = np.nan
    if intensity455 < minRefInt455:
        if (minRefInt455 - intensity455) <= minIntensityRange455:
            intensity455 = minRefInt455
        else:
            return np.nan
    elif intensity455 > maxRefInt455:
        if (intensity455 - maxRefInt455) <= maxIntensityRange455:
            intensity455 = maxRefInt455
        else:
            return np.nan
    if intensity730 < minRefInt730:
        if (minRefInt730 - intensity730) <= minIntensityRange730:
            intensity730 = minRefInt730
        else:
            return np.nan
    elif intensity730 > maxRefInt730:
        if (intensity730 - maxRefInt730) <= maxIntensityRange730:
            intensity730 = maxRefInt730
        else:
            return np.nan
    if (intensity455 > maxIntensityRatioRejectedT * maxRefInt455) and (
            intensity730 > maxIntensityRatioRejectedT * maxRefInt730):
        return np.nan
    ''' NOUVEAUX TEST '''
    if intensity455 < minRefInt455 + marge_restriction and intensity455 > minRefInt455:
        intensity455 = minRefInt455
    if intensity455 > maxRefInt455 - marge_restriction and intensity455 < maxRefInt455:
        intensity455 = maxRefInt455
    if intensity730 < minRefInt730 + marge_restriction and intensity730 > minRefInt730:
        intensity730 = minRefInt730
    if intensity730 > maxRefInt730 - marge_restriction and intensity730 < maxRefInt730:
        intensity730 = maxRefInt730
    ''' '''

    if (intensity455 >= minRefInt455) and (intensity455 <= maxRefInt455) and (intensity730 >= minRefInt730) and (
            intensity730 <= maxRefInt730):
        normalizeIntensity455 = (intensity455 - minRefInt455) / (maxRefInt455 - minRefInt455)
        normalizeIntensity730 = (intensity730 - minRefInt730) / (maxRefInt730 - minRefInt730)

        thicknessList455 = []

        n0 = 1.0
        n1 = polymerIndice455
        n2 = siliciumIndice455
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity455 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 455
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList455.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList455.append(h2)

                phase1 += np.pi
                phase2 += np.pi

        thicknessList730 = []

        n0 = 1.0
        n1 = polymerIndice730
        n2 = siliciumIndice730
        r01 = (n0 - n1) / (n0 + n1)
        r12 = (n1 - n2) / (n1 + n2)
        Rmin = (r01 - r12) * (r01 - r12) / (1 + r01 * r01 * r12 * r12 - 2 * r01 * r12)
        Rmax = (r01 + r12) * (r01 + r12) / (1 + r01 * r01 * r12 * r12 + 2 * r01 * r12)
        R = 0.0

        if Rmax < Rmin:
            val = Rmax
            Rmin = Rmax
            Rmax = val

        R = Rmin + normalizeIntensity730 * (Rmax - Rmin)

        if R < 0:
            R = 0

        if R < 1:
            llambda = 730
            a = r01 * r01 + r12 * r12
            b = 1 + r01 * r01 * r12 * r12
            c = 2 * r01 * r12
            d = (b * R - a) / (c * (1.0 - R))

            phase1 = np.nan
            phase2 = np.nan

            if d > 1:
                d = 1
            elif d < -1:
                d = -1

            phase1 = 0.5 * np.arccos(d)
            phase2 = np.pi - phase1

            while True:
                h1 = llambda * phase1 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h1) or (h1 > maxThickness)):
                    break
                thicknessList730.append(h1)

                h2 = llambda * phase2 / (2 * np.pi * n1 * np.cos(angle))
                if (np.isnan(h2) or (h2 > maxThickness)):
                    break

                if (h1 != h2):
                    thicknessList730.append(h2)

                phase1 += np.pi
                phase2 += np.pi

            index455 = -1
            index730 = -1
            delta = -1

            for i in range(len(thicknessList455)):
                h1 = thicknessList455[i]
                for j in range(len(thicknessList730)):
                    h2 = thicknessList730[j]
                    if ((delta == -1) or (delta > abs(h2 - h1))):
                        delta = abs(h2 - h1)
                        index455 = i
                        index730 = j
            if (delta >= 0):
                weight455 = ResidualWeight + np.sin(
                    (np.pi / ((maxRefInt455 - marge_restriction) - (minRefInt455 + marge_restriction))) * (
                                intensity455 - (minRefInt455 + marge_restriction)))
                weight730 = ResidualWeight + np.sin(
                    (np.pi / ((maxRefInt730 - marge_restriction) - (minRefInt730 + marge_restriction))) * (
                                intensity730 - (minRefInt730 + marge_restriction)))
                totalWeight = weight455 + weight730
                res = (weight455 * thicknessList455[index455] + weight730 * thicknessList730[index730]) / totalWeight
    return res


# res = computeThicknessweighted(38.3,80,15,75,35,95,0.5)

for acquisition in liste_acquisitions:

    data_plaque = []

    bins = []
    bin_edge = MinThickness
    while bin_edge <= MaxThickness:
        bins.append(bin_edge)
        bin_edge += ThicknessBinSize
    if bins[-1] != MaxThickness:
        bins.append(MaxThickness)

    for lettre in lignes:
        for chiffre in colonnes:
            print(lettre + chiffre)
            # lettre = 'A'
            # chiffre = '1'
            chemin_csv_puits = chemin_dossier_acquisition_avec_csv + acquisition + '\\Images\\' + lettre + chiffre + '_455.tif_' + lettre + chiffre + '_730.tif_Reflecto_Results.csv'
            data = np.genfromtxt(chemin_csv_puits, delimiter=';', encoding='utf-8', skip_header=1)
            thickness = []
            ecarts_puits = []
            # print('\n\n' + lettre + chiffre)
            for i in range(len(data)):
                intensity455 = data[i, 15]
                intensity730 = data[i, 19]
                minRefInt455 = data[i, 9]
                maxRefInt455 = data[i, 10]
                minRefInt730 = data[i, 11]
                maxRefInt730 = data[i, 12]
                # epaisseur_calculee = computeThicknessweightedrestrictedabsolute(intensity455,intensity730,minRefInt455,maxRefInt455,minRefInt730,maxRefInt730,0,300)
                # epaisseur_calculee = computeThicknessweightedrestrictedabsolute_02(intensity455,intensity730,minRefInt455,maxRefInt455,minRefInt730,maxRefInt730,0,50)
                # epaisseur_calculee = computeThicknessweightedrestricted_02(intensity455,intensity730,minRefInt455,maxRefInt455,minRefInt730,maxRefInt730,0,0.09)
                epaisseur_calculee, ecart_calcule = computeThicknessweightedrestricted(intensity455, intensity730,
                                                                                       minRefInt455, maxRefInt455,
                                                                                       minRefInt730, maxRefInt730, 0,
                                                                                       0.20)
                # epaisseur_calculee = computeThicknessweighted(intensity455,intensity730,minRefInt455,maxRefInt455,minRefInt730,maxRefInt730,0)
                # epaisseur_calculee = computeThickness(intensity455,intensity730,minRefInt455,maxRefInt455,minRefInt730,maxRefInt730)
                thickness.append(epaisseur_calculee)
                ecarts_puits.append(ecart_calcule)
            thickness = np.asarray(thickness)
            ecarts_puits = np.asarray(ecarts_puits)

            # chemin_enregistrement_csv = chemin_dossier_acquisition_avec_csv + acquisition + '\\' + lettre + chiffre + '.csv'
            # np.savetxt(chemin_enregistrement_csv, thickness, delimiter=';', fmt='%f', comments='')

            epaisseur_puits = np.nan
            histogramme = []
            if len(thickness) >= MinAreaCount:
                for i in range(len(bins) - 1):
                    temp = []
                    for j in range(len(thickness)):
                        if thickness[j] >= bins[i] and thickness[j] < bins[i + 1]:
                            temp.append(thickness[j])
                    if len(temp) >= MinThicknessFrequency:
                        histogramme.append(temp)
                valid_thickness = []
                for i in range(len(histogramme)):
                    for j in range(len(histogramme[i])):
                        valid_thickness.append(histogramme[i][j])
                valid_thickness = np.asarray(valid_thickness)
                if ThicknessDefinition == 'Median':
                    epaisseur_puits = np.median(valid_thickness)

            data_plaque.append(
                [lettre + chiffre, str(epaisseur_puits), str(np.mean(ecarts_puits)), str(np.max(ecarts_puits))])

    data_plaque = np.asarray(data_plaque)
    chemin_enregistrement_csv_plaque = chemin_dossier_acquisition_avec_csv + acquisition + '\\Images\\compute_thickness_restr_20.csv'
    np.savetxt(chemin_enregistrement_csv_plaque, data_plaque, delimiter=';', fmt='%s', comments='')


def reconstruction_moy_pond(chemin_acquisition, n4, n7):
    lignes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    colonnes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']

    """ Paramètres Réflecto """
    minIntensityRange455 = 150
    maxIntensityRange455 = 250
    minIntensityRange730 = 150
    maxIntensityRange730 = 250
    maxIntensityRatioRejectedT = 0.97
    angle = 0
    maxThickness = 200
    # polymerIndice730 = 1.49
    polymerIndice730 = n7
    siliciumIndice730 = 3.65
    # polymerIndice455 = 1.48
    polymerIndice455 = n4
    siliciumIndice455 = 4.55

    """ Paramètres Statistiques """
    MinAreaCount = 2
    MinThickness = 20
    MaxThickness = 210
    ThicknessBinSize = 3
    MinThicknessFrequency = 2
    ThicknessDefinition = 'Median'

    data_plaque = []

    bins = []
    bin_edge = MinThickness
    while bin_edge <= MaxThickness:
        bins.append(bin_edge)
        bin_edge += ThicknessBinSize
    if bins[-1] != MaxThickness:
        bins.append(MaxThickness)

    for lettre in lignes:
        for chiffre in colonnes:
            # lettre = 'A'
            # chiffre = '1'
            chemin_csv_puits = chemin_acquisition + '\\' + lettre + chiffre + '_455.tif_' + lettre + chiffre + '_730.tif_Reflecto_Results.csv'
            data = np.genfromtxt(chemin_csv_puits, delimiter=';', encoding='utf-8', skip_header=1)
            thickness = []
            for i in range(len(data)):
                intensity455 = data[i, 15]
                intensity730 = data[i, 19]
                minRefInt455 = data[i, 9]
                maxRefInt455 = data[i, 10]
                minRefInt730 = data[i, 11]
                maxRefInt730 = data[i, 12]
                epaisseur_calculee = computeThicknessweighted(intensity455, intensity730, minRefInt455, maxRefInt455,
                                                              minRefInt730, maxRefInt730, 0)
                # epaisseur_calculee = computeThickness(intensity455,intensity730,minRefInt455,maxRefInt455,minRefInt730,maxRefInt730)
                thickness.append(epaisseur_calculee)
            thickness = np.asarray(thickness)
            chemin_enregistrement_csv = chemin_acquisition + '\\' + lettre + chiffre + '.csv'
            np.savetxt(chemin_enregistrement_csv, thickness, delimiter=';', fmt='%f', comments='')

            epaisseur_puits = np.nan
            histogramme = []
            if len(thickness) >= MinAreaCount:
                for i in range(len(bins) - 1):
                    temp = []
                    for j in range(len(thickness)):
                        if thickness[j] >= bins[i] and thickness[j] < bins[i + 1]:
                            temp.append(thickness[j])
                    if len(temp) >= MinThicknessFrequency:
                        histogramme.append(temp)
                valid_thickness = []
                for i in range(len(histogramme)):
                    for j in range(len(histogramme[i])):
                        valid_thickness.append(histogramme[i][j])
                valid_thickness = np.asarray(valid_thickness)
                if ThicknessDefinition == 'Median':
                    epaisseur_puits = np.median(valid_thickness)

            data_plaque.append([lettre + chiffre, str(epaisseur_puits)])

    data_plaque = np.asarray(data_plaque)
    chemin_enregistrement_csv_plaque = chemin_acquisition + '\\Synthese_Reflecto_MoyPond_data.csv'
    np.savetxt(chemin_enregistrement_csv_plaque, data_plaque, delimiter=';', fmt='%s', comments='')
