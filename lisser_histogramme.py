import rasterio
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
import os
import glob
import fnmatch

lignes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
colonnes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']

acquisition_path = 'F:\\optimisation_plate_type_BG_NF\\Cas_ZC_rfl_v1_SCBG250423-11_ZC22_gamme_large_16052025\\expo\\Images'
marge_X = 25
marge_Y = 25
nb_sous_zones_X = 3
nb_sous_zones_Y = 2

window_length = 651

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


def load_tif_image(image_path):
    with rasterio.open(image_path) as src:
        image = src.read(1)  # Lecture de la première bande
    return image


def plot_histogram(image, x_min, x_max, y_min, y_max, window_length):
    roi = image[y_min:y_max, x_min:x_max]

    histogram, bin_edges = np.histogram(roi, bins=4096, range=(0, 4096))

    # Lisser l'histogramme avec un filtre Savitzky-Golay
    polyorder = 3
    smoothed_histogram = savgol_filter(histogram, window_length, polyorder)

    # plt.plot(bin_edges[0:-1], histogram)
    # plt.plot(bin_edges[0:-1], smoothed_histogram)
    # plt.title("Histogramme lissé")
    # plt.xlabel("Valeur des pixels")
    # plt.ylabel("Fréquence")
    # plt.show()

    max_index = np.argmax(smoothed_histogram)
    # max_index = np.argmax(histogram)
    max_value = bin_edges[max_index]

    return max_value


def calcul_liste_intensite(image_path, marge_X, marge_Y, nb_sous_zones_X, nb_sous_zones_Y, window_length):
    image = load_tif_image(image_path)
    taille_image = image.shape
    liste_intensite = []
    dimension_X_sous_zone = int((taille_image[1] - 2 * marge_X * taille_image[1] / 100) / nb_sous_zones_X)
    dimension_Y_sous_zone = int((taille_image[0] - 2 * marge_Y * taille_image[0] / 100) / nb_sous_zones_Y)
    for i in range(nb_sous_zones_Y):
        for j in range(nb_sous_zones_X):
            y_min = int(taille_image[0] * marge_Y / 100) + i * dimension_Y_sous_zone
            y_max = int(taille_image[0] * marge_Y / 100) + (i + 1) * dimension_Y_sous_zone
            x_min = int(taille_image[1] * marge_X / 100) + j * dimension_X_sous_zone
            x_max = int(taille_image[1] * marge_X / 100) + (j + 1) * dimension_X_sous_zone

            max_abscissa = plot_histogram(image, x_min, x_max, y_min, y_max, window_length)
            liste_intensite.append(max_abscissa)
    return liste_intensite


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


intensites_acquisition = []
data_plaque = []

bins = []
bin_edge = MinThickness
while bin_edge <= MaxThickness:
    bins.append(bin_edge)
    bin_edge += ThicknessBinSize
if bins[-1] != MaxThickness:
    bins.append(MaxThickness)

for lettre in range(len(lignes)):
    for chiffre in range(len(colonnes)):
        puits = lignes[lettre] + colonnes[chiffre]
        thickness = []

        # calcul des 2 listes d'intensité 455 et 730
        motif_730 = os.path.join(acquisition_path, puits + '_730(*.tif')
        chemin_image_730 = glob.glob(motif_730)
        motif_455 = os.path.join(acquisition_path, puits + '_455(*.tif')
        chemin_image_455 = glob.glob(motif_455)
        if chemin_image_730 and chemin_image_455:
            liste_455 = calcul_liste_intensite(chemin_image_455[0], marge_X, marge_Y, nb_sous_zones_X, nb_sous_zones_Y,
                                               window_length)
            liste_730 = calcul_liste_intensite(chemin_image_730[0], marge_X, marge_Y, nb_sous_zones_X, nb_sous_zones_Y,
                                               window_length)
            intensites_acquisition.append([liste_455, liste_730])

        # lecture du csv correspondant au puits pour récupérer les valeurs d'intensité de calibration
        chemin_csv_puits = acquisition_path + '\\' + puits + '_455.tif_' + puits + '_730.tif_Reflecto_Results.csv'
        data = np.genfromtxt(chemin_csv_puits, delimiter=';', encoding='utf-8', skip_header=1)
        for i in range(len(data)):
            intensity455 = liste_455[i]
            intensity730 = liste_730[i]
            minRefInt455 = data[i, 9]
            maxRefInt455 = data[i, 10]
            minRefInt730 = data[i, 11]
            maxRefInt730 = data[i, 12]
            epaisseur_calculee = computeThicknessweighted(intensity455, intensity730, minRefInt455, maxRefInt455,
                                                          minRefInt730, maxRefInt730, 0)
            # epaisseur_calculee = computeThickness(intensity455,intensity730,minRefInt455,maxRefInt455,minRefInt730,maxRefInt730)
            thickness.append(epaisseur_calculee)
        thickness = np.asarray(thickness)

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

        data_plaque.append([puits, str(epaisseur_puits)])

        print(puits)

data_plaque = np.asarray(data_plaque)
chemin_enregistrement_csv_plaque = acquisition_path + '\\lisse_histo_savgol_' + str(window_length) + '.csv'
np.savetxt(chemin_enregistrement_csv_plaque, data_plaque, delimiter=';', fmt='%s', comments='')