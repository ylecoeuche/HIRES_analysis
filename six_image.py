##Author: Yannick Lecoeuche
##Purpose: Returns a spectral plot given the start wavelength and wavelength range

import numpy as np
from astropy.io import fits as pf
import matplotlib
from matplotlib import pyplot as plt
import lmfit
from lmfit import minimize, Parameters, report_fit, fit_report
from IPython.display import Image

wave_values_b = pf.open('keck_bwav.fits')[0].data
wave_values_r = pf.open('keck_rwav.fits')[0].data
wave_values_i = pf.open('keck_iwav.fits')[0].data

six_start = [3958, 4298, 5160, 5259, 5889, 6553]
six_range = [20,   20,   50,   20,   20,   20]

def six_image(spec_file):

    sections = ['b', 'r', 'i']

    b_file = 'b' + str(spec_file[1:])
    r_file = 'r' + str(spec_file[1:])
    i_file = 'i' + str(spec_file[1:])

    image_b = pf.open(b_file)[0].data
    image_r = pf.open(r_file)[0].data
    image_i = pf.open(i_file)[0].data

    for i in range (3):
        wavelength = eval('wave_values_' + sections[i])
        image = eval('image_' + sections[i])
        order_max = wavelength.shape[0]
        wavelength_max = wavelength.shape[1]

        for s in range(order_max):
            order_length = max(wavelength[s]) - min(wavelength[s])
            wave_center = wavelength[s] - min(wavelength[s]) - order_length/2.
            norm_spectrum = image[s]/max(image[s])

            coef = np.polyfit(wave_center, norm_spectrum, 5)

            curve_fit = np.poly1d(coef)

            image[s] = norm_spectrum/curve_fit(wave_center)

    print ('Any empty sections in the graph are not covered in the spectrum.')

    plt.plot([1,2,3])

    for g in range(6):

        start_value = six_start[g]
        wavelength_range = six_range[g]

        sections = ['b', 'r', 'i']
        end_value = start_value + wavelength_range

        for i in range (3):
            wavelength = eval('wave_values_' + sections[i])
            image = eval('image_' + sections[i])
            order_max = wavelength.shape[0]
            wavelength_max = wavelength.shape[1]

            if i == 0:
                if start_value < wavelength[0,0]:
                    start_file = 0
                    start_order = 0
                    start_wavelength = 0

            if i == 2:
                if end_value > wavelength[order_max - 1,wavelength_max - 1]:
                    end_file = 2
                    end_order = order_max - 1
                    end_wavelength = wavelength_max - 1


            if wavelength[0, 0] <= start_value < wavelength[order_max - 1, wavelength_max - 1]:
                start_file = i
                for n in range(order_max):
                    if wavelength[n, 0] <= start_value < wavelength[n, wavelength_max - 1]:
                        start_order = n
                        for m in range(wavelength_max):
                            if abs(start_value - wavelength[n, m]) < .02:
                                start_wavelength = m
                for n in range(order_max - 1):
                    if wavelength[n, wavelength_max - 1] <= start_value < wavelength[n+1,0]:
                        start_order = n + 1
                        start_wavelength = 0

            if wavelength[0, 0] <= end_value < wavelength[order_max - 1, wavelength_max - 1]:
                end_file = i
                for n in range(order_max):
                    if wavelength[n, 0] <= end_value < wavelength[n, wavelength_max - 1]:
                        end_order = n
                        for m in range(wavelength_max):
                            if abs(end_value - wavelength[n, m]) < .02:
                                end_wavelength = m
                        break
                for n in range(order_max - 1):
                    if wavelength[n, wavelength_max - 1] <= end_value < wavelength[n+1,0]:
                        end_order = n
                        end_wavelength = wavelength_max - 1

        for i in range (2):
            wavelength1 = eval('wave_values_' + sections[i])
            image1 = eval('image_' + sections[i])
            wavelength2 = eval('wave_values_' + sections[i+1])
            image2 = eval('image_' + sections[i+1])
            order_max1 = wavelength1.shape[0]
            wavelength_max1 = wavelength1.shape[1]
            order_max2 = wavelength2.shape[0]
            wavelength_max2 = wavelength2.shape[1]

            if wavelength1[order_max1-1, wavelength_max1-1] <= start_value < wavelength2[0, 0]:
                start_file = i + 1
                start_order = 0
                start_wavelength = 0

            if wavelength1[order_max1-1, wavelength_max1-1] <= end_value < wavelength2[0, 0]:
                end_file = i
                end_order = order_max1 - 1
                end_wavelength = wavelength_max1 - 1


        file_range = end_file - start_file
        order_range = end_order - start_order

        plt.subplot(3, 2, g + 1)

        if file_range == 0 and order_range == 0:
            wavelength = eval('wave_values_' + sections[start_file])
            image = eval('image_' + sections[start_file])
            x = wavelength[start_order, start_wavelength:end_wavelength]
            y = image[start_order, start_wavelength:end_wavelength]
            plt.plot(x,y)
            axes = plt.gca()
            axes.set_xlim(start_value, end_value)
            axes.set_ylim(0, 2)
            plt.xlabel('Wavelength (A)')
            plt.ylabel('Brightness')

        if file_range == 0 and order_range != 0:
            wavelength = eval('wave_values_' + sections[start_file])
            image = eval('image_' + sections[start_file])
            x_start = wavelength[start_order, start_wavelength:]
            y_start = image[start_order, start_wavelength:]
            x_end = wavelength[end_order, :end_wavelength]
            y_end = image[end_order, :end_wavelength]
            plt.plot(x_start, y_start)
            plt.plot(x_end, y_end)

            for i in range(start_order + 1, end_order):
                x = wavelength[i]
                y = image[i]
                plt.plot(x,y)

            axes = plt.gca()
            axes.set_xlim(start_value, end_value)
            axes.set_ylim(0, 2)
            plt.xlabel('Wavelength (A)')
            plt.ylabel('Brightness')

        if file_range != 0:
            wavelength1 = eval('wave_values_' + sections[start_file])
            image1 = eval('image_' + sections[start_file])
            order_max1 = wavelength1.shape[0]
            wavelength2 = eval('wave_values_' + sections[end_file])
            image2 = eval('image_' + sections[end_file])
            x_start = wavelength1[start_order, start_wavelength:]
            y_start = image1[start_order, start_wavelength:]
            x_end = wavelength2[end_order, :end_wavelength]
            y_end = image2[end_order, :end_wavelength]
            plt.plot(x_start, y_start)
            plt.plot(x_end, y_end)

            for i in range(start_order + 1, order_max1):
                x = wavelength1[i]
                y = image1[i]
                plt.plot(x,y)

            for l in range(0, end_order):
                x = wavelength2[l]
                y = image2[l]
                plt.plot(x,y)

            for s in range(start_file + 1, end_file):
                wavelength = eval('wave_values_' + sections[s])
                image = eval('image_' + sections[s])
                order_max = wavelength.shape[0]
                for t in range(0, order_max):
                    x = wavelength[t]
                    y = image[t]
                    plt.plot(x,y)

            axes = plt.gca()
            axes.set_xlim(start_value, end_value)
            axes.set_ylim(0, 2)
            plt.xlabel('Wavelength (A)')
            plt.ylabel('Brightness')


    plt.tight_layout()
    plt.show()
