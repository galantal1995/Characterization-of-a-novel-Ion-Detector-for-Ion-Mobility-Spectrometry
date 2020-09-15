# If you would like to use sensitivity, make the sensitivity variable equal to True:
sensitivity = False   # With one row measuerements leave this variable False.

import numpy as np
import os
import matplotlib.pyplot as plt
import math
from tkinter import *
import tkinter.filedialog
import matplotlib


Iref = 0.586e-9
given_voltage = 6  # in kV
intt = 1700  # in microsec

onloop = 2
loops_per_on = 10 # 10 loops on 10 loops off

#Opening window
#############################################################################################################
def charges(Apix):
    given_voltage = 5.5  # in kV
    e = 1.6e-19  # Coluomb
    Q = 1 / e
    Aref = 1.69e6  # microm2
    Ipix = Apix * (Iref / Aref)  # nA
    e_per_sec = Q * Ipix
    electrons = e_per_sec * (intt * 10 ** -6)  # in one integration time on one pixel
    return electrons

number_charges = np.zeros(shape=(64,64))
number_charges[:,0:40]=charges(11*11)
number_charges[0:32,40:64]=charges(13*13)
number_charges[32:64,32:64]=charges(15*15)

class path_of_data():
    def ask_for_path(self):
        self.file = tkinter.filedialog.askopenfilename(parent=my_window)
        my_window.destroy()

    def directory(self):
        return self.file

path = path_of_data()
my_window = Tk()
my_window.title("Measurement")
my_window.geometry('500x300+300+300')
Label(my_window, text='Select the path of the measurement data: ', font=3).grid(row=1, column=1)

Button(my_window, text='Print path', command=path.ask_for_path, width=8, font=3).grid(row=1, column=2, padx=3, pady=3)

my_window.grid_rowconfigure(0, weight=1)
my_window.grid_rowconfigure(3, weight=1)
my_window.grid_columnconfigure(0, weight=1)
my_window.grid_columnconfigure(3, weight=1)

my_window.mainloop()

filename = path.directory()
npzfile = np.load(filename)

###########################################################################################################
# Input values, loaded fro the zip file that contains the measured data
v_count = 19     # uV pnp.array per count
t_count = 6.25   # us per count for integration time

series_timestamp = npzfile['time_value_list']
series_cds_rcf   = npzfile['series_cds_rcf'] * v_count
dark_frame       = npzfile['dark_frame'] * v_count
# These op thing are avoidable
num_loops = npzfile['loops']
num_frames = npzfile['frames_to_load'] - 1
intt = int(npzfile['tint'] * t_count) # in microsec
startrow = npzfile['startrow']
numrows = npzfile['numrows']
numcols = 64

# Dark frame correction using
for frame in range(0,num_loops * num_frames):
   series_cds_rcf[:,:,frame] = series_cds_rcf[:,:,frame] - dark_frame[:,:]

#########################################################
# Used variables
average_rcf = np.zeros(shape=(numrows,numcols,num_loops))   # Creating the arrays with the average values
average_noise = np.zeros(shape=(numrows, numcols, onloop))
average_signal = np.zeros(shape=(numrows, numcols, onloop))
average_noise_each_loop = np.zeros(shape=(numrows, numcols))
average_signal_each_loop = np.zeros(shape=(numrows, numcols))
delta_V =  np.zeros(shape=(numrows, numcols))    #Average of signal minus average of noise (absolute value)
##########################################################

for h in range (numrows):
    for i in range (num_loops):
        for j in range (numcols):
            average_rcf[h,j,i] = np.average(series_cds_rcf[h,j,i*num_frames:(i+1)*num_frames])

for h in range(numrows):
    for i in range (numcols):
        for j in range (onloop):
            average_noise[h,i,j] = np.average(average_rcf[h,i,j*loops_per_on*2:j*loops_per_on*2+loops_per_on])

for h in range(numrows):
    for i in range (numcols):
        for j in range (onloop):
            average_signal[h,i,j] = np.average(average_rcf[h,i,j*loops_per_on*2+loops_per_on:(j+1)*loops_per_on*2])

for h in range(numrows):
    for i in range(numcols):
        average_noise_each_loop[h,i] = np.average(average_noise[h,i,:])

for h in range(numrows):
    for i in range(numcols):
        average_signal_each_loop[h,i] = np.average(average_signal[h,i,:])

for i in range(numrows):
    for j in range(numcols):
        average_noise_each_loop[i,j] = np.average(series_cds_rcf[i,j,0:1000])
        average_signal_each_loop[i,j] = np.average(series_cds_rcf[i,j,1101:1901])

for h in range(numrows):
    for i in range (numcols):
        delta_V[h,i] = average_signal_each_loop[h,i] - average_noise_each_loop[h,i]

if sensitivity == True:
    delta_V = delta_V/number_charges

if numrows>60:  # this is just to delete the last row when the full chip is displayed.
    for i in range(numcols):
        delta_V[63,i]=0

def onclick(event):
    try:
        x_axis = int(round((event.xdata)))-1
        y_axis = int(round((event.ydata)))-startrow-1
        y_axis_plt = int(numrows-1-y_axis)
        x_axis_plot = np.arange(1, num_loops + 1)
        y_axis_plot = average_rcf[y_axis,x_axis,:]
        fig2 = plt.figure(2)
        plt.clf()

        ax = fig2.add_subplot(1, 1, 1)
        plt.plot(np.arange(1, np.size(series_cds_rcf, 2) + 1), series_cds_rcf[y_axis_plt, x_axis, :])
        ax.set_title("Full spectrum in row %d, column %d" %(y_axis+1+startrow, x_axis+1), fontsize=16)
        ax.set_xlabel("Number of frames", fontsize=14)
        ax.set_ylabel("Signal [μV]", fontsize=14)
        noise_line = ax.axhline(average_noise_each_loop[y_axis_plt, x_axis], xmin=0, xmax=num_loops, color='green',
                                 linestyle="--", label='Avr. Noise')
        signal_line = ax.axhline(average_signal_each_loop[y_axis_plt, x_axis], xmin=0, xmax=num_loops, color='red',
                                  linestyle="--", label='Avr. (Sig+Noise)')
        delta_V_line = ax.plot([0, 0], [average_noise_each_loop[y_axis_plt, x_axis], average_signal_each_loop[y_axis_plt, x_axis]],
                                color='black', linestyle="--", label='Avr. (Sig+Noise) - Avr. Noise')

        plt.text(1100, average_signal_each_loop[y_axis_plt, x_axis]-5000, '%d ions'%(int(number_charges[y_axis_plt, x_axis])), fontsize=16, fontweight='bold', color='red')
        plt.text(3100, average_signal_each_loop[y_axis_plt, x_axis]-5000, '%d ions'%(int(number_charges[y_axis_plt, x_axis])), fontsize=16, fontweight='bold', color='red')
        plt.legend(loc='center right')
        plt.grid(True)

        plt.show()
    except:
         print('Please, click only on the color map.')

fig = plt.figure(1)
my_cmap = matplotlib.cm.get_cmap('YlOrRd')
#my_cmap.set_under('gray')
#my_cmap.set_over('black')
plt.imshow(delta_V, cmap=my_cmap, extent=[.5,64.5,startrow+.5,startrow+numrows+.5], aspect='auto')
clb = plt.colorbar(mappable=None, cax=None, ax=None)
clb.set_label('Average signal - Average noise [μV]', fontsize=16)
if sensitivity == True:
    clb.set_label('Sensitivity [μV/charge]', fontsize=16)  # IF SENSITIVITY IS DESIRED TO BE CALCULATED, USE THIS LINE
plt.xlabel('Number of Columns', fontsize=16)
plt.ylabel('Number of Rows', fontsize=16)

plt.text(9, 16, 'T1\nSm', fontsize=14, fontweight='bold', color='black')
plt.text(17, 16, 'T2\nSm', fontsize=14, fontweight='bold', color='black')
plt.text(25, 16, 'T3\nSm', fontsize=14, fontweight='bold', color='black')
plt.text(33, 16, 'T4\nLrg', fontsize=14, fontweight='bold', color='black')
plt.text(41, 16, 'T5\nLrg', fontsize=14, fontweight='bold', color='black')
plt.text(49, 16, 'T6\nLrg', fontsize=14, fontweight='bold', color='black')
plt.text(57, 16, 'T6\nLrg', fontsize=14, fontweight='bold', color='black')

plt.text(9, 48, 'T7\nSm', fontsize=14, fontweight='bold', color='black')
plt.text(17, 48, 'T8\nSm', fontsize=14, fontweight='bold', color='black')
plt.text(25, 48, 'T9\nSm', fontsize=14, fontweight='bold', color='black')
plt.text(33, 48, 'T10\nSm', fontsize=14, fontweight='bold', color='black')
plt.text(41, 48, 'T11\nMed', fontsize=14, fontweight='bold', color='black')
plt.text(49, 48, 'T12\nMed', fontsize=14, fontweight='bold', color='black')
plt.text(57, 48, 'T12\nMed', fontsize=14, fontweight='bold', color='black')

plt.connect('button_press_event', onclick)  # Here I connect the onclick function to this colormap.
wm = plt.get_current_fig_manager()  # To display the figure in maximum size.
wm.window.state('zoomed')

plt.show()
