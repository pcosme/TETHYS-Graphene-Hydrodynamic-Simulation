# import required libraries
import numpy as np
import matplotlib.pyplot as plt


# Read file
file_name = "preview_1D_S=40.00vF=10.00vis=0.00l=0.00.dat"
data = np.loadtxt(file_name)
time = data[:, 0]
denL = data[:, 1]
velL = data[:, 2]
den0 = data[:, 3]
vel0 = data[:, 4]

dt = time[1]-time[0]
fs = int(1.0/dt)
nfft = 2048
tmax = time[-1]

#spectrum, freqs, im_spectrum = plt.magnitude_spectrum(denL-np.mean(denL),Fs=fs,scale='linear')
#Pxx_psd, freqs_psd = plt.psd(denL-np.mean(denL),NFFT=nfft,noverlap=nfft/2,Fs=fs,scale_by_freq=False)

fig, (ax1,ax2) = plt.subplots(nrows=2)
ax1.plot(time, denL)
ax1.set_xlim(0.0, tmax)
ax1.set_ylabel('n(L)/n0')
Pxx, freqs, bins, im = ax2.specgram(denL-np.mean(denL), Fs=fs, NFFT=nfft, noverlap=nfft/2, mode='psd', scale='dB', detrend='linear')
ax2.set_xlim(0.0, tmax)
ax2.set_xlabel('time')
ax2.set_ylabel('freq.')
#plt.show()

fig_name = "fig.jpeg"
plt.savefig(fig_name, quality=60)
plt.close("all")


#from PyEMD import EMD, Visualisation
#emd = EMD()
#emd.emd(denL)
#imfs, res = emd.get_imfs_and_residue()
#vis = Visualisation()
#vis.plot_imfs(imfs=imfs, residue=res, t=time, include_residue=True)
#vis.plot_instant_freq(time, imfs=imfs)
#vis.show()
