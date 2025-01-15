# This script is used to capture many important data from each xyz trajectory

import numpy as np
from scipy import signal
import os, sys
from itertools import islice
from scipy import fftpack
import glob
import time
import MDAnalysis as mda    

def smooth(x,window_len=11,window='hamming'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

class MD_Analysis:
    def __init__(self, xyz_filename, dtfs = 2, nframe_max=114514):
        self.xyz_filename = xyz_filename
        self.dtfs = dtfs
        self.dtau = dtfs * 1e-15 / 2.418884326e-17
        self.nframe_max = nframe_max
        self.natoms = 0
        self.labels = []
        self.traj = []
        self.nmolecule = 0
        self.load_xyz(self.xyz_filename)
        # After read xyz file, now we calculate different properties
        
    def load_xyz(self, xyz_filename):
        time_start = time.time()
        with open(xyz_filename, 'r') as myfile:
            self.natoms = int(myfile.readline().strip()) - 2
        self.nmolecule = self.natoms // 5
        try:
            print("Try to load file from %s.npy" %xyz_filename)
            self.traj = np.load("%s.npy" %xyz_filename)
            print(np.shape(self.traj))
            print("Loaded file from %s.npy" %xyz_filename)
        except:
            data_raw = mda.Universe(xyz_filename)
            traj = data_raw.trajectory
            nframes = len(traj)
            natoms = np.shape(np.array(traj[0]))[0]
            trajectory_data = np.zeros((natoms, 3, nframes))
            for ts in range(nframes):
                #trajectory_data[:,:,ts] = np.array(traj[ts])
                trajectory_data[:,:,ts] = traj[ts]._pos.copy()
            print("There are %d frames in %s" %(nframes, xyz_filename))
            self.traj = trajectory_data[:self.natoms,:,:]
            print(np.shape(self.traj))
            np.save(xyz_filename + ".npy", self.traj)
        time_end = time.time()
        print(f'read traj cost time = {time_end-time_start:.2f} s')

    def cacl_dipole_traj(self, which_molecule=None):
        print("Calculating dipole of water")
        cC, cH = -0.24, 0.06
        if which_molecule == None:
            C_traj  = self.traj[0:self.nmolecule*5:5, :, :]
            H1_traj = self.traj[1:self.nmolecule*5:5, :, :]
            H2_traj = self.traj[2:self.nmolecule*5:5, :, :]
            H3_traj = self.traj[3:self.nmolecule*5:5, :, :]
            H4_traj = self.traj[4:self.nmolecule*5:5, :, :]
        else:
            C_traj  = self.traj[5*which_molecule, :, :]
            H1_traj = self.traj[5*which_molecule+1, :, :]
            H2_traj = self.traj[5*which_molecule+2, :, :]
            H3_traj = self.traj[5*which_molecule+3, :, :]
            H4_traj = self.traj[5*which_molecule+4, :, :]
        self.dipole_water_traj = np.sum(C_traj * cC + H1_traj * cH + H2_traj * cH + H3_traj * cH + H4_traj * cH, axis=0)
        print("Calculating the autocorrelation function...")
        self.dacf_x = self.auto_correlation_function_simple(self.dipole_water_traj[0,:])
        self.dacf_y = self.auto_correlation_function_simple(self.dipole_water_traj[1,:])
        self.dacf_z = self.auto_correlation_function_simple(self.dipole_water_traj[2,:])
        self.dacf_tot = self.dacf_x + self.dacf_y + self.dacf_z
        self.dacf_time_fs = np.linspace(0.0, self.dtfs*(self.dacf_x.size -1), self.dacf_x.size)
        print("Calculating the FFT of autocorrelation function")
        self.dacf_x_freq, self.dacf_x_sp = self.fft3(self.dacf_x)
        self.dacf_y_freq, self.dacf_y_sp = self.fft3(self.dacf_y)
        self.dacf_z_freq, self.dacf_z_sp = self.fft3(self.dacf_z)
        self.dacf_tot_freq, self.dacf_tot_sp = self.fft3(self.dacf_tot)

    def cacl_angle_dynamic_traj(self, which_molecule=None, begintime=0, timelength=5):
        print("Calculating angles of ch4")
        length = int((1000 * timelength) / self.dtfs)
        begin  = int((1000 * begintime) / self.dtfs)
        
        if which_molecule == None:
            C_traj  = self.traj[0:self.nmolecule*5:5, :, begin:begin+length]
            H1_traj = self.traj[1:self.nmolecule*5:5, :, begin:begin+length]
            H2_traj = self.traj[2:self.nmolecule*5:5, :, begin:begin+length]
            H3_traj = self.traj[3:self.nmolecule*5:5, :, begin:begin+length]
            H4_traj = self.traj[4:self.nmolecule*5:5, :, begin:begin+length]
        else:
            C_traj  = self.traj[5*which_molecule, :, begin:begin+length]
            H1_traj = self.traj[5*which_molecule+1, :, begin:begin+length]
            H2_traj = self.traj[5*which_molecule+2, :, begin:begin+length]
            H3_traj = self.traj[5*which_molecule+3, :, begin:begin+length]
            H4_traj = self.traj[5*which_molecule+4, :, begin:begin+length]

        def get_angle(C_traj, H1_traj, H2_traj):
            CH1 = H1_traj - C_traj
            CH2 = H2_traj - C_traj
            cross = np.einsum('ijk,ijk->ik', CH1, CH2, optimize=True)
            norm1 = np.sqrt(np.sum(np.abs(CH1)**2, axis=1))
            norm2 = np.sqrt(np.sum(np.abs(CH2)**2, axis=1))
            cosangle = np.einsum('ij,ij,ij->ij', cross, 1/norm1, 1/norm2, optimize=True)
            angle = np.arccos(cosangle)
            return angle
        
        shape = np.shape(C_traj)
        self.angle_ch4_traj = np.zeros((6,shape[0],shape[2]))
        self.angle_ch4_traj[0,:,:] = get_angle(C_traj, H1_traj, H2_traj)
        self.angle_ch4_traj[1,:,:] = get_angle(C_traj, H1_traj, H3_traj)
        self.angle_ch4_traj[2,:,:] = get_angle(C_traj, H1_traj, H4_traj)
        self.angle_ch4_traj[3,:,:] = get_angle(C_traj, H2_traj, H3_traj)
        self.angle_ch4_traj[4,:,:] = get_angle(C_traj, H2_traj, H4_traj)
        self.angle_ch4_traj[5,:,:] = get_angle(C_traj, H3_traj, H4_traj)
        
        time_aac_start = time.time()
        print("Calculating the autocorrelation function...")
        aacf_dynamic = np.zeros((6,shape[0],shape[2] // 2))
        for imole in range(shape[0]):
            for iangle in range(6):
                aacf_dynamic[iangle,imole,:] = self.auto_correlation_function_simple(self.angle_ch4_traj[iangle,imole,:])
        self.aacf_dynamic = aacf_dynamic
        self.aacf_dynamic_time_fs = np.linspace(0.0, self.dtfs*(shape[2] // 2 - 1), shape[2] // 2)
        time_aac_end = time.time()
        print(f'calculation the autocorrelation function cost time = {time_aac_end-time_aac_start:.2f} s')
        
        time_fft_start = time.time()
        print("Calculating the FFT of autocorrelation function")
        fft_aacf_dynamic_sp = np.zeros((6,shape[0],shape[2] // 2))
        for imole in range(shape[0]):
            for iangle in range(6):
                self.aacf_dynamic_tot_freq, fft_aacf_dynamic_sp[iangle,imole,:] = self.fft3(self.aacf_dynamic[iangle,imole,:])

        mole_dynamic_sp = np.mean(fft_aacf_dynamic_sp, axis=0)
        aacf_dynamic_tot_sp = None
        if aacf_dynamic_tot_sp == None : aacf_dynamic_tot_sp = np.mean(mole_dynamic_sp, axis=0)
        time_fft_end = time.time()
        print(f'calculation the autocorrelation function cost time = {time_fft_end-time_fft_start:.2f} s')
        return aacf_dynamic_tot_sp

    def cacl_raman_traj(self, which_molecule=None):
        print("Calculating angles of ch4")
        if which_molecule == None:
            C_traj  = self.traj[0:self.nmolecule*5:5, :, :]
            H1_traj = self.traj[1:self.nmolecule*5:5, :, :]
            H2_traj = self.traj[2:self.nmolecule*5:5, :, :]
            H3_traj = self.traj[3:self.nmolecule*5:5, :, :]
            H4_traj = self.traj[4:self.nmolecule*5:5, :, :]
        else:
            C_traj  = self.traj[5*which_molecule, :, :]
            H1_traj = self.traj[5*which_molecule+1, :, :]
            H2_traj = self.traj[5*which_molecule+2, :, :]
            H3_traj = self.traj[5*which_molecule+3, :, :]
            H4_traj = self.traj[5*which_molecule+4, :, :]

        def get_angle(C_traj, H1_traj, H2_traj):
            CH1 = H1_traj - C_traj
            CH2 = H2_traj - C_traj
            cross = np.einsum('ijk,ijk->ik', CH1, CH2, optimize=True)
            norm1 = np.sqrt(np.sum(np.abs(CH1)**2, axis=1))
            norm2 = np.sqrt(np.sum(np.abs(CH2)**2, axis=1))
            cosangle = np.einsum('ij,ij,ij->ij', cross, 1/norm1, 1/norm2, optimize=True)
            angle = np.arccos(cosangle)
            return angle - (107.6 / 180) * 3.14159265359
        
        def get_length(C_traj, H_traj):
            CH = H_traj - C_traj
            norm = np.sqrt(np.sum(np.abs(CH)**2, axis=1))
            return norm - 1.101
        
        shape = np.shape(C_traj)
        self.angle_ch4_traj = np.zeros((6,shape[0],shape[2]))
        self.angle_ch4_traj[0,:,:] = get_angle(C_traj, H1_traj, H2_traj)
        self.angle_ch4_traj[1,:,:] = get_angle(C_traj, H1_traj, H3_traj)
        self.angle_ch4_traj[2,:,:] = get_angle(C_traj, H1_traj, H4_traj)
        self.angle_ch4_traj[3,:,:] = get_angle(C_traj, H2_traj, H3_traj)
        self.angle_ch4_traj[4,:,:] = get_angle(C_traj, H2_traj, H4_traj)
        self.angle_ch4_traj[5,:,:] = get_angle(C_traj, H3_traj, H4_traj)
        
        self.length_ch4_traj = np.zeros((4,shape[0],shape[2]))
        self.length_ch4_traj[0,:,:] = get_length(C_traj, H1_traj)
        self.length_ch4_traj[1,:,:] = get_length(C_traj, H2_traj)
        self.length_ch4_traj[2,:,:] = get_length(C_traj, H3_traj)
        self.length_ch4_traj[3,:,:] = get_length(C_traj, H4_traj)

        def get_s2a(angle_traj):
            s2a_traj  = np.zeros_like(angle_traj[0,:,:])
            s2a_traj += 2 * angle_traj[0,:,:]
            s2a_traj -= angle_traj[1,:,:]
            s2a_traj -= angle_traj[2,:,:]
            s2a_traj -= angle_traj[3,:,:]
            s2a_traj -= angle_traj[4,:,:]
            s2a_traj += 2 * angle_traj[5,:,:]
            return s2a_traj * (12 ** (-0.5))
        
        def get_s2b(angle_traj):
            s2b_traj  = np.zeros_like(angle_traj[0,:,:])
            s2b_traj += angle_traj[1,:,:]
            s2b_traj -= angle_traj[2,:,:]
            s2b_traj -= angle_traj[3,:,:]
            s2b_traj += angle_traj[4,:,:]
            return s2b_traj * 0.5
        
        def get_s3x(length_traj):
            s3x_traj  = np.zeros_like(length_traj[0,:,:])
            s3x_traj += length_traj[0,:,:]
            s3x_traj -= length_traj[1,:,:]
            s3x_traj += length_traj[2,:,:]
            s3x_traj -= length_traj[3,:,:]
            return s3x_traj * 0.5
        
        def get_s3y(length_traj):
            s3y_traj  = np.zeros_like(length_traj[0,:,:])
            s3y_traj += length_traj[0,:,:]
            s3y_traj -= length_traj[1,:,:]
            s3y_traj -= length_traj[2,:,:]
            s3y_traj += length_traj[3,:,:]
            return s3y_traj * 0.5
        
        def get_s3z(length_traj):
            s3z_traj  = np.zeros_like(length_traj[0,:,:])
            s3z_traj += length_traj[0,:,:]
            s3z_traj += length_traj[1,:,:]
            s3z_traj -= length_traj[2,:,:]
            s3z_traj -= length_traj[3,:,:]
            return s3z_traj * 0.5
        
        def get_s4x(angle_traj):
            s4x_traj  = np.zeros_like(angle_traj[0,:,:])
            s4x_traj += angle_traj[4,:,:]
            s4x_traj -= angle_traj[1,:,:]
            return s4x_traj * (2 ** (-0.5))
        
        def get_s4y(angle_traj):
            s4y_traj  = np.zeros_like(angle_traj[0,:,:])
            s4y_traj += angle_traj[3,:,:]
            s4y_traj -= angle_traj[2,:,:]
            return s4y_traj * (2 ** (-0.5))
        
        def get_s4z(angle_traj):
            s4z_traj  = np.zeros_like(angle_traj[0,:,:])
            s4z_traj += angle_traj[5,:,:]
            s4z_traj -= angle_traj[0,:,:]
            return s4z_traj * (2 ** (-0.5))

        self.s1_ch4_traj  = 0.5 * np.sum(self.length_ch4_traj, axis=0)
        self.s2a_ch4_traj = get_s2a(self.angle_ch4_traj)
        self.s2b_ch4_traj = get_s2b(self.angle_ch4_traj)
        self.s3x_ch4_traj = get_s3x(self.length_ch4_traj)
        self.s3y_ch4_traj = get_s3y(self.length_ch4_traj)
        self.s3z_ch4_traj = get_s3z(self.length_ch4_traj)
        self.s4x_ch4_traj = get_s4x(self.angle_ch4_traj)
        self.s4y_ch4_traj = get_s4y(self.angle_ch4_traj)
        self.s4z_ch4_traj = get_s4z(self.angle_ch4_traj)

        # Molecular Physics, 64:6, 1061-1071, DOI: 10.1080/00268978800100713
        pe  = 14.945
        p1  = 13.802 # in angstrom^(-1)
        p11 = 11.9   # in angstrom^(-2)
        p22 = 2.42   # in rad^(-2)
        p33 = 11.09  # in angstrom^(-2)
        p44 = 1.5    # in rad^(-2)
        p34 = 0.85   # in angstrom^(-1) * rad^(-1)

        time_aac_start = time.time()
        print("Calculating the autocorrelation function...")
        pacf_s1  = np.zeros((shape[0],shape[2] // 2))
        pacf_s2a = np.zeros((shape[0],shape[2] // 2))
        pacf_s2b = np.zeros((shape[0],shape[2] // 2))
        pacf_s3x = np.zeros((shape[0],shape[2] // 2))
        pacf_s3y = np.zeros((shape[0],shape[2] // 2))
        pacf_s3z = np.zeros((shape[0],shape[2] // 2))
        pacf_s4x = np.zeros((shape[0],shape[2] // 2))
        pacf_s4y = np.zeros((shape[0],shape[2] // 2))
        pacf_s4z = np.zeros((shape[0],shape[2] // 2))
        for imole in range(shape[0]):
            pacf_s1[imole,:]  = self.auto_correlation_function_simple(self.s1_ch4_traj[imole,:])
            pacf_s2a[imole,:] = self.auto_correlation_function_simple(self.s2a_ch4_traj[imole,:])
            pacf_s2b[imole,:] = self.auto_correlation_function_simple(self.s2b_ch4_traj[imole,:])
            pacf_s3x[imole,:] = self.auto_correlation_function_simple(self.s3x_ch4_traj[imole,:])
            pacf_s3y[imole,:] = self.auto_correlation_function_simple(self.s3y_ch4_traj[imole,:])
            pacf_s3z[imole,:] = self.auto_correlation_function_simple(self.s3z_ch4_traj[imole,:])
            pacf_s4x[imole,:] = self.auto_correlation_function_simple(self.s4x_ch4_traj[imole,:])
            pacf_s4y[imole,:] = self.auto_correlation_function_simple(self.s4y_ch4_traj[imole,:])
            pacf_s4z[imole,:] = self.auto_correlation_function_simple(self.s4z_ch4_traj[imole,:])
            
        self.pacf_time_fs = np.linspace(0.0, self.dtfs*(shape[2] // 2 - 1), shape[2] // 2)
        time_aac_end = time.time()
        print(f'calculation the autocorrelation function cost time = {time_aac_end-time_aac_start:.2f} s')
        
        time_fft_start = time.time()
        print("Calculating the FFT of autocorrelation function")
        fft_pacf_s1_sp  = np.zeros((shape[0],shape[2] // 2))
        fft_pacf_s2a_sp = np.zeros((shape[0],shape[2] // 2))
        fft_pacf_s2b_sp = np.zeros((shape[0],shape[2] // 2))
        fft_pacf_s3x_sp = np.zeros((shape[0],shape[2] // 2))
        fft_pacf_s3y_sp = np.zeros((shape[0],shape[2] // 2))
        fft_pacf_s3z_sp = np.zeros((shape[0],shape[2] // 2))
        fft_pacf_s4x_sp = np.zeros((shape[0],shape[2] // 2))
        fft_pacf_s4y_sp = np.zeros((shape[0],shape[2] // 2))
        fft_pacf_s4z_sp = np.zeros((shape[0],shape[2] // 2))
        for imole in range(shape[0]):
            self.pacf_tot_freq, fft_pacf_s1_sp[imole,:] = self.fft3(pacf_s1[imole,:])
            _, fft_pacf_s2a_sp[imole,:]  = self.fft3(pacf_s2a[imole,:])
            _, fft_pacf_s2b_sp[imole,:]  = self.fft3(pacf_s2b[imole,:])
            _, fft_pacf_s3x_sp[imole,:]  = self.fft3(pacf_s3x[imole,:])
            _, fft_pacf_s3y_sp[imole,:]  = self.fft3(pacf_s3y[imole,:])
            _, fft_pacf_s3z_sp[imole,:]  = self.fft3(pacf_s3z[imole,:])
            _, fft_pacf_s4x_sp[imole,:]  = self.fft3(pacf_s4x[imole,:])
            _, fft_pacf_s4y_sp[imole,:]  = self.fft3(pacf_s4y[imole,:])
            _, fft_pacf_s4z_sp[imole,:]  = self.fft3(pacf_s4z[imole,:])

        self.pacf_s1_sp  = np.mean(fft_pacf_s1_sp, axis=0)
        self.pacf_s2a_sp = np.mean(fft_pacf_s2a_sp, axis=0)
        self.pacf_s2b_sp = np.mean(fft_pacf_s2b_sp, axis=0)
        self.pacf_s3x_sp = np.mean(fft_pacf_s3x_sp, axis=0)
        self.pacf_s3y_sp = np.mean(fft_pacf_s3y_sp, axis=0)
        self.pacf_s3z_sp = np.mean(fft_pacf_s3z_sp, axis=0)
        self.pacf_s4x_sp = np.mean(fft_pacf_s4x_sp, axis=0)
        self.pacf_s4y_sp = np.mean(fft_pacf_s4y_sp, axis=0)
        self.pacf_s4z_sp = np.mean(fft_pacf_s4z_sp, axis=0)
        
        time_fft_end = time.time()
        print(f'calculation the autocorrelation function cost time = {time_fft_end-time_fft_start:.2f} s')

    def auto_correlation_function_fft(self, x):
        corr = signal.fftconvolve(x, x[::-1], mode='same')
        corr = corr[corr.size // 2: ]
        return corr / corr[0] * np.mean(x * x)

    def auto_correlation_function_simple(self, x):
        n = x.size
        if n % 2 == 0:
            x_shifted = np.zeros(n*2)
        else:
            x_shifted = np.zeros(n*2-1)
        x_shifted[n//2 : n//2+n] = x
        # Convolute the shifted array with the flipped array, which is equivalent to performing a correlation
        autocorr_full = (signal.fftconvolve(x_shifted, x[::-1], mode='same')[-n:]/ np.arange(n, 0, -1))
        # Truncate the autocorrelation array
        autocorr = autocorr_full[0:n//2]
        return autocorr

    def fft(self, x):
        #sp = np.fft.fft(x)
        #freq_au = 2.0 * np.pi * np.fft.fftfreq(np.size(x), self.dtau)
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        #freq_cminverse = freq_au * 219474.63
        #return freq_cminverse[0:sp.size//2], sp[0:sp.size//2]
        lineshape = fftpack.dct(x, type=1)
        freq_au = np.linspace(0, 0.5/self.dtfs * 1e15, len(x))
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        freq_cminverse = freq_au / (100.0 * 299792458.0)
        return freq_cminverse, lineshape

    def fft3(self, x ):
        # Adding zeros to the end of x
        lineshape = fftpack.dct(x, type=1)
        freq_au = np.linspace(0, 0.5/self.dtfs * 1e15, len(x))
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        freq_cminverse = freq_au / (100.0 * 299792458.0)
        # Calculate spectra
        #field_description =  freq_au**2
        field_description =  freq_au**2
        spectra = lineshape * field_description
        return freq_cminverse, spectra
        #return freq_cminverse[0:spectra.size//2], spectra[0:spectra.size//2]

    def output_dipole_autocorrelation(self, which_molecule=None):
        if which_molecule is None:
            local_filename = "%s.dac.txt" %self.xyz_filename
        else:
            local_filename = "%s.dac_%d.txt" %(self.xyz_filename, which_molecule)
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated dipole autocorrelation for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_dipole_traj(which_molecule=which_molecule)
            # output data
            data = np.zeros((np.size(self.dacf_x), 10))
            data[:, 0] = self.dacf_time_fs
            data[:, 1] = self.dacf_x
            data[:, 2] = self.dacf_y
            data[:, 3] = self.dacf_z
            data[:, 4] = self.dacf_tot
            data[:, 5] = self.dacf_x_freq
            data[:, 6] = smooth(self.dacf_x_sp)
            data[:, 7] = smooth(self.dacf_y_sp)
            data[:, 8] = smooth(self.dacf_z_sp)
            data[:, 9] = smooth(self.dacf_tot_sp)
            comments = "# dacf_time_fs, dacf_x, dacf_y, dacf_z, dacf_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments)
    
    def output_dynamic_angle_autocorrelation(self, which_molecule=None, simulation_time=20, timelength=5):
        for begintime in range(simulation_time-timelength+1):
            print(f'calculating angle autocorrelation and fft within time interval [{begintime},{begintime+timelength}]')
            if which_molecule is None:
                local_filename = "%s.aac_%d.txt" % (self.xyz_filename, begintime)
            else:
                local_filename = "%s.aac_%d_%d.txt" %(self.xyz_filename, begintime, which_molecule)
            if os.path.isfile(local_filename + "dd"):
                print("Have calculated angle autocorrelation for %s, skipping..." %self.xyz_filename)
            else:
                aacf_dynamic_tot_sp = self.cacl_angle_dynamic_traj(which_molecule=which_molecule, begintime=begintime, timelength=timelength)
                # output data
                aac_len = np.shape(aacf_dynamic_tot_sp)[0]
                data = np.zeros((aac_len, 4))
                data[:, 0] = self.aacf_dynamic_time_fs
                data[:, 1] = np.mean(np.mean(self.aacf_dynamic,axis=0),axis=0)
                data[:, 2] = self.aacf_dynamic_tot_freq
                data[:, 3] = smooth(aacf_dynamic_tot_sp)
                comments = "# dacf_time_fs, aacf_tot, freq, sp_tot"
                np.savetxt(local_filename, data, comments=comments)

    def output_symmetry_coordinate_autocorrelation(self, which_molecule=None):
        if which_molecule is None: 
            local_filename = "%s.pac.txt" %self.xyz_filename
            coord_filename = "%s.scoord.txt" %self.xyz_filename
        else: 
            local_filename = "%s.pac_%d.txt" %(self.xyz_filename, which_molecule)
            coord_filename = "%s.scoord_%d.txt" %(self.xyz_filename, which_molecule)
        if os.path.isfile(coord_filename + "dd"): print("Have calculated angle autocorrelation for %s, skipping..." %self.xyz_filename)
        else: 
            self.cacl_raman_traj(which_molecule=which_molecule)
            
            pac = np.zeros((len(self.pacf_tot_freq) ,10))
            pac[:, 0] = smooth(self.pacf_s1_sp)
            pac[:, 1] = smooth(self.pacf_s2a_sp)
            pac[:, 2] = smooth(self.pacf_s2b_sp)
            pac[:, 3] = smooth(self.pacf_s3x_sp)
            pac[:, 4] = smooth(self.pacf_s3y_sp)
            pac[:, 5] = smooth(self.pacf_s3z_sp)
            pac[:, 6] = smooth(self.pacf_s4x_sp)
            pac[:, 7] = smooth(self.pacf_s4y_sp)
            pac[:, 8] = smooth(self.pacf_s4z_sp)
            pac[:, 9] = self.pacf_tot_freq
            np.savetxt(local_filename, pac)

            def get_dot(traj):
                return np.einsum('ij,ij->ij', traj, traj)
            
            coord = np.zeros((np.shape(self.s1_ch4_traj)[1] ,9))
            coord[:, 0] = np.mean(get_dot(self.s1_ch4_traj), axis=0)
            coord[:, 1] = np.mean(get_dot(self.s2a_ch4_traj), axis=0)
            coord[:, 2] = np.mean(get_dot(self.s2b_ch4_traj), axis=0)
            coord[:, 3] = np.mean(get_dot(self.s3x_ch4_traj), axis=0)
            coord[:, 4] = np.mean(get_dot(self.s3y_ch4_traj), axis=0)
            coord[:, 5] = np.mean(get_dot(self.s3z_ch4_traj), axis=0)
            coord[:, 6] = np.mean(get_dot(self.s4x_ch4_traj), axis=0)
            coord[:, 7] = np.mean(get_dot(self.s4y_ch4_traj), axis=0)
            coord[:, 8] = np.mean(get_dot(self.s4z_ch4_traj), axis=0)
            np.savetxt(coord_filename, coord)

if __name__ == "__main__":
    paths = sys.argv[1:]
    for i, path in enumerate(paths):
        filenames = glob.glob("%s/simu_*.xc.xyz" %path)
        print(len(filenames))
        for filename in filenames:
            a = MD_Analysis(xyz_filename=filename)
            #a.output_dipole_autocorrelation(which_molecule=None)
            #a.output_dynamic_angle_autocorrelation(which_molecule=None, simulation_time=20, timelength=5)
            a.output_symmetry_coordinate_autocorrelation(which_molecule=None)
