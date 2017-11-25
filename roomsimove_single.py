'''
ROOMSIMOVE_SINGLE Compute shoebox room filters for a single source
###########################################################################
# Copyright 2003 Douglas R. Campbell
# Copyright 2008-2016 Emmanuel Vincent
# Copyright 2017 Sunit Sivasankaran
# This software is a python version of the stripped-down version of the Roomsim toolbox version
# 3.3 by Douglas R. Campbell ,
# The matlab function for the stripped down version can be found here:
Roomsimove, http://homepages.loria.fr/evincent/software/Roomsimove.zip
# This code is distributed under the terms of the GNU Public License version 3
# (http://www.gnu.org/licenses/gpl.txt)
# 
# If you find it useful, please cite the following reference:
###########################################################################
One  difference between the matlab version and this code is that 
RT60 value is assumed to be same for all frequencies.

Tested for sampling rate of 16000 Hz. 

Usage:
=========

As standalone file:
------------------
    python roomsimove_single.py config_file source_pos_x source_pos_y source_pos_z output_file

    The help options will also give the details
    python roomsimove_single.py -h

As a module:
------------
    using config_file
    -----------------
    import roomsimove_single
    sim_rir = roomsimove_single.RoomSim.init_from_config_file(config_file)
    source_pos = [1, 1, 1]
    rir = sim_rir.create_rir(source_pos)

    using default values of absorption coeffecients
    -----------------------------------------------
    import roomsimove_single
    rt60 = 0.5 # in seconds
    room_dim = [4.2, 3.4, 5.2] # in meters
    absorption = roomsimove_single.rt60_to_absorption(room_dim, rt60)
    room = roomsimove_single.Room(room_dim, abs_coeff=absorption)
    mic_pos = [2, 2, 2] # in  meters
    mic1 = roomsimove_single.Microphone(mic_pos, 1,  \
            orientation=[0.0, 0.0, 0.0], direction='omnidirectional')
    mic_pos = [2, 2, 1] # in  meters
    mic2 = roomsimove_single.Microphone(mic_pos, 2,  \
            orientation=[0.0, 0.0, 0.0], direction='cardioid')
    mics = [mic1, mic2]
    sampling_rate = 16000
    sim_rir = roomsimove_single.RoomSim(sampling_rate, room, mics, RT60=rt60)
    source_pos = [1, 1, 1] # in  meters
    rir = sim_rir.create_rir(source_pos)

Appyling RIR to data
-------------------
    import olafilt
    import sounfile as sf
    # Assuming single channel data
    [data, fs] = sf.read(wav_file)
    reverb_data = olafilt.olafilt(rir,data)
'''

import argparse
import numpy as np
from scipy.interpolate import interp1d
import scipy.signal as scipy_sig


def do_everything(room_dim, mic_positions, source_pos, rt60):
    absorption = rt60_to_absorption(room_dim, rt60)
    room = Room(room_dim, abs_coeff=absorption)
    mics = []
    for idx, mic in enumerate(mic_positions):
        temp_mic = Microphone(mic, idx,  \
            orientation=[0.0, 0.0, 0.0], direction='omnidirectional')
        mics.append(temp_mic)
    sim_rir = RoomSim(16000, room, mics, RT60=rt60)
    rir = sim_rir.create_rir(source_pos)
    return rir

def get_rt60(F_abs, room_size, A):
    '''
    Get RT 60 given the room characteristics
    '''
    m_air = 6.875e-4*(F_abs.T/1000)**(1.7)
    # attenuation factors for one metre travelled in air
    room_size = np.array(room_size)
    atten_air = np.exp(-0.5*m_air).T
    Lx = room_size[0]
    Ly = room_size[1]
    Lz = room_size[2]
    #Volume of room m^3
    V_room=Lx*Ly*Lz
    area_xz=Lx*Lz
    area_yz=Ly*Lz
    area_xy=Lx*Ly
    total_area = 2*(area_xz+area_yz+area_xy)# Total area of shoebox room surfaces
    # Effective absorbing area of room surfaces at each frequency
    Se=area_yz*(A[0]+A[1])+area_xz*(A[2]+A[3])+area_xy*(A[5]+A[4])
    a_bar=Se/total_area # Mean absorption of each room surface
    # Norris-Eyring estimate adjusted for air absorption
    RT60=0.1611*V_room/(4*m_air.T*V_room-total_area*np.log(1-a_bar))
    return RT60

def rt60_to_absorption(room_obj_dim, rt60):
    '''
    Norris-Eyring formula %%
     Converts a given reverberation time into a single absorption coefficient for all surfaces 
    '''
    room_vol = np.prod(room_obj_dim)
    area_xz=room_obj_dim[0] * room_obj_dim[2]
    area_yz=room_obj_dim[1] * room_obj_dim[2]
    area_xy=room_obj_dim[0] * room_obj_dim[1]
    total_area =2*(area_xz+area_yz+area_xy); # Total area of shoebox room surfaces
    absorption = 1-np.exp(-0.1611*room_vol/(total_area*rt60))
    return absorption

class Microphone(object):
    '''
        Deal with a single microphone
    '''
    def __init__(self, pos, id_val,  \
            orientation=[0.0, 0.0, 0.0], direction='omnidirectional'):
        self.x_pos = pos[0] 
        self.y_pos = pos[1]
        self.z_pos = pos[2]
        self.pos = pos
        self._id = str(id_val)
        self.orientation = orientation
        self.direction = direction

class Room(object):
    '''
    Room characteristics
    '''
    def __init__(self, dim, F_abs=None, abs_coeff=None):
        self.x_val = dim[0]
        self.y_val = dim[1]
        self.z_val = dim[2]
        self.room_size = np.array(dim)
        self.freq_dep_absorption = {}
        if F_abs is None:
            self.freq_dep_absorption['F_abs'] = np.array([125, 250, 500, 1000, 2000, 4000, 8000])
        else:
            self.freq_dep_absorption['F_abs'] = np.array(F_abs)
        if abs_coeff is None:
            self.__set_absorption()
        else:
            if isinstance(abs_coeff, float) or isinstance(abs_coeff, int):
                self.__set_absorption(abs_val=abs_coeff)
            else:
                self.freq_dep_absorption['Ax1'] = np.array(abs_coeff[0])
                self.freq_dep_absorption['Ax2'] = np.array(abs_coeff[1])
                self.freq_dep_absorption['Ay1'] = np.array(abs_coeff[2])
                self.freq_dep_absorption['Ay2'] = np.array(abs_coeff[3])
                self.freq_dep_absorption['Az1'] = np.array(abs_coeff[4])
                self.freq_dep_absorption['Az2'] = np.array(abs_coeff[5])

    def __set_absorption(self, abs_val=0.671):
        self.freq_dep_absorption['Ax1'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Ax2'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Ay1'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Ay2'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Az1'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Az2'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))


class Config(object):
    '''
    Interface to read config files and put it to the right objects
    '''
    def __init__(self, config_file):
        self._file = config_file
        self.config = {}
        with open(config_file) as fid:
            for line in fid:
                line = line.strip()
                if line.startswith('%') or line == '':
                # This is a comment. Ignore
                    continue
                temp = line.split()
                try:
                    self.config[temp[0]] = [float(temp_) for temp_ in temp[1:]]
                except:
                    self.config[temp[0]] = [temp_ for temp_ in temp[1:]]
        self.config['Fs'] = int(self.config['Fs'][0])
        dict_keys = self.config.keys()
        self.sp_keys = [ke for  ke in dict_keys if ke.startswith('sp')]
        self.sd_keys = [ke for  ke in dict_keys if ke.startswith('sd')]
        self.so_keys = [ke for  ke in dict_keys if ke.startswith('so')]
        self.__verify_config()


    def __verify_config(self):
        assert 'room_size' in self.config, 'room_size not found in config'
        assert 'F_abs' in self.config, 'F_abs not found in config'
        assert 'Ax1' in self.config, 'Ax1 not found in config'
        assert 'Ax2' in self.config, 'Ax2 not found in config'
        assert 'Ay1' in self.config, 'Ay1 not found in config'
        assert 'Ay2' in self.config, 'Ay2 not found in config'
        assert 'Az1' in self.config, 'Az1 not found in config'
        assert 'Az2' in self.config, 'Az2 not found in config'
        assert 'sp1' in self.config, 'sp1 not found in config'
        assert 'sd1' in self.config, 'sd1 not found in config'
        assert 'so1' in self.config, 'so1 not found in config'
        assert len(self.sp_keys) == len(self.sd_keys) == len(self.so_keys), \
            'sp, sd and so are not of same length'

    def create_room_et_mic_objects(self):
        room_size = [float(_) for _ in self.config['room_size']]
        F_abs = [float(_) for _ in self.config['F_abs']]
        Ax1 = [float(_) for _ in self.config['Ax1']]
        Ax2 = [float(_) for _ in self.config['Ax2']]
        Ay1 = [float(_) for _ in self.config['Ay1']]
        Ay2 = [float(_) for _ in self.config['Ay2']]
        Az1 = [float(_) for _ in self.config['Az1']]
        Az2 = [float(_) for _ in self.config['Az2']]
        room = Room(room_size, F_abs, [Ax1, Ax2, Ay1, Ay2, Az1, Az2])
        mics = []
        for mic_idx in range(len(self.sp_keys)):
            mic_idx += 1
            _xp, _yp, _zp = self.config['sp'+str(mic_idx)]
            orientation = self.config['so'+str(mic_idx)]
            direction = self.config['sd'+str(mic_idx)][0].replace("'",'')
            mics.append(Microphone([_xp, _yp, _zp], mic_idx,\
                                  orientation = orientation, direction = direction))
        return[self.config['Fs'], room, mics]

class RoomSim(object):
    '''
    Class to handle RIR creation:
        Input
        -----
        room_config : Roomconfig object

    '''

    def __init__(self, fs, room, mics, RT60=None):
        self._do_init(fs, room, mics, RT60)
        self.verify_positions()

    def verify_positions(self):
        '''
        Method to verify if all the microphones are inside the room
        '''
        for mic in self.mics:
            assert mic.x_pos < self.room.x_val,\
                    mic._id+' x position is outside the room'
            assert mic.y_pos < self.room.y_val,\
                    mic._id+' y position is outside the room'
            assert mic.z_pos < self.room.z_val,\
                    mic._id+' z position is outside the room'


    @classmethod
    def init_from_config_file(cls, room_config_file, RT60=None):
        '''
        constructor to read config file and initialize an instance
        '''
        config = Config(room_config_file)
        sample_rate, room, mics = config.create_room_et_mic_objects()
        obj = cls(sample_rate, room, mics, RT60)
        return obj

    def _do_init(self, fs, room, mics, RT60):
        self.sampling_rate = fs
        self.room = room
        self.mics = mics
        mic_count = 0
        for mic in self.mics:
            mic_count += 1
            mic._id = str(mic_count)
        self.channels = len(mics)
        self.room_size = room.room_size
        self.F_abs = room.freq_dep_absorption['F_abs']
        Ax1 = room.freq_dep_absorption['Ax1']
        Ax2 = room.freq_dep_absorption['Ax2']
        Ay1 = room.freq_dep_absorption['Ay1']
        Ay2 = room.freq_dep_absorption['Ay2']
        Az1 = room.freq_dep_absorption['Az1']
        Az2 = room.freq_dep_absorption['Az2']
        self.A = np.array([Ax1, Ax2, Ay1, Ay2, Az1, Az2])
        self.A = self.A[:, self.F_abs<=self.sampling_rate/2.0]
        self.F_abs = self.F_abs[self.F_abs<=self.sampling_rate/2.0]
        if self.F_abs[0] != 0:
            self.A = np.vstack((self.A.T[0], self.A.T)).T
            self.F_abs = np.hstack((0, self.F_abs))
        if self.F_abs[-1] != self.sampling_rate/2.0:
            self.A = np.vstack((self.A.T, self.A.T[-1]))
            self.F_abs = np.hstack((self.F_abs, self.sampling_rate/2.0))

        self.tm_sensor = np.zeros((self.channels, 3, 3))
        self.sensor_xyz = np.zeros((self.channels, 3))
        self.sensor_off = np.zeros((self.channels, 3))
        for idx, mic in enumerate(self.mics):
            self.sensor_xyz[idx, :] = mic.pos
            self.sensor_off[idx, :] = mic.orientation
            self.tm_sensor[idx, :, :] = self.__create_tm(\
                                    self.__create_psi_theta_phi(mic.orientation))
        if RT60 is None:
            self.RT60 = get_rt60(self.F_abs, self.room_size, self.A)
        else:
            self.RT60 = np.array([RT60] * len(self.F_abs))


    def create_rir(self, source_xyz, source_off=None, source_dir=None):
        '''
        Create the RIR
        source_xyz : list containing xyz position of the source
        source_off: 3 x 1 list representing the source orientation (azimuth,
        elevation, roll)
        source_dir: source directivity np txt file of dimension 181 x 361
        '''

        source_xyz = np.array(source_xyz)
        if source_dir is None:
            # omnidirectional
            source_dir = np.ones((181,361))
        else:
            source_dir = np.loadtxt(source_dir)
        if source_off is None:
            source_off = np.zeros(source_xyz.shape)

        [c_psi, s_psi, c_theta, s_theta, c_phi, s_phi] = \
            self.__create_psi_theta_phi(source_off)
        tm_source = self.__create_tm([c_psi, s_psi, c_theta, s_theta, c_phi, s_phi])
        Two_pi = 2*np.pi
        sampling_period = 1.0/self.sampling_rate
        nyquist = self.sampling_rate/2.0 # Half sampling frequency
        Fs_c = self.sampling_rate/343.0 # Samples per metre
        # Reflection order and impulse response length
        # H_length = longest reverberation time in samples (rounded down to integer)
        H_length = np.floor(np.max(self.RT60)*self.sampling_rate)
        range_ = H_length/Fs_c # H_length in metres
        Lx = self.room_size[0]
        Ly = self.room_size[1]
        Lz = self.room_size[2]
        order_x = np.ceil(range_/(2*Lx)); #  Number in +x direction
        order_y = np.ceil(range_/(2*Ly)); #  Number in +y direction
        order_z = np.ceil(range_/(2*Lz)); #  Number in +z direction
        #Maximum number of image sources
        n_isources = int((2*order_x+1)*(2*order_y+1)*(2*order_z+1)*8)
        delay_s = Fs_c*np.sqrt(np.sum((source_xyz.T-self.sensor_xyz)**2, axis=1))
        # Ensure H_length > 200 points so that a full CIPIC or MIT HRIR can be viewed
        H_length = int(np.max((H_length, np.ceil(np.max(np.max(delay_s)))+200)))
        # Interpolation filter for fractional delays
        N_frac = 32 # Order of FIR fractional delay filter
        Tw = N_frac*sampling_period # Window duration (seconds)
        Two_pi_Tw = Two_pi/Tw # Compute here for efficiency
        # Filter time window NB column vector of length (N_frac+1) symmetrical about t=0
        t = np.arange(-Tw/2, Tw/2+sampling_period, sampling_period)
        # Column vector of zero values for post-padding
        pad_frac = np.zeros((N_frac,1))
        # Second order high-pass IIR filter to remove DC buildup (nominal -4dB cut-off at 20 Hz)
        w = 2*np.pi*20
        r1 = np.exp(-w*sampling_period)
        r2 =  np.exp(-w*sampling_period)
        b1 = -(1+r2)
        b2 = np.copy(r2) #Numerator coefficients (fix zeros)
        a1 = 2*r1*np.cos(w*sampling_period)
        a2 = -r1*r1 #Denominator coefficients (fix poles)
        HP_gain = (1-b1+b2)/(1+a1-a2) #Normalisation gain
        b_HP = [1, b1, b2]/HP_gain
        a_HP = [1,-a1,-a2]
        # Further constants
        Two_Lx = 2*self.room_size[0] # Twice Length (Depth)
        Two_Ly = 2*self.room_size[1] # Twice Width
        Two_Lz = 2*self.room_size[2] # Twice Height
        #codes the eight permutations of x+/-xp, y+/-yp, z+/-zp
        #(the source to receiver vector components) where [-1 -1 -1] identifies the parent source.
        isource_ident = np.array([[-1, -1, -1],
                                  [-1, -1, +1],
                                  [-1, +1, -1],
                                  [-1, +1, +1],
                                  [+1, -1, -1],
                                  [+1, -1, +1],
                                  [+1, +1, -1],
                                  [+1, +1, +1]])
        # Includes/excludes bx, by, bz depending on 0/1 state.
        surface_coeff = np.array([[0, 0, 0],
                                  [0, 0, 1],
                                  [0, 1, 0],
                                  [0, 1, 1],
                                  [1, 0, 0],
                                  [1, 0, 1],
                                  [1, 1, 0],
                                  [1, 1, 1]])
        qq = surface_coeff[:,0] #  for bx1
        jj = surface_coeff[:,1] #  for by1
        kk = surface_coeff[:,2] #  for bz1
        F_abs_N = self.F_abs/nyquist # Normalise the standard absorption frequency range for surfaces, (0 to 1) = (0 to Fs/2)
        N_refl = int(2*np.round(nyquist/self.F_abs[1])) # Required length of FIR filter modelling impulse response of surface(+air)
        Half_I = int(N_refl/2) # Half length of FIR filter model
        Half_I_plusone = Half_I+1 # Half length shift required for FIR filter model of surface impulse response
        # Compute the (N_refl+1) point column vector Hann window
        window = 0.5*(1 - np.cos(2*np.pi*np.arange(0, N_refl+1).T/N_refl))
        #Image locations and impulse responses
        isource_xyz = np.zeros((3, n_isources)) # image source co-ordinates
        RR = len(self.F_abs);
        refl = np.zeros((RR, n_isources)) # surface reflection impulse amplitude (MEMORY CRITICAL POINT)
        xx = isource_ident[:,0]*source_xyz[0] # partial x coord of image.
        yy = isource_ident[:,1]*source_xyz[1] # partial y coord of image.
        zz = isource_ident[:,2]*source_xyz[2] # partial z coord of image.
        xx_yy_zz = np.array([xx, yy, zz])
        n_images=-1; #number of significant images of each parent source
        # Frequency dependent surface reflection and coordinates and distance for each image
        B = np.sqrt(1-self.A);
        bx1, bx2, by1, by2, bz1, bz2 = B

        for n in np.arange(-order_x, order_x+1, 1):
            bx2_abs_n = bx2**np.abs(n) # Compute here for efficiency
            Two_n_Lx = n*Two_Lx # Compute here for efficiency
            for l in np.arange(-order_y, order_y+1, 1):
                bx2y2_abs_nl = bx2_abs_n*(by2**np.abs(l)) # Compute here for efficiency
                Two_l_Ly=l*Two_Ly # Compute here for efficiency
                for m in np.arange(-order_z, order_z+1, 1):
                    # Compute here for efficiency
                    bx2y2z2_abs_nlm=bx2y2_abs_nl*(bz2**np.abs(m))
                    Two_m_Lz=m*Two_Lz # Compute here for efficiency
                    # Concatenate here for efficiency
                    Two_nlm_Lxyz = [Two_n_Lx, Two_l_Ly, Two_m_Lz]
                    for permu in np.arange(8):
                        n_images=n_images+1 #Accumulate count of the image sources
                        # calculate xyz coordinates of image source n_images
                        isource_xyz[:,n_images] = Two_nlm_Lxyz - xx_yy_zz[:,permu]
                        # minimal delay to sensors in samples
                        delay=np.min(Fs_c*np.sqrt(np.sum(\
                                    (isource_xyz[:, n_images] - \
                                     self.sensor_xyz)**2, axis=1)));
                        # compute only for image sources within impulse response length
                        if delay <= H_length:
                            refl[:,n_images] = bx1**np.abs(n-qq[permu])*\
                                                by1**np.abs(l-jj[permu])*\
                                                bz1**np.abs(m-kk[permu])*\
                                                bx2y2z2_abs_nlm
                            # (NB refl always +ve for air to surface, otherwise need abs here)
                            if np.sum(refl[:,n_images]) < 1E-6:
                                # Delete image sources with a sum of reflection coeffs below 1*10^-6 i.e. -120dB
                                n_images=n_images-1
                        else:
                            # Delete image sources with a delay > impulse response length H_length
                            n_images=n_images-1
        # Complete impulse response for the source
        n_images = n_images + 1
        isource_xyz = isource_xyz[:, :n_images] # Re-Allocate array for image source co-ordinates (discard trailing zero values)
        refl = refl[:, :n_images] # Re-Allocate array for surface reflection impulse amplitude (discard trailing zero values)
        H = np.zeros((H_length, self.channels))
        m_air = 6.875e-4*(self.F_abs/1000)**(1.7)
        # attenuation factors for one metre travelled in air
        temp_count = 0
        atten_air = np.exp(-0.5*m_air).T
        for mic in self.mics:
            # Get the sensor direction-dependent impulse responses
            sensor_dir = mic.direction
            sensor_dir = np.loadtxt(sensor_dir+'.txt')
            sensor_No = int(mic._id)-1
            # for each of the n_images image sources
            for idx_image in np.arange(n_images):
                b_refl = refl[:, idx_image]
                # Position vector from sensor_No to source(idx_image)
                xyz = isource_xyz[:, idx_image]-self.sensor_xyz[sensor_No, :]
                # Distance (m) between image source(idx_image) and sensor_No
                dist = np.sqrt(np.sum(xyz**2))
                # Include effect of distance (ie. 1/R) attenuation
                b_refl = b_refl/dist
                # Include the absorption due to air
                b_refl = b_refl*(atten_air**dist)
                # Estimate the values of reflection coefficient at the linear
                # interpolated grid points
                b_refl_func = interp1d(F_abs_N, b_refl)
                b_refl = b_refl_func(1.0/Half_I*np.arange(Half_I+1))
                # Half spectrum of data b_refl is now made conjugate-symmetric
                #about Nyquist frequency, and last data point
                #discarded to make periodic spectrum corresponding to a real data sequence.
                b_refl = np.hstack((b_refl, b_refl[::-1][1:-1]))
                # Transform surface data from frequency response to impulse response.
                # IFFT to calculate impulse response column vector of length N_refl samples
                h_refl = np.real(np.fft.ifft(b_refl, N_refl))
                # Make the impulse realisable (half length shift) and Hann window it
                h_refl = window*np.hstack((h_refl[Half_I_plusone-1:N_refl], h_refl[:Half_I_plusone]))
                # For primary sources, and image sources with impulse response peak magnitudes >= -100dB (1/100000)
                if (n_images==1) or  np.max(np.abs(h_refl[:Half_I_plusone])) >= 1E-5:
                    # Fractional delay filter
                    delay = Fs_c*dist; # delay in samples = (Samples per metre)*Distance
                    rdelay = np.round(delay); # Extract integer delay (concatenated later with impulse response)
                    t_Td = t-(delay-rdelay)*sampling_period; # Take account of fractional delay  -0.5 < D < +0.5 sample period
                    hsf=.5*(1+np.cos(Two_pi_Tw*t_Td))*np.sinc(self.sampling_rate*t_Td); # Compute delayed filter impulse response for sensor
                    # Convolve channel signals
                    sig_to_conv = np.vstack((h_refl.reshape(len(h_refl), 1), pad_frac))
                    sig_to_conv = sig_to_conv.reshape(len(sig_to_conv),)
                    h = scipy_sig.lfilter(hsf, 1, sig_to_conv)
                    len_h=len(h); # length of impulse response modelling image source response
                    adjust_delay = int(rdelay - np.ceil(len_h/2.0)) # Half length shift to remove delay due to impulse response
                    # Sensor filter
                    # position vector from each sensor location to each image source in sensor axes system
                    xyz_source = np.dot(self.tm_sensor[sensor_No, :, :], xyz)
                    # Distance (m) between sensor_No and proj of image source on xy plane
                    hyp = np.sqrt(xyz_source[0]**2+xyz_source[1]**2);
                    elevation = np.arctan(xyz_source[2]/(hyp+np.finfo(float).eps)); # Calculate -pi/2 <= elevation <= +pi/2 rads
                    azimuth = np.arctan2(xyz_source[1],xyz_source[0]); # Calculate -pi <= azimuth <= +pi rad
                    e_index = int(np.round(elevation*180/np.pi)+90)
                    a_index = int(np.round(azimuth*180/np.pi)+180)
                    sensor_ir=[sensor_dir[e_index,a_index]]
                    #h=scipy_sig.lfilter(sensor_ir,1,np.hstack((h, np.zeros((len(sensor_ir)-1,1)))))
                    h = scipy_sig.lfilter(sensor_ir, 1, h)
                    # Source filter
                    # position vector from each image source location to each sensor in source axes system
                    xyz_sensor = -1 * np.dot(tm_source, xyz)
                    # Distance (m) between image source and proj of sensor_No on xy plane
                    hyp = np.sqrt(xyz_sensor[0]**2 + xyz_sensor[1]**2)
                    # Calculate -pi/2 <= elevation <= +pi/2 rads
                    elevation=np.arctan(xyz_sensor[2]/(hyp+np.finfo(float).eps))
                    # Calculate -pi <= azimuth <= +pi rad
                    azimuth=np.arctan2(xyz_sensor[1],xyz_sensor[0])
                    e_index = int(np.round(elevation*180/np.pi)+90)
                    a_index = int(np.round(azimuth*180/np.pi)+180)
                    source_ir = [source_dir[e_index, a_index]]
                    #h = scipy_sig.lfilter(source_ir,1,np.hstack((h, np.zeros((len(source_ir)-1,1)))))
                    h = scipy_sig.lfilter(source_ir, 1, h)
                    len_h = len(h);
                    #Accumulate the impulse responses from each image source within an array of length H_length
                    start_index_Hp = max(adjust_delay+(adjust_delay >= 0), 0)
                    stop_index_Hp = min(adjust_delay+len_h, H_length)
                    start_index_h = max(-adjust_delay, 0)
                    stop_index_h = start_index_h + (stop_index_Hp - start_index_Hp)
                    #print(temp_count, start_index_Hp, stop_index_Hp, start_index_h, stop_index_h)
                    temp_count += 1
                    if stop_index_h < 0:
                        continue
                    #Add whole or part of impulse response
                    H[start_index_Hp:stop_index_Hp, sensor_No] = H[start_index_Hp:stop_index_Hp, sensor_No] + h[start_index_h:stop_index_h];
            #High-pass filtering
            H[:, sensor_No] = scipy_sig.lfilter(b_HP, a_HP, H[:, sensor_No])
        return H


    def __create_psi_theta_phi(self, source_off):
        c_psi = np.cos(np.pi/180*source_off[0])
        s_psi = np.sin(np.pi/180*source_off[0])
        c_theta = np.cos(-np.pi/180*source_off[1])
        s_theta = np.sin(-np.pi/180*source_off[1])
        c_phi = np.cos(np.pi/180*source_off[2])
        s_phi = np.sin(np.pi/180*source_off[2])
        return [c_psi, s_psi, c_theta, s_theta, c_phi, s_phi]

    def __create_tm(self, psi_theta_phi):
        c_psi, s_psi, c_theta, s_theta, c_phi, s_phi = psi_theta_phi
        tm_source = np.array([[c_theta*c_psi, \
                        c_theta*s_psi, \
                        -s_theta], \
                   [s_phi*s_theta*c_psi-c_phi*s_psi, \
                        s_phi*s_theta*s_psi+c_phi*c_psi, \
                        s_phi*c_theta], \
                   [c_phi*s_theta*c_psi+s_phi*s_psi, \
                        c_phi*s_theta*s_psi-s_phi*c_psi, \
                        c_phi*c_theta]])
        return tm_source


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='Config file')
    parser.add_argument('source_pos_x', help='Source x pos')
    parser.add_argument('source_pos_y', help='Source y pos')
    parser.add_argument('source_pos_z', help='Source z pos')
    parser.add_argument('out_file', help='File to write the RIR')
    args = parser.parse_args()
    source_pos = [float(args.source_pos_x), \
                    float(args.source_pos_y),\
                    float(args.source_pos_z)]
    sim_rir = RoomSim.init_from_config_file(args.config)
    rir = sim_rir.create_rir(source_pos)
    np.savetxt(args.out_file, rir)
