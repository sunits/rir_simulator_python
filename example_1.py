import roomsimove_single
import ipdb
import soundfile as sf
import olafilt
import numpy as np


rt60 = 0.4 # in seconds
room_dim = [4.2, 3.4, 5.2] # in meters
mic_pos1 = [2, 2, 2] # in  meters
mic_pos2 = [2, 2, 1] # in  meters
source_pos = [1, 1, 1] # in  meters
sampling_rate = 16000

mic_positions = [mic_pos1, mic_pos2]
rir = roomsimove_single.do_everything(room_dim, mic_positions, source_pos, rt60)

## Alternative method for more control
# absorption = roomsimove_single.rt60_to_absorption(room_dim, rt60)
# room = roomsimove_single.Room(room_dim, abs_coeff=absorption)
# mic1 = roomsimove_single.Microphone(mic_pos1, 1,  \
#         orientation=[0.0, 0.0, 0.0], direction='omnidirectional')
# mic2 = roomsimove_single.Microphone(mic_pos2, 2,  \
#         orientation=[0.0, 0.0, 0.0], direction='cardioid')
# mics = [mic1, mic2]
# sim_rir = roomsimove_single.RoomSim(sampling_rate, room, mics, RT60=rt60)
# rir = sim_rir.create_rir(source_pos)
ipdb.set_trace()

[data, fs] = sf.read('data.wav',always_2d=True)
data =  data[:,0]
data_rev_ch1 = olafilt.olafilt(rir[:,0], data)
data_rev_ch2 = olafilt.olafilt(rir[:,1], data)
data_rev = np.array([data_rev_ch1, data_rev_ch2])
sf.write('data_rev.wav', data_rev.T, fs)
