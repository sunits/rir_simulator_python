import olafilt
import roomsimove_single as rs
import numpy as np
import utils
import ipdb


class RandomRIR(object):
    """
    Generate a random room, microphone and source  position and generate the corresponding RIR. 
    
    # Arguments
        sampling_rate: Sampling rate of the RIR
        max_rt_60: Maximum value of RT60 in seconds. Actual RT60 is random between [0.1, max_rt_60]
        min_room_di, max_room_dim: Min and Maximum value of the room dim. 
                Room dimensions are random picks between [min_room_dim, max_room_dim]

    # Usage
    rir_if = RandomRIR(sampling_rate=16000)
    src = [np.random.rand(10000), np.random.rand(10000)]
    rev_sig = rir_if.reverberate(src)
    
    """
    def __init__(self, sampling_rate, max_rt_60=0.5, min_room_dim=3, max_room_dim=5):
        self.sampling_rate = sampling_rate
        self.max_rt_60 = max_rt_60
        self.max_room_dim = max_room_dim
        self.min_room_dim = min_room_dim
    
    def create_rir(self, src_cnt, mic_cnt=1):
        room_dim = utils.create_new_room(self.min_room_dim, self.max_room_dim)
        room = rs.Room(room_dim.dim)
        rt60 = utils.generate_rt60(0.1, self.max_rt_60)
        all_ele = []
        all_mics = []
        for mic_id in np.arange(mic_cnt):
            mic_pos = utils.new_element_pos(room_dim, all_ele)
            mic = rs.Microphone(mic_pos.dim, 2,  \
                    orientation=[0.0, 0.0, 0.0], direction='cardioid')
            all_mics.append(mic)
            all_ele.append(mic_pos)
        all_srcs = []
        for mic_id in np.arange(src_cnt):
            src_pos = utils.new_element_pos(room_dim, all_ele)
            all_srcs.append(src_pos)
            all_ele.append(src_pos)
        all_rir = []
        sim_rir = rs.RoomSim(self.sampling_rate, room, all_mics, RT60=rt60)
        for src in all_srcs:
            rir = sim_rir.create_rir(src.dim)
            all_rir.append(rir)
        return all_rir
    
    def reverberate(self, src_list, mic_cnt=1):
        """
        Create the RIR with random values and convolves with sources
        # Arguments:
            src_list: wav signals for different sources
            mic_cnt: Number of micrphones
        """
        src_cnt = len(src_list)
        rirs = self.create_rir(src_cnt, mic_cnt=mic_cnt)
        rev_sig = []
        for src_idx, src_rir in enumerate(rirs):
            src_ch = [] # multiple channels
            for mic_src_rir in src_rir.T:
                data_rev = olafilt.olafilt(mic_src_rir, src_list[src_idx])
                src_ch.append(data_rev)
            src_ch = np.stack(src_ch, 1)
            rev_sig.append(src_ch)
        return rev_sig


if __name__ == '__main__':
    rir_if = RandomRIR(sampling_rate=16000)
    src = [np.random.rand(10000), np.random.rand(10000)]
    ipdb.set_trace()
    rev_sig = rir_if.reverberate(src)
