# Python RIR Simulator
Room impulse response simulator using python

Copyright 2003 Douglas R. Campbell

Copyright 2008-2016 Emmanuel Vincent

Copyright 2017 Sunit Sivasankaran

 This software is a python version of the stripped-down version of the Roomsim toolbox version 3.3 by Douglas R. Campbell.
 The matlab function for the stripped down version as well as RIR generation for moving sources can be found here:
Roomsimove, https://www.irisa.fr/metiss/members/evincent/Roomsimove.zip
 
If you find the code useful, please cite the following reference:
Room Impulse Response Generator, https://github.com/sunits/rir_simulator_python


One  difference between the matlab version and this code is that 
RT60 value is assumed to be same for all frequencies.

Tested for sampling rate of 16000 Hz. 

## Usage:

### Generate random Room Impulse Responses
    import roomsimove_single
    # Create an interface. RIRs are not yet generated
    rir_if = roomsimove_single.RandomRIR(sampling_rate=16000)
    # Supports multiple sources. Each source is placed in a random position inside the same room.
    src = [np.random.rand(10000), np.random.rand(10000)] 
    # Support multiple microphones. Position of each device is also random.
    # Minimum distance between each element (mics or sources) is 0.5 meters
    rev_sig = rir_if.reverberate(src, mic_cnt=2)
    # rev_sig contains a list of reverberate sources. 
    # Each element in the list is of dimension [src_len x mic_cnt]


### As standalone file:
    python roomsimove_single.py config_file source_pos_x source_pos_y source_pos_z output_file

    The help options will also give the details
    python roomsimove_single.py -h
    
### As a module:
    using config_file
    -----------------
    import roomsimove_single
    sim_rir = roomsimove_single.RoomSim.init_from_config_file(config_file)
    source_pos = [1, 1, 1]
    rir = sim_rir.create_rir(source_pos)

    using default values of absorption coeffecients
    -----------------------------------------------
    import roomsimove_single
    room_dim = [4.2, 3.4, 5.2]
    room = roomsimove_single.Room(room_dim)
    mic_pos = [2, 2, 2]
    mic1 = roomsimove_single.Microphone(mic_pos, 1,  \
            orientation=[0.0, 0.0, 0.0], direction='omnidirectional'):
    mic_pos = [2, 2, 1]
    mic2 = roomsimove_single.Microphone(mic_pos, 2,  \
            orientation=[0.0, 0.0, 0.0], direction='cardioid'):
    mics = [mic1, mic2]
    sample_rate = 16000
    sim_rir = roomsimove_single.RoomSim(sample_rate, room, mics, RT60=0.3)
    source_pos = [1, 1, 1]
    rir = sim_rir.create_rir(source_pos)

    # Applying RIR to data
    import fftfilt
    import soundfile as sf
    # Assuming single channel data
    [data, fs] = sf.read(wav_file)
    reverb_data = fftfilt.fftfilt(rir,data)




    
