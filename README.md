rir_simulator_python
--------------------
Room impulse response simulator using python

Copyright 2003 Douglas R. Campbell

Copyright 2008-2016 Emmanuel Vincent

Copyright 2017 Sunit Sivasankaran

 This software is a python version of the stripped-down version of the Roomsim toolbox version 3.3 by Douglas R. Campbell.
 The matlab function for the stripped down version as well as RIR generation for moving sources can be found here:
Roomsimove, http://homepages.loria.fr/evincent/software/Roomsimove.zip
 
If you find the code useful, please cite the following reference:
Room Impulse Response Generator, https://github.com/sunits/rir_simulator_python


One  difference between the matlab version and this code is that 
RT60 value is assumed to be same for all frequencies.

Tested for sampling rate of 16000 Hz. 

Usage:
-------

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
    room_dim = [4.2, 3.4, 5.2]
    room = roomsimove_single.Room(room_dim)
    mic_pos = [2, 2, 2]
    mic1 = roomsimove_single.Microphone(mic_pos, 1,  \
            orientation=[0.0, 0.0, 0.0], direction='omnidirectional'):
    mic_pos = [2, 2, 1]
    mic2 = roomsimove_single.Microphone(mic_pos, 2,  \
            orientation=[0.0, 0.0, 0.0], direction='cardioid'):
    mics = [mic1, mic2]
    sampling_rate = 16000
    sim_rir = roomsimove_single.RoomSim(sample_rate, room, mics, RT60=300)
    source_pos = [1, 1, 1]
    rir = sim_rir.create_rir(source_pos)

Appyling RIR to data
-------------------
    import fftfilt
    import sounfile as sf
    # Assuming single channel data
    [data, fs] = sf.read(wav_file)
    reverb_data = fftfilt.fftfilt(rir,data)



    
