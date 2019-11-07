import numpy as np

class CoOrdinates(object):
    '''
    x y z value container
    '''
    def __init__(self, pos_x, pos_y, pos_z):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.pos_z = pos_z
        self.dim = [self.pos_x, self.pos_y, self.pos_z]

    def __str__(self):
        return 'x:'+str(self.pos_x)+\
                ' y:'+str(self.pos_y)+\
                ' z:'+str(self.pos_z)

    def dim_to_str(self):
        '''
        Convert dimension to string. Usefull for file name creation
        '''
        return '_'.join([str(ele) for ele in self.dim])


class RoomSpec(CoOrdinates):
    '''
    A position holder
    '''
    def __init__(self, pos_x, pos_y, pos_z):
        super(RoomSpec, self).__init__(pos_x, pos_y, pos_z)

    def verify_if_inside_room(self, pos, delta=0.0):
        '''
            Check if a position is inside a room.
            delta is the buffer
        '''
        if pos.pos_x < 0 or pos.pos_y < 0 or pos.pos_z <0:
            return False
        if pos.pos_x + delta < self.pos_x and \
                pos.pos_y + delta < self.pos_y and \
                pos.pos_z + delta < self.pos_z:
            return True
        return False

    def verify_if_inside_room_with_pos(self, pos_x, pos_y, pos_z, delta=0.0):
        '''
            Check if a position is inside a room.
            delta is the buffer
        '''
        if pos_x < 0 or pos_y < 0 or pos_z <0:
            return False
        if pos_x + delta <= self.pos_x and \
                pos_y + delta <= self.pos_y and \
                pos_z + delta <= self.pos_z:
            return True
        return False

    def verify_if_inside_room_multi_pos(self, pos_list, delta=0.0):
        '''
        Verify if points are inside room for multiple points
        '''
        status = []
        for pos in pos_list:
            status.append(self.verify_if_inside_room(pos, delta))
        return status


class Position(CoOrdinates):
    '''
    A position holder
    '''
    def __init__(self, pos_x, pos_y, pos_z):
        super(Position, self).__init__(pos_x, pos_y, pos_z)


    def compute_distance(self, pos):
        '''
        Compute distance to other position
        '''
        return np.sqrt(np.sum(
            (self.pos_x - pos.pos_x)**2+
            (self.pos_y - pos.pos_y)**2+
            (self.pos_z - pos.pos_z)**2))

    def compute_distance_with_val(self, pos_x, pos_y, pos_z):
        '''
        Compute distance to other position
        '''
        return np.sqrt(np.sum(
            (self.pos_x - pos_x)**2+
            (self.pos_y - pos_y)**2+
            (self.pos_z - pos_z)**2))

    def compute_multiple_distance(self, pos_list):
        '''
        Compute distances for multiple positions
        '''
        distance = []
        for pos in pos_list:
            distance.append(self.compute_distance(pos))
        return np.array(distance)

def new_element_pos(room_dim, other_mics=None, dist_from_wall=0.5):
    '''
        Create a new microphone position. Make sure it is far of from the
        other microphone positions
    '''
    if other_mics is not None and len(other_mics) == 0:
        other_mics = None
    pos_min = dist_from_wall
    x_pos_max = room_dim.pos_x - dist_from_wall
    y_pos_max = room_dim.pos_y - dist_from_wall
    z_pos_max = room_dim.pos_z - dist_from_wall
    x_space = x_pos_max-pos_min
    y_space = y_pos_max-pos_min
    z_space = z_pos_max-pos_min

    assert x_space > 0 and y_space > 0 and z_space > 0,\
            'No space in the room to keep a element. \
            Change ROOM settings in the config file'

    pos_x = x_space * np.random.rand() + pos_min
    pos_y = y_space * np.random.rand() + pos_min
    pos_z = z_space * np.random.rand() + pos_min
    pos = Position(pos_x, pos_y, pos_z)
    if other_mics is not None:
        if np.sum(pos.compute_multiple_distance(other_mics) < \
                0.5) > 0:
            pos = new_element_pos(room_dim, other_mics)
    return pos


def create_new_room(room_max, room_min):
    '''
    Decide soom width
    '''
    room_range = room_max - room_min 
    x_val = np.random.rand() * room_range + room_min 
    y_val = np.random.rand() * room_range + room_min 
    z_val = np.random.rand() * room_range + room_min 
    room_obj = RoomSpec(x_val, y_val, z_val)
    return room_obj


def generate_rt60(low_rt, high_rt):
    '''
        Get a random RT60
    '''
    rt60 = (high_rt - low_rt) * np.random.rand() + low_rt
    return rt60
