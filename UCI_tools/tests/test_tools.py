from mock import patch
import numpy as np
import numpy.testing as npt
import unittest

import UCI_tools.tools as tools

########################################################################

class TestLoadFIREData( unittest.TestCase ):
    '''Test suite for loading FIRE data.
    '''

    ########################################################################

    def test_simple( self ):
        '''This test just ensures that the data can be loaded at all.'''

        sim_dir = './tests/test_data/downsampled_sim_data/fire_sim/output'
        snapshot = 600
        particle_type = 'PartType0'

        snapshot = tools.read_snapshot_simple(
            sim_dir = sim_dir,
            snapshot = snapshot,
            particle_type = particle_type,
        )

        assert 'Density' in snapshot.columns
        assert 'Coordinates0' in snapshot.columns

    ########################################################################

    def test_error( self ):
        '''This function tests if the code sends the proper error when it breaks.'''

        assert False