#from mock import patch
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
        '''
        This function tests if the code sends the proper error when it breaks.
        '''

        desired_msg = (
            'Cannot find snapshot at specified locations. '
            + 'Locations searched:\n'
            + 'this_dir_does_not_exist/snapdir_600\n'
            + 'this_dir_does_not_exist/snapshot_600.hdf5'
        )
        with self.assertRaises( IOError ) as error:
            tools.read_snapshot_simple(
                sim_dir = 'this_dir_does_not_exist',
                snapshot = 600,
                particle_type = 'PartType4',
            )

        self.assertEqual( str( error.exception ), desired_msg )

#################################################################

class TestMisc( unittest.TestCase  ):
    ''' Testing nonloaded data
    '''
	
    ########################################################################

    def test_sft_to_ages( self ):
        npt.assert_allclose(tools.sft_to_ages(1), 0, atol = .15) #snapshot 600
        npt.assert_allclose(
                tools.sft_to_ages(0.8550955),
                13.79874688-11.69441659,
                atol = .15
            ) #snapshot 500
        npt.assert_allclose(
                tools.sft_to_ages(0.6958599),
                13.79874688-9.19969494,
                atol = .15
            ) #snapshot 400
        npt.assert_allclose(
                ools.sft_to_ages(0.5366242),
                13.79874688-6.58906279,
                atol = .15
            ) #snapshot 300
        npt.assert_allclose(
                tools.sft_to_ages(0.3777778),
                13.79874688-4.04069309,
                atol = .15
            ) #snapshot 200
        npt.assert_allclose(
                tools.sft_to_ages(0.2187500),
                13.79874688-1.81321181,
                atol = .15
            ) #snapshot 100
        npt.assert_allclose(
                tools.sft_to_ages(0.0100000),
                13.79874688-0.01780470,
                atol = .15
            ) #snapshot 0	
