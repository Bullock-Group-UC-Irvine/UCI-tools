import glob
import h5py
import os
import pandas as pd

################################################################################

def read_snapshot_simple( sim_dir, snapshot, particle_type ):
    '''Minimal possible code for opening up FIRE data.
    Args:
        sim_dir (str):
            Directory containing all the simulation output,
            e.g. 'main_storage_dir/m12i_r7100/output'

        snapshot (int):
            The snapshot number you want to open, e.g. 600

        particle_type (str):
            The particle type you want to load, e.g.
            'PartType4' usually loads stars.
    '''

    # Check if the directory for this particular snapshot exists
    # os.path.join is the safer version of file_dir + '/' + file_name
    snapshot_dir = os.path.join( sim_dir, 'snapdir_{:03d}'.format( snapshot ) )
    snapshot_file_if_no_dir = os.path.join( sim_dir, 'snapshot_{:03d}.hdf5'.format( snapshot ) )
    if os.path.isdir( snapshot_dir ):
        filepaths = glob.glob( os.path.join( snapshot_dir, '*.hdf5' ) )
    # Sometimes the data is stored in single hdf5 files.
    elif os.path.isfile( snapshot_file_if_no_dir ):
        filepaths = [ snapshot_file_if_no_dir, ]
    # Raise an error if not found.
    else:
        raise IOError(
            'Cannot find snapshot at specified locations. ' + \
            'Locations searched:\n{}\n{}'.format(
                snapshot_dir,
                snapshot_file_if_no_dir
            )
        )

    # The snapshot may be broken up into multiple pieces,
    # so we loop through them.
    snapshot_pieces = []
    for filepath in filepaths:
    
        # Load the file itself
        f = h5py.File( filepath, 'r' )

        # Loop through each property, e.g. temperature, density, etc.
        data = {}
        for column in f[particle_type].keys():

            # Actually get the data
            data_in_col = f[particle_type][column][...]

            # Check the shape of the data.
            # If it is a multidimensional array,
            # such as the coordinates, then separate into multiple
            # 1D arrays
            if len( data_in_col.shape ) > 1:
                for i in range( data_in_col.shape[1] ):
                    data[column + str(i)] = data_in_col[:,i]
            # Regular case
            else:
                data[column] = data_in_col
                
        # Turn each piece into a dataframe
        p_i = pd.DataFrame( data )
        snapshot_pieces.append( p_i )

    # And then combine the dataframes!
    p = pd.concat( snapshot_pieces )
                
    return p
