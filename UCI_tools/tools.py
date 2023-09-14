import csv
import h5py
import numpy as np
import pandas as pd

try:
    '''
    The user should have a paths.py script defining the idiosyncratic
    locations of folders which the methods below need to use.
    E.g. The contents of the user's paths.py might be as follows:
        data = '$HOME/project/data/'
        paper = '$HOME/project/paper' 
    This is necessary because each user and each project will have a
    different file structure
    '''
    import paths
except:
    pass

def save_var_latex(key, value, fname='data.txt'):
    '''
    Save a data point do a file in paths.paper so the user's LaTeX paper can
    use it.

    Parameters
    ----------
    key: str
        The key that the LaTeX paper will use to look up this data point.
    value: object
        The value of the data point the user wants to save for the LaTeX paper
        to use. The code is written so any object (e.g. float, int, str) should
        work, but it's probably best to use a str.
    fname: str, default 'data.txt'
        The file in paths.paper to which the code should save the prediction 
        and its uncertainty.

    Returns
    -------
    None
    '''

    dict_var = {}

    file_path = paths.paper + fname 

    try:
        with open(file_path, newline="") as file:
            reader = csv.reader(file)
            for row in reader:
                dict_var[row[0]] = row[1]
    except FileNotFoundError:
        pass

    # The code I took this from was written to handle saving multiple datum at
    # once, so that's why it's written this way even though we only have one
    # key, value pair. I'm keeping the code as in in case we want to expand its
    # functionality in the future.
    dict_var[key] = value

    with open(file_path, "w") as f:
        for key in dict_var.keys():
            f.write(f"{key},{dict_var[key]}\n")

    return None

def save_prediction(string, y, dy, fname='data.txt'):
    '''
    Save a prediction and its uncertainty to `fname` to be used by a Latex 
    file.

    Parameters
    ----------
    string: str 
	The identifier for the variable

    y: str 
        The value of the target variable

    dy: np.ndarray of str
	The uncertainty in the variable. The user *must* provide it as an
        np.ndarray, even if it's only one single symmetric uncertainty.

        E.g. If the error is symmetric, the user would provide something like
             `dy=np.array(['0.1'])`.
             If asymmetric, `dy=np.array(['0.1', '0.2'])`.

    fname: str, default 'data.txt'
        The file in paths.paper to which the code should save the prediction 
        and its uncertainty.

    Returns
    -------
    None
    ''' 

    save_var_latex(string, y)
    if len(dy)==1:
        save_var_latex('d'+string, dy[0])
        save_var_latex('d'+string+'_plus', 'n/a')
        save_var_latex('d'+string+'_minus', 'n/a')
    elif len(dy)==2:
        save_var_latex('d'+string, 'n/a')
        save_var_latex('d'+string+'_plus', dy[0])
        save_var_latex('d'+string+'_minus', dy[1])
    else:
        raise ValueError('Margin of error has a maximum of two elements.')
    return None

def read_snapshot_simple( filepath, particle_type='PartType0' ):
    f = h5py.File( filepath, 'r' )

    data = {}
    for column in f[particle_type].keys():
        data_in_col = f[particle_type][column][...]
        if len( data_in_col.shape ) > 1:
            for i in range( data_in_col.shape[1] ):
                data[column + str(i)] = data_in_col[:,i]
        else:
            data[column] = data_in_col
            
    p = pd.DataFrame( data )
            
    return p

def fe_over_h_ratios(mfrac,he_frac,fe_frac):
    h_frac=1-mfrac-he_frac
    #...some constants                                                                                           
    sun_fe_h_frac = 0.0030/91.2
    mass_h = 1.0084 # in Atomic Mass Units                                                                                              
    mass_fe= 55.845 # in Atomic Mass Units
    #...Need to convert mass fractions to number fractions                                                                                                
    fe_h_num = (fe_frac/h_frac)*(mass_h/mass_fe)
    #...Abundance ratio                                                                                                                            
    ab_fe_h = np.asarray(np.log10(fe_h_num/sun_fe_h_frac))
    return ab_fe_h

def sft_to_ages(sft):
    from astropy.cosmology import Planck13 
    z = (1/sft)-1
    ages = np.array((Planck13.lookback_time(z)))
    return ages

def calc_cyl_vels(v_vecs_rot, coords_rot):
    '''
    Put velocities into cylindrical coordinates given Cartesian velocites and 
    coordinates. These velocity and coordinate arguments should be centered and
    rotated so their z components align with the total angular momentum of the
    galaxy's stars (i.e. so their z component is perpendicular to the disc).
    
    Parameters
    ----------
    v_vecs_rot: np.ndarray, shape=(number of particles, 3)
        Velocity vectors in Cartesian coordinates rotated so their z-axis 
        aligns with the angular
        momentum vector of the stars.
    coords_rot: np.ndarray, shpae=(number of particles, 3)
        Position vectors in Cartesian coordinates rotated so their z-axis
        aligns with the angular momentum vector of the stars.
    
    Returns
    -------
    d: dict
        Dictionary containing the velocity information for the particles.
        It contains the folowing key, value pairs.
        'v_vec_disc': Only x and y components of velocity
        'coord_disc': Only x and y components of coordinates
        'v_dot_rhat': The r (scalar) component of velocity, where 
            d['coord_disc'] gives
            the r vector (can be negative if the projection of velocity along r
            points in the opposite direction of r)
        'v_r_vec': The projection of velocity along the r vector, expressed in
            Cartesian x and y components
        'v_phi_vec': The projection of velocity along the phi vector, expressed
            in Cartesian x and y components
        'v_phi_mag': The magnitude of the projection of velocity along the phi
            vector (always positive)
        'v_dot_phihat': The phi (scalar) component of velocity (can be negative
            if the projection of velocity along phi points in the opposite
            direction of phi)
    '''

    # Initialize a dictionary into which we put particle kenematic properties
    d = {} 
    #xy component of velocity
    d['v_vec_disc'] = v_vecs_rot[:,:2] 
    #xy component of coordinates
    d['coord_disc'] = coords_rot[:,:2] 

    ###################################################################
    ## Find the projection of velocity onto the xy vector (i.e. v_r) 
    ## v_r = v dot r / r^2 * r 
    ###################################################################
    vdotrs = np.sum(d['v_vec_disc'] \
                   * d['coord_disc'], axis=1)
    rmags = np.linalg.norm(d['coord_disc'], axis=1)
    d['v_dot_rhat'] = vdotrs / rmags

    #v dot r / r^2
    vdotrs_r2 = (vdotrs \
        / np.linalg.norm(d['coord_disc'], axis=1) **2.)
    #Need to reshape v dot r / r^2 so its shape=(number of particles,1)
    #so we can mutilpy those scalars by the r vector
    vdotrs_r2 = vdotrs_r2.reshape(len(d['coord_disc']), 1)
    d['v_r_vec'] = vdotrs_r2 * d['coord_disc']
    ###################################################################

    #v_phi is the difference of the xy velocity and v_r
    d['v_phi_vec'] = d['v_vec_disc'] - d['v_r_vec']
    d['v_phi_mag'] = np.linalg.norm(d['v_phi_vec'],axis=1)

    ###################################################################
    ## Finding the vphi in \vec{vphi} = vphi * \hat{vphi} 
    ## (i.e. v dot phi )
    ## Note v_phi_mag = |v dot phi|, and we want to determine whether
    ## v dot phi is positive or negative.
    ###################################################################
    rxvs = np.cross(d['coord_disc'], 
                    d['v_vec_disc']) #r cross v
    # v dot phi is positive if the z component of r cross v is 
    # positive,
    # because this means its angular momentum is in the same direction
    # as the disc's.
    signs = rxvs/np.abs(rxvs) 
    d['v_dot_phihat'] = d['v_phi_mag']*signs
    ###################################################################

    return d
