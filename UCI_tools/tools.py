import csv
import h5py
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
