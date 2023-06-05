import h5py
import pandas as pd

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


