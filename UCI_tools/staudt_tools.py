def init_df(mass_class=12):
    '''
    Initialize a dataframe specifying information about the file names and
    resolutions of simulated galaxies of
    a given mass class.

    Parameters
    ----------
    mass_class: {9, 10, 11, 12, 'all'}
        The log10(mass / M_sun) mass class of galaxies in which the user is
        interested. The user can also specify 'all' if they want the dataframe
        to contain information for all masses classes from 9 through 12.

    Returns
    -------
    df: pd.DataFrame
    '''
    import pandas as pd
    import numpy as np

    if not isinstance(mass_class,(int,str)):
        raise ValueError('mass_class should be an integer or \'all\'.')
    suffixes_m9 = ['']
    # It's necessary to have different suffix lists for cropped and not cropped
    # because for cropped, I separate the RJ, RR, RL pair folders.
    suffixes_m9cropped = suffixes_m9.copy() 
    ress_m9 = [250]
    gal_names_m9 = ['m9']
    host_keys_m9 = ['host.index']
    mass_classs_m9 = [9]

    '''
    suffixes_m10 = ['b','c','d','e','f','g','h','i','j','k','l','m','q','v','w',
                'xa','xb','xc','xd','xe','xf','xg','xh','xi',
                'y','z']
    ress_m10 = [500]*12 + [250]*2 + [2100] + [4000]*9 + [2100]*2
    '''
    suffixes_m10 = ['q','v','w','y','z']
    suffixes_m10cropped = suffixes_m10.copy()
    #ress_m10 = [500]*12 + [250]*2 + [2100] + [4000]*9 + [2100]*2
    ress_m10 = [250]*2 + [2100]*3
    gal_names_m10=['m10' + g for g in suffixes_m10]
    host_keys_m10 = ['host.index']*len(suffixes_m10)
    mass_classs_m10 = [10]*len(suffixes_m10)
    assert len(suffixes_m10) == len(ress_m10) == len(host_keys_m10) \
           == len(gal_names_m10) == len(mass_classs_m10)

    #suffixes_m11 = ['a','b','c','d','e','f','g','h','i','q','v']
    #ress_m11 = [2100]*3 + [7100]*2 + [17000]*2 + [7100]*4
    #f and g don't have halo files, so I'm skipping them for now.
    suffixes_m11 = ['a','b','c','d','e','h','i','q','v']
    suffixes_m11cropped = suffixes_m11.copy()
    ress_m11 = [2100]*3 + [7100]*2 + [7100]*4
    gal_names_m11 = ['m11' + g for g in suffixes_m11]
    host_keys_m11 = ['host.index']*len(suffixes_m11)
    mass_classs_m11 = [11]*len(suffixes_m11)
    assert len(suffixes_m11) == len(ress_m11) == len(host_keys_m11) \
           == len(gal_names_m11) == len(mass_classs_m11)

    suffixes_m12 = ['b','c','f','i','m','r','w','z',
                '_elvis_RomeoJuliet',
                '_elvis_RomeoJuliet',
                '_elvis_RomulusRemus',
                '_elvis_RomulusRemus',
                '_elvis_ThelmaLouise',
                '_elvis_ThelmaLouise']
    suffixes_m12cropped = ['b','c','f','i','m','r','w','z',
                '_elvis_Romeo',
                '_elvis_Juliet',
                '_elvis_Romulus',
                '_elvis_Remus',
                '_elvis_Thelma',
                '_elvis_Louise']
    ress_m12 = [7100,7100,7100,7100,7100,7100,7100,4200,
                3500,3500,4000,4000,4000,4000]
    gal_names_m12=['m12' + g for g in suffixes_m12[:8]] + ['Romeo','Juliet',
                                                       'Romulus','Remus',
                                                       'Thelma','Louise']
    host_keys_m12 = ['host.index']*8+['host.index','host2.index']*3
    mass_classs_m12 = [12]*len(suffixes_m12)
    assert len(suffixes_m12) == len(ress_m12) == len(host_keys_m12) \
           == len(gal_names_m12) == len(mass_classs_m12)

    if mass_class in [9,10,11,12]:
        lcls=locals()
        exec('suffixes = suffixes_m'+str(mass_class),globals(),lcls)
        exec('suffixes_cropped = suffixes_m'+str(mass_class)+'cropped',
             globals(),lcls)
        exec('ress = ress_m'+str(mass_class),globals(),lcls)
        exec('gal_names = gal_names_m'+str(mass_class),globals(),lcls)
        exec('host_keys = host_keys_m'+str(mass_class),globals(),lcls)
        exec('mass_classs = mass_classs_m'+str(mass_class),globals(),lcls)
        suffixes=lcls['suffixes']
        suffixes_cropped = lcls['suffixes_cropped']
        ress=lcls['ress']
        gal_names=lcls['gal_names']
        host_keys=lcls['host_keys']
        mass_classs=lcls['mass_classs']
    elif mass_class=='all':
        '''
        suffixes = suffixes_m9 + suffixes_m10 + suffixes_m11 + suffixes_m12
        suffixes_cropped = suffixes_m9cropped + suffixes_m10cropped \
                           + suffixes_m11cropped + suffixes_m12cropped
        ress = ress_m9 + ress_m10 + ress_m11 + ress_m12
        gal_names = gal_names_m9 + gal_names_m10 + gal_names_m11 + gal_names_m12
        host_keys = host_keys_m9 + host_keys_m10 + host_keys_m11 + host_keys_m12
        '''
        #Removing m9 and m10 for now
        suffixes = suffixes_m11 + suffixes_m12
        suffixes_cropped = suffixes_m11cropped + suffixes_m12cropped
        ress = ress_m11 + ress_m12
        gal_names = gal_names_m11 + gal_names_m12
        host_keys = host_keys_m11 + host_keys_m12
        mass_classs = mass_classs_m11 + mass_classs_m12
    else:
        raise ValueError('mass_class must be 9,10,11,12 or \'all\'')
        
    df=pd.DataFrame(index=gal_names,
                    data=np.transpose([suffixes, suffixes_cropped, 
                                       ress, host_keys]),
                    columns=['fsuffix','fsuffix_cropped','res','host_key'])
    #Need to add mass_classs *after* creating df because we need it to have an 
    #int dtype.
    df['mass_class']=mass_classs 

    return df
