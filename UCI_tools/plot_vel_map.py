def get_m12_path_olti(sim_name, host_idx, snap):
    '''
    Generate the path to simulation files in Olti's directory.

    Parameters
    ----------
    sim_name: {
        'm12b_res7100',
        'm12c_res7100',
        'm12_elvis_RomeoJuliet_res3500',
        'm12_elvis_RomulusRemus_res4000',
        'm12_elvis_ThelmaLouise_res4000',
        'm12f_res7100',
        'm12i_res7100',
        'm12m_res7100',
        'm12r_res7100'
    }
        Name of the simulation to load.
    host_idx: int
        The index of the host to analyze in the given simulation. For Latte
        runs, this should always be 0. For Elvis pairs, the index will either
        be 0 or 1 depending on which of the the pair the user wants to analyze.
    snap: str
        The snapshot number to load. The snapshot number should be in
        string format with three digits.

    Returns
    -------
    path: str
        The path to the file being analyzed.
    '''
    import os
    path = os.path.join(
        '/DFS-L/DATA/cosmo/grenache/omyrtaj/analysis_data/metaldiff/',
        sim_name,
        'id_jnet_jzjc_jjc_'
            + snap
            + '_host'
            + str(host_idx)
            + '_20kpc_rockstar_centers_metaldiff.hdf5'
    )
    return path

def load_m12_data_olti(sim_path, snap):
    '''
    Load Olti's data for use with `UCI_tools.plot_vel_map.plot`

    Parameters
    ----------
    sim_path: str 
        The path to the simulation the user wants to analyze. The user could
        use UCI_tools.plot_vel_map.get_m12_path_olti to easily generate a path
        that leads to Olti's files, or they could supply their own path.
        Another option would be for the user to write their own `get_m12_path`
        method in UCI_tools.plot_vel_map and make a pull request so everyone
        has it.
    snap: str
        The snapshot number corresponding to `sim_path`. It should be in string
        format with three digits. The code uses this to
        display the correct look-back time in the plot.

    Returns
    -------
    data_out: dict
        A dictionary containing the data that `UCI_tools.plot_vel_map.plot`
        requires. Keys are as follows:
            'pos_gas': The centered, unrotated position vectors of the
                simulation's gas particles in physical kpc
            'vel_gas': The unrotated velocity vectors in Cartesian coordinates
                of the gas particles, 
                relative to
                the host center
            'jnet_gas': The net specific angular momentum vector of all the 
                gas within 20 kpc of the host center
            'temp': The temperature of the gas particles in Kelvin
            'mass_gas': The mass of each gas particle in physical units of 
                M_sun
            'pos_star': The centered, unrotated position vector of each star
                particle in the simulation in physical kpc
            'vel_star': The unrotated velocity vector in Cartesian coordinates
                of each star particle,
                relative to the host center
            'sft': The time in Gyr since the formation of each star particle,
                relative to the given snapshot.
            'jnet_star': The net specific angular momentum vector of all the
                stars within 20 kpc of the host center
    '''
    import h5py
    import numpy as np
    from . import paths

    snapshot_times = np.loadtxt(
        paths.snap_times
    )
    a = float(snapshot_times[int(snap)][1])

    data_out = {}
    with h5py.File(sim_path, 'r') as data:
        host_center = np.array(data['host_center'])
        host_vel = np.array(data['host_velocity'])

        # Load gas data
        data_out['pos_gas'] = a * (
            np.array(data['gas_coord_unrotated']) - host_center
        )
        data_out['vel_gas'] = np.array(data['gas_vel_unrotated']) - host_vel
        data_out['jnet_gas'] = np.array(data['jnet_gas'])
        data_out['temp'] = np.array(data['gas_temp'])
        data_out['mass_gas'] = np.array(data['mass_gas'])
        
        # Load star data
        data_out['pos_star'] = a * (
            np.array(data['star_coord_unrotated']) - host_center
        )
        data_out['vel_star'] = np.array(data['star_vel_unrotated']) - host_vel
        data_out['sft'] = np.array(data['sft_Gyr'])
        data_out['jnet_star'] = np.array(data['jnet_young_star'])
    return data_out

def plot(
        sim_path, 
        display_name,
        snap,
        gas_num=50,
        star_num=20,
        save_plot=True):
    '''
    Plot the v_y velocity map of gas and young stars for a given simulation and
    return the data for those two maps. 
    
    The grid orientation of the returned maps 
    follows the standard matrix convention specified in the 
    `matplotlib.axes.Axes.pcolormesh` documentation; They have shape 
    (nrows, ncolumns) with the column number as X and the row number as Y. 
    
    Note
    that this orientation is the transpose of the `H` output of 
    `np.histogram2d`,
    which has X data along the 0 axis and Y data 
    along the 1 axis.
    If one were to print out `H`, 
    visually,
    one might expect the X-axis to run horizontally along `H` and the 
    Y-axis to
    run vertically along `H`. However, the opposite is true. Additionally, 
    `matplotlib.axes.Axes.pcolormesh`
    plots the inputted `C` mesh arry with the column number as X and the row 
    number as Y.
    Therefore, we must provide the transpose of `H` to
    `pcolormesh`.

    Parameters
    ----------
    sim_path: str 
        The path to the simulation the user wants to analyze. The user could
        use UCI_tools.plot_vel_map.get_m12_path_olti to easily generate a path
        that leads to Olti's files, or they could supply their own path.
        Another option would be for the user to write their own `get_m12_path`
        method in UCI_tools.plot_vel_map and make a pull request so everyone
        has it.
    display_name: str 
        Simulation name to show in the plot.
    snap: str
        The snapshot number corresponding to `sim_path`. It should be in string
        format with three digits. The code uses this to
        display the correct look-back time in the plot.
    gas_num: int, default 50
        The minimum number of gas particles a 2d histogram bin must have
        for it to be included.
    star_num: int, default 20
        The minimum number of star particles a 2d histogram bin must have for
        it to be included.
    save_plot: bool, default True
        Whether to save the plot to disk. If True, the code will save the plot
        in the `figures` directory specified in the user's 
        config.ini file in
        their home directory. 

    Returns
    -------
    v_y_colormesh_gas: np.ndarray
        The gas velocity colormap data for use with 
        `matplotlib.axes.Axis.pcolormesh`. X data is along the 1 axis. Z data
        is along the 0 axis. The user can directly input this into
        `pcolormesh`.
    v_y_colormesh_star: np.ndarray
        The young-star velocity colormap data for use with 
        `matplotlib.axes.Axis.pcolormesh`. X data is along the 1 axis. Z data
        is along the 0 axis. The user can directly input this into
        `pcolormesh`.
    '''

    import os
    import h5py
    import numpy as np

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    from . import paths

    from . import rotate_galaxy
    from .rotate_galaxy import calculate_ang_mom, cal_rotation_matrix

    snapshot_times = np.loadtxt(
        paths.snap_times
    )
    time = float(snapshot_times[int(snap)][3])
    lbt = np.abs(time - 13.8)

    data = load_m12_data_olti(sim_path, snap)

    aux = data['temp'] < 1e4
    temp = data['temp'][aux]
    pos_gas = data['pos_gas'][aux]
    vel_gas = data['vel_gas'][aux]
    mass_gas = data['mass_gas'][aux]
    jnet_gas = calculate_ang_mom(mass_gas, pos_gas, vel_gas)

    r_matrix_gas = cal_rotation_matrix(jnet_gas, np.array((0.0, 0.0, 1.0)))
    pos_gas = rotate_galaxy.rotate(pos_gas, r_matrix_gas)
    vel_gas = rotate_galaxy.rotate(vel_gas, r_matrix_gas)

    v_gas = np.linalg.norm(vel_gas, axis=1)
    v_max = 220
    aux = v_gas <= v_max
    vel_gas = vel_gas[aux]
    pos_gas = pos_gas[aux]

    v_x_gas = vel_gas[:, 0]
    v_y_gas = vel_gas[:, 1]  # Use for colormap
    v_z_gas = vel_gas[:, 2]

    x_gas = pos_gas[:, 0]
    y_gas = pos_gas[:, 1]
    z_gas = pos_gas[:, 2]

    nbins_gas = 100
    # Create 2D histogram for gas
    hist_gas, x_edges_gas, z_edges_gas = np.histogram2d(
        x_gas,
        z_gas,
        bins=nbins_gas
    )
    hist_gas += 1  # Avoid log(0)

    # Apply the mask to keep bins with at least `gas_num` gas particles
    mask_gas = hist_gas >= gas_num

    #**********************************************************************
    # Bin the v_y values for gas and create a colormap based on the average 
    # v_y in each bin

    # Get indices of the x and z locations into which each particle falls
    x_bin_indices_gas = np.digitize(x_gas, x_edges_gas) - 1
    z_bin_indices_gas = np.digitize(z_gas, z_edges_gas) - 1
    v_y_colormap_gas = np.zeros_like(hist_gas)
    count_map_gas = np.zeros_like(hist_gas)

    for i in range(len(x_gas)):
        if (
                0 <= x_bin_indices_gas[i] < nbins_gas 
                and 0 <= z_bin_indices_gas[i] < nbins_gas):
            v_y_colormap_gas[
                x_bin_indices_gas[i],
                z_bin_indices_gas[i]
            ] += v_y_gas[i]
            count_map_gas[x_bin_indices_gas[i], z_bin_indices_gas[i]] += 1

    # Avoid division by zero:
    count_map_gas[count_map_gas == 0] = 1
    # Finally, calculating the avg v_y in each bin:
    v_y_colormap_gas /= count_map_gas
    #**********************************************************************

    # Apply the mask to remove bins with fewer than `gas_num` gas particles
    v_y_colormap_gas = np.where(mask_gas, v_y_colormap_gas, np.nan)

    young_mask = (data['sft'] <= (lbt + 0.5))
    pos_star = data['pos_star'][young_mask]
    vel_star = data['vel_star'][young_mask]

    r_matrix_star = cal_rotation_matrix(
        data['jnet_star'],
        np.array((0.0, 0.0, 1.0))
    )
    pos_star = rotate_galaxy.rotate(pos_star, r_matrix_star)
    vel_star = rotate_galaxy.rotate(vel_star, r_matrix_star)

    v_star = np.linalg.norm(vel_star, axis=1)
    aux = v_star <= v_max
    vel_star = vel_star[aux]
    pos_star = pos_star[aux]

    v_x_star = vel_star[:, 0]
    v_y_star = vel_star[:, 1]  # Use for colormap
    v_z_star = vel_star[:, 2]

    x_star = pos_star[:, 0]
    y_star = pos_star[:, 1]
    z_star = pos_star[:, 2]

    nbins_star = 100
    # Create 2D histogram for stars
    hist_star, x_edges_star, z_edges_star = np.histogram2d(
        x_star,
        z_star,
        bins=nbins_star
    )
    hist_star += 1  # Avoid log(0)

    # Apply the mask to keep bins with at least `star_num` star particles
    mask_star = hist_star >= star_num

    #**************************************************************************
    # Bin the v_y values for stars and create a colormap based on the average 
    # v_y in each bin
    x_bin_indices_star = np.digitize(x_star, x_edges_star) - 1
    z_bin_indices_star = np.digitize(z_star, z_edges_star) - 1
    v_y_colormap_star = np.zeros_like(hist_star)
    count_map_star = np.zeros_like(hist_star)

    for i in range(len(x_star)):
        if (
                0 <= x_bin_indices_star[i] < nbins_star 
                and 0 <= z_bin_indices_star[i] < nbins_star):
            v_y_colormap_star[
                x_bin_indices_star[i], 
                z_bin_indices_star[i]
            ]  += v_y_star[i]
            count_map_star[x_bin_indices_star[i], z_bin_indices_star[i]] += 1

    # Avoid division by zero:
    count_map_star[count_map_star == 0] = 1
    # Finally, calculating the avg v_y in each bin:
    v_y_colormap_star /= count_map_star
    #**************************************************************************

    # Apply the mask to remove bins with fewer than `star_num` star particles
    v_y_colormap_star = np.where(mask_star, v_y_colormap_star, np.nan)

    # Set up the figure and axis
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    # Get global vmin and vmax for the colormap based on both gas and stars
    #vmin = min(np.nanmin(v_y_colormap_gas), np.nanmin(v_y_colormap_star))
    vmax = max(np.nanmax(v_y_colormap_gas), np.nanmax(v_y_colormap_star))
    vmin = -1*vmax

    # The `H` output of `np.histogram2d` has x data along the 0 axis and y data 
    # along the 1 axis.
    # On one hand, that makes sense. However, if one prints out `H`, visually,
    # one might expect the x-axis to run horizontally along `H` and the y-axis 
    # to
    # run vertically along `H`. This is not the case. Additionally, 
    # `pcolormesh`
    # expects the column number as x and the row number as y 
    # Therefore, we must provide the transpose of `H` to
    # `pcolormesh`.
    v_y_colormesh_gas = v_y_colormap_gas.T
    v_y_colormesh_star = v_y_colormap_star.T

    # Plot gas with colormap based on v_y_gas
    pcol_gas = ax[0].pcolormesh(
        x_edges_gas,
        z_edges_gas,
        v_y_colormesh_gas,
        cmap=plt.cm.seismic_r,
        vmin=vmin,
        vmax=vmax
    )

    # Plot stars with colormap based on v_y_star
    pcol_star = ax[1].pcolormesh(
        x_edges_star,
        z_edges_star,
        v_y_colormesh_star,
        cmap=plt.cm.seismic_r,
        vmin=vmin,
        vmax=vmax
    )

    # Add colorbars for both plots
    fig.colorbar(pcol_gas, ax=ax[0], label=r'Gas LOS Velocity [kms$^{-1}]$')
    fig.colorbar(pcol_star, ax=ax[1], label=r'Star LOS Velocity [kms$^{-1}]$')

    # Add labels and text
    ax[0].text(
        0.9,
        0.9,
        'Gas',
        transform=ax[0].transAxes,
        fontsize=16,
        color='k',
        ha='right'
    )
    ax[1].text(
        0.9,
        0.9,
        'Young stars',
        transform=ax[1].transAxes,
        fontsize=16,
        color='k',
        ha='right'
    )

    ax[0].text(
        0.1,
        0.9,
        '{0}'.format(display_name),
        transform=ax[0].transAxes,
        color='k',
        fontsize=16
    )
    ax[0].text(
        0.1,
        0.85,
        'LBT = ' + str(np.round(lbt, 2)) + ' Gyr',
        transform=ax[0].transAxes,
        color='k',
        fontsize=14
    )

    # Axis labels
    ax[0].set_ylabel('Z [kpc]', fontsize=16)
    for i in range(len(ax)):
        ax[i].set_xlabel('X [kpc]', fontsize=16)
        ax[i].tick_params(
            axis='x',
            direction='in',
            pad=10,
            which='both',
            top=True,
            bottom=True,
            color='k',
            length = 6,
            width = 1.3
        )
        ax[i].tick_params(
            axis='y',
            direction='in',
            pad=10,
            which='both',
            left=True,
            right=True,
            color='k',
            length = 6,
            width = 1.3
        )
        ax[i].set_xlim(-17.5, 17.5)
        ax[i].set_ylim(-17.5, 17.5)

    # Tight layout and spacing
    plt.tight_layout()

    if save_plot:
        plt.savefig(os.path.join(
            paths.figures, 
            'vel_map_{0}_snap{1}.png'.format(display_name.lower(), snap)
        ))
    plt.show()

    return v_y_colormesh_gas, v_y_colormesh_star 
