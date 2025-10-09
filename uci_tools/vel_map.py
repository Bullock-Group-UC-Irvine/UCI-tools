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

def load_m12_data_olti(sim_path, snap, xmax=None, zmax=None):
    '''
    Load Olti's data for use with `uci_tools.vel_map.plot`. You can also
    use this with your own data at `sim_path` as long as it's an hdf5 file with
    the following `h5py.Dataset`s:
        'gas_coord_unroated'
        'gas_vel_unrotated'
        'jnet_gas'
        'gas_temp'
        'mass_gas'
        'star_coord_unrotated'
        'star_vel_unrotated',
        'sft_Gyr',
        'jnet_young_star'

    Parameters
    ----------
    sim_path: str 
        The path to the simulation the user wants to analyze. The user could
        use uci_tools.vel_map.get_m12_path_olti to easily generate a path
        that leads to Olti's files, or they could supply their own path.
        Another option would be for the user to write their own `get_m12_path`
        method in uci_tools.vel_map and make a pull request so everyone
        has it.
    snap: str
        The snapshot number corresponding to `sim_path`. It should be in string
        format with three digits. The code uses this to determine the scale
        factor and thereby physical distances.
    xmax: float, default None
        The absolute value of the rotated x-axis distance from the host center
        that a particle must be at or below for the code to include it.
    zmax: float, default None
        The absolute value of the rotated z-axis distance from the host center
        that a particle must be at or below for the code to include it.

    Returns
    -------
    pos_star: np.ndarray, shape (N_stars, 3)
        The centered, rotated position vectors of the
        simulation's star particles in Cartesian coordinates in physical kpc.
        The function rotated them
        so the the x- and y-axes are in the plane of the disc by aligning 
        the z-axis 
        with the net angular momentum of the young stars.
    vel_star: np.ndarray, shape (N_stars, 3)
        The rotated velocity vectors of the star particle in Cartesian
        coordinates
        relative to
        the host center. The function rotated them
        so the the x- and y-axes are in the plane of the disc by aligning 
        the z-axis 
        with the net angular momentum of the cold gas (T <= 1e4 K)
    pos_gas: np.ndarray, shape (N_gas, 3)
        The centered, rotated position vectors of the
        simulation's gas particles in Cartesian coordinates in physical kpc.
        The function rotated them
        so the the x- and y-axes are in the plane of the disc by aligning 
        the z-axis 
        with the net angular momentum of the cold gas 
        (T <= 1e4 K)
    vel_gas: np.ndarray, shape (N_gas, 3)
        The rotated velocity vectors of the gas particle in Cartesian
        coordinates
        relative to
        the host center. The function rotated them
        so the the x- and y-axes are in the plane of the disc by aligning 
        the z-axis 
        with the net angular momentum of the cold gas (T <= 1e4 K)
    '''
    import h5py
    import numpy as np
    from . import config 
    from . import rotate_galaxy

    if xmax is None:
        xmax = np.inf
    if zmax is None:
        zmax = np.inf

    snapshot_times = np.loadtxt(
        config.config['uci_tools_paths']['snap_times']
    )
    time = float(snapshot_times[int(snap)][3])
    lbt = np.abs(time - 13.8)
    a = float(snapshot_times[int(snap)][1])

    data = {}
    with h5py.File(sim_path, 'r') as f:
        host_center = np.array(f['host_center'])
        host_vel = np.array(f['host_velocity'])

        # Load gas data
        pos_gas = a * (
            np.array(f['gas_coord_unrotated']) - host_center
        )
        vel_gas = np.array(f['gas_vel_unrotated']) - host_vel
        jnet_gas = np.array(f['jnet_gas'])
        temp = np.array(f['gas_temp'])
        mass_gas = np.array(f['mass_gas'])
        
        # Load star data
        pos_star = a * (
            np.array(f['star_coord_unrotated']) - host_center
        )
        vel_star = np.array(f['star_vel_unrotated']) - host_vel
        sft = np.array(f['sft_Gyr'])
        jnet_star = np.array(f['jnet_young_star'])

    aux = temp < 1e4
    temp = temp[aux]
    pos_gas = pos_gas[aux]
    vel_gas = vel_gas[aux]
    mass_gas = mass_gas[aux]
    jnet_gas = rotate_galaxy.calculate_ang_mom(mass_gas, pos_gas, vel_gas)

    # Rotation matrix
    r_matrix_gas = rotate_galaxy.cal_rotation_matrix(
        jnet_gas,
        np.array((0.0, 0.0, 1.0))
    )
    pos_gas = rotate_galaxy.rotate(pos_gas, r_matrix_gas)
    vel_gas = rotate_galaxy.rotate(vel_gas, r_matrix_gas)

    gas_in_x = np.abs(pos_gas[:, 0]) <= xmax
    gas_in_z = np.abs(pos_gas[:, 2]) <= zmax

    v_gas = np.linalg.norm(vel_gas, axis=1)
    v_max = 220
    aux = v_gas <= v_max
    vel_gas = vel_gas[aux & gas_in_x & gas_in_z]
    pos_gas = pos_gas[aux & gas_in_x & gas_in_z]

    young_mask = (sft <= (lbt + .5))
    pos_star = pos_star[young_mask]
    vel_star = vel_star[young_mask]

    # Rotation matrix
    r_matrix_star = rotate_galaxy.cal_rotation_matrix(
        jnet_star,
        np.array((0.0, 0.0, 1.0))
    )
    pos_star = rotate_galaxy.rotate(pos_star, r_matrix_star)
    vel_star = rotate_galaxy.rotate(vel_star, r_matrix_star)

    stars_in_x = np.abs(pos_star[:, 0]) <= xmax
    stars_in_z = np.abs(pos_star[:, 2]) <= zmax

    v_star = np.linalg.norm(vel_star, axis=1)
    aux = v_star <= v_max
    vel_star = vel_star[aux & stars_in_x & stars_in_z]
    pos_star = pos_star[aux & stars_in_x & stars_in_z]

    return pos_star, vel_star, pos_gas, vel_gas 

def plot(
        pos_star,
        vel_star,
        pos_gas,
        vel_gas,
        display_name,
        snap,
        horiz_axis=0,
        vert_axis=2,
        res=100,
        gas_num=50,
        star_num=20,
        save_plot=False,
        show_plot=True):
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
    pos_gas: np.ndarray, shape (N_gas, 3)
        Cartesian position vectors relative to the host center, in physical
        kpc, of the
        gas particles whose velocities the user wants to map.
        The resulting velocity map assumes the line of sight is down the 1 axis
        (or y-axis). If the user provides rotated vectors with their z-axis
        aligned with the galaxy's net angular momentum, this is analagous to
        looking into the disc.
    vel_gas: np.ndarray, shape (N_gas, 3)
        Cartesian velocity vectors relative to the host center, in physical
        km/s, of the gas particles whose
        velocities the user wants to map.
        The resulting velocity map assumes the line of sight is down the 1 axis
        (or y-axis). Therefore, this function maps the 1-axis velocities. If
        the user provides rotated vectors with their z-axis
        aligned with the galaxy's net angular momentum, this is analagous to
        looking into the disc.
    pos_star: np.ndarray, shape (N_stars, 3)
        Cartesian position vectors relative to the host center, in physical
        kpc, of the star
        particles whose velocities the user wants to map.
        The resulting velocity map assumes the line of sight is down the 1 axis
        (or y-axis). If the user provides rotated vectors with their z-axis
        aligned with the galaxy's net angular momentum, this is analagous to
        looking into the disc.
    vel_star: np.ndarray, shape (N_stars, 3)
        Cartesian velocity vectors relative to the host center, in physical
        km/s, of the gas particles whose
        velocities the user wants to map.
        The resulting velocity map assumes the line of sight is down the 1 axis
        (or y-axis). Therefore, this function maps the 1-axis velocities. If
        the user provides rotated vectors with their z-axis
        aligned with the galaxy's net angular momentum, this is analagous to
        looking into the disc.
    display_name: str 
        Simulation name to show in the plot.
    snap: str
        The snapshot number corresponding to the data. It should be in string
        format with three digits. The code uses this to
        display the correct look-back time in the plot.
    horiz_axis: int, default 0
        The index of the axis to show horizontally on the plot. The default is 
        the 0
        or x-axis.
    vert_axis: int, default 2
        The index of the axis to show vertically on the plot. The default is
        the 2 or z-axis.
    res: int, default 100
        Resolution: the number of pixels (i.e. bins) in each axis.
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
    show_plot: bool, default True
        Whether to display the velocity map.

    Returns
    -------
    velmap_gas: np.ndarray
        The gas velocity colormap data for use with 
        `matplotlib.axes.Axis.pcolormesh`. X data is along the 1 axis. Z data
        is along the 0 axis. The user can directly input this into
        `pcolormesh`.
    velmap_star: np.ndarray
        The young-star velocity colormap data for use with 
        `matplotlib.axes.Axis.pcolormesh`. X data is along the 1 axis. Z data
        is along the 0 axis. The user can directly input this into
        `pcolormesh`.
    x_edges_gas
    z_edges_gas
    x_edges_star
    z_edges_star
    quadmesh_gas
    quadmesh_star
    '''

    import os
    import h5py
    import numpy as np

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    import astropy.cosmology as cosmo
    import astropy

    from . import config

    snapshot_times = np.loadtxt(
        config.config['uci_tools_paths']['snap_times']
    )
    time = float(snapshot_times[int(snap)][3])
    lbt = np.abs(time - 13.8)

    axes = [0, 1, 2]
    # Determine which axis the user did not specify as a projection axis 
    los_axis = np.setdiff1d(axes, [horiz_axis, vert_axis])
    if len(los_axis) > 1:
        # There should only be one line-of-sight axis
        raise ValueError('Something is wrong with the axis specifications')
    los_axis = los_axis[0]

    # Only the `los_axis`-axis velocity is necessary.
    v_y_gas = vel_gas[:, los_axis]  # Use for colormap
    #v_x_gas = vel_gas[:, horiz_axis]
    #v_z_gas = vel_gas[:, vert_axis]

    # The function looks down the `los_axis`-axis and so does not use its
    # positional
    # information.
    x_gas = pos_gas[:, horiz_axis]
    #y_gas = pos_gas[:, los_axis]
    z_gas = pos_gas[:, vert_axis]

    # Create 2D histogram for gas
    hist_gas, x_edges_gas, z_edges_gas = np.histogram2d(
        x_gas,
        z_gas,
        bins=res
    )
    # Need to reverse z_edges so that z_edges[0] is the highest z. This is what
    # matplotlib.pyplot.imshow expects, and it makes sense when one considers
    # that when printing the mesh data, high z should be on the top of the
    # matrix.
    z_edges_gas = z_edges_gas[::-1]

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
                0 <= x_bin_indices_gas[i] < res 
                and 0 <= z_bin_indices_gas[i] < res):
            v_y_colormap_gas[
                x_bin_indices_gas[i],
                z_bin_indices_gas[i]
            ] += v_y_gas[i]
            count_map_gas[x_bin_indices_gas[i], z_bin_indices_gas[i]] += 1

    # Avoid division by zero:
    count_map_gas[count_map_gas == 0] = 1
    # Finally, calculating the avg v_y in each bin:
    v_y_colormap_gas /= count_map_gas

    # Apply the mask to remove bins with fewer than `gas_num` gas particles
    v_y_colormap_gas = np.where(mask_gas, v_y_colormap_gas, np.nan)
    #**********************************************************************

    # Only the `los_axis`-axis velocity is necessary.
    v_y_star = vel_star[:, los_axis]  # Use for colormap
    #v_x_star = vel_star[:, horiz_axis]
    #v_z_star = vel_star[:, vert_axis]

    # The function looks down the `los_axis`-axis and so does not use its
    # positional
    # information.
    x_star = pos_star[:, horiz_axis]
    #y_star = pos_star[:, los_axis]
    z_star = pos_star[:, vert_axis]

    # Create 2D histogram for stars
    hist_star, x_edges_star, z_edges_star = np.histogram2d(
        x_star,
        z_star,
        bins=res
    )
    # Need to reverse z_edges so that z_edges[0] is the highest z. This is what
    # matplotlib.pyplot.imshow expects, and it makes sense when one considers
    # that when printing the mesh data, high z should be on the top of the
    # matrix.
    z_edges_star = z_edges_star[::-1]

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
                0 <= x_bin_indices_star[i] < res 
                and 0 <= z_bin_indices_star[i] < res):
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
    velmap_gas = v_y_colormap_gas.T
    velmap_star = v_y_colormap_star.T

    # Plot gas with colormap based on v_y_gas
    quadmesh_gas = ax[0].pcolormesh(
        x_edges_gas,
        z_edges_gas,
        velmap_gas,
        cmap=plt.cm.seismic_r,
        vmin=vmin,
        vmax=vmax
    )

    # Plot stars with colormap based on v_y_star
    quadmesh_star = ax[1].pcolormesh(
        x_edges_star,
        z_edges_star,
        velmap_star,
        cmap=plt.cm.seismic_r,
        vmin=vmin,
        vmax=vmax
    )

    if show_plot or save_plot:
        # Add colorbars for both plots
        fig.colorbar(
            quadmesh_gas,
            ax=ax[0],
            label=r'Gas LOS Velocity [kms$^{-1}]$'
        )
        fig.colorbar(
            quadmesh_star,
            ax=ax[1],
            label=r'Star LOS Velocity [kms$^{-1}]$'
        )

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
            'Stars',
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
        axis_labels = ['$x$', '$y$', '$z$']
        ax[0].set_ylabel(
            '{0} [kpc]'.format(axis_labels[vert_axis]),
            fontsize=16
        )
        for i in range(len(ax)):
            ax[i].set_xlabel(
                '{0} [kpc]'.format(axis_labels[horiz_axis]),
                fontsize=16
            )
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

        # Tight layout and spacing
        plt.tight_layout()

        if save_plot:
            plt.savefig(os.path.join(
                config.config['uci_tools_paths']['output'], 
                'vel_map_{0}_snap{1}.png'.format(display_name.lower(), snap)
            ))
        if show_plot:
            plt.show()
    else:
        plt.close()

    return (
        velmap_gas,
        velmap_star,
        x_edges_gas,
        z_edges_gas,
        x_edges_star,
        z_edges_star,
        quadmesh_gas,
        quadmesh_star,
    )
