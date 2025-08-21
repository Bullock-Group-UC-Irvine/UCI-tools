import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from .rotate_galaxy import calculate_ang_mom, calculate_rotation_matrix

snap = ['153']

def get_m12_path(simname, snap):
    '''
    Parameters
    ----------
    simname: {
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
    snap: str
        The snapshot number to load. The snapshot number should be in
        string format.

    Returns
    -------
    path: str
        The path to the file being analyzed.
    '''
    import os
    path = os.path.join(
        '/DFS-L/DATA/cosmo/grenache/omyrtaj/analysis_data/metaldiff/',
        simname,
        'id_jnet_jzjc_jjc_'
            + snap
            + '_host0_20kpc_rockstar_centers_metaldiff.hdf5'
    )
    return path

def plot(
        simname, 
        snap,
        gas_num=50,
        star_num=20):
    '''
    Parameters
    ----------
    simname: {
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
    snap: str
        The snapshot number to load. The snapshot number should be in
        string format.
    gas_num: int
        The minimum number of gas particles a 2d histogram bin must have
        for it to be included.
    star_num: int
        The minimum number of star particles a 2d histogram bin must have for
        it to be included.
    '''

    from . import paths

    snapshot_times = np.loadtxt(
        '/DFS-L/DATA/cosmo/grenache/omyrtaj/fofie/snapshot_times.txt'
    )
    time = float(snapshot_times[int(snap)][3])
    a = float(snapshot_times[int(snap)][1])
    lbt = np.abs(time - 13.8)

    # Load gas data
    path = get_m12_path(simname, snap)
    data = h5py.File(
        path,
        'r'
    )

    host_center = np.array(data['host_center'])
    host_vel = np.array(data['host_velocity'])
    pos_gas = a * (np.array(data['gas_coord_unrotated']) - host_center)
    vel_gas = np.array(data['gas_vel_unrotated']) - host_vel
    jnet_gas = np.array(data['jnet_gas'])
    temp = np.array(data['gas_temp'])
    aux = temp < 1e4
    temp = temp[aux]
    pos_gas = pos_gas[aux]
    vel_gas = vel_gas[aux]
    mass_gas = np.array(data['mass_gas'])[aux]
    jnet_gas = calculate_ang_mom(mass_gas,pos_gas,vel_gas)

    r_matrix_gas = cal_rotation_matrix(jnet_gas, np.array((0.0, 0.0, 1.0)))
    pos_gas = rotate_matrix(pos_gas, r_matrix_gas)
    vel_gas = rotate_matrix(vel_gas, r_matrix_gas)

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

    # Bin the v_y values for gas and create a colormap based on the average v_y in each bin
    x_bin_indices_gas = np.digitize(x_gas, x_edges_gas) - 1
    z_bin_indices_gas = np.digitize(z_gas, z_edges_gas) - 1
    v_y_colormap_gas = np.zeros_like(hist_gas)
    count_map_gas = np.zeros_like(hist_gas)

    for i in range(len(x_gas)):
        if 0 <= x_bin_indices_gas[i] < nbins_gas and 0 <= z_bin_indices_gas[i] < nbins_gas:
            v_y_colormap_gas[x_bin_indices_gas[i], z_bin_indices_gas[i]] += v_y_gas[i]
            count_map_gas[x_bin_indices_gas[i], z_bin_indices_gas[i]] += 1

    # Avoid division by zero
    count_map_gas[count_map_gas == 0] = 1
    v_y_colormap_gas /= count_map_gas

    # Apply the mask to remove bins with fewer than 1000 gas particles
    v_y_colormap_gas = np.where(mask_gas, v_y_colormap_gas, np.nan)

    # Load star data
    pos_star = a * (np.array(data['star_coord_unrotated']) - host_center)
    vel_star = np.array(data['star_vel_unrotated']) - host_vel
    sft = np.array(data['sft_Gyr'])
    young_mask = (sft <= (lbt + 0.5))
    pos_star = pos_star[young_mask]
    vel_star = vel_star[young_mask]
    jnet_star = np.array(data['jnet_young_star'])

    r_matrix_star = cal_rotation_matrix(jnet_star, np.array((0.0, 0.0, 1.0)))
    pos_star = rotate_matrix(pos_star, r_matrix_star)
    vel_star = rotate_matrix(vel_star, r_matrix_star)

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
    hist_star, x_edges_star, z_edges_star = np.histogram2d(x_star, z_star, bins=nbins_star)
    hist_star += 1  # Avoid log(0)

    # Apply the mask to keep bins with at least `star_num` star particles
    mask_star = hist_star >= star_num

    # Bin the v_y values for stars and create a colormap based on the average v_y in each bin
    x_bin_indices_star = np.digitize(x_star, x_edges_star) - 1
    z_bin_indices_star = np.digitize(z_star, z_edges_star) - 1
    v_y_colormap_star = np.zeros_like(hist_star)
    count_map_star = np.zeros_like(hist_star)

    for i in range(len(x_star)):
        if 0 <= x_bin_indices_star[i] < nbins_star and 0 <= z_bin_indices_star[i] < nbins_star:
            v_y_colormap_star[x_bin_indices_star[i], z_bin_indices_star[i]] += v_y_star[i]
            count_map_star[x_bin_indices_star[i], z_bin_indices_star[i]] += 1

    # Avoid division by zero
    count_map_star[count_map_star == 0] = 1
    v_y_colormap_star /= count_map_star

    # Apply the mask to remove bins with fewer than 1000 star particles
    v_y_colormap_star = np.where(mask_star, v_y_colormap_star, np.nan)

    # Set up the figure and axis
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    # Get global vmin and vmax for the colormap based on both gas and stars
    #vmin = min(np.nanmin(v_y_colormap_gas), np.nanmin(v_y_colormap_star))
    vmax = max(np.nanmax(v_y_colormap_gas), np.nanmax(v_y_colormap_star))
    vmin = -1*vmax

    # Plot gas with colormap based on v_y_gas
    pcol_gas = ax[0].pcolormesh(x_edges_gas, z_edges_gas, v_y_colormap_gas.T, cmap=plt.cm.seismic_r, vmin=vmin, vmax=vmax)

    # Plot stars with colormap based on v_y_star
    pcol_star = ax[1].pcolormesh(x_edges_star, z_edges_star, v_y_colormap_star.T, cmap=plt.cm.seismic_r, vmin=vmin, vmax=vmax)

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
        0.95,
        'Thelma',
        transform=ax[0].transAxes,
        color='k',
        fontsize=16
    )
    ax[0].text(
        0.1,
        0.9,
        'LBT = ' + str(np.round(lbt[0], 2)) + ' Gyr',
        transform=ax[0].transAxes,
        color='k',
        fontsize=14
    )

    # Axis labels
    ax[0].set_ylabel('Z [kpc]', fontsize=16)
    for i in range(len(ax)):
        ax[i].set_xlabel('X [kpc]', fontsize=16)
        ax[i].tick_params(axis='x', direction='in', pad=10, which='both', top=True, bottom=True, color='k', length = 6, width = 1.3)
        ax[i].tick_params(axis='y', direction='in', pad=10, which='both', left=True, right=True, color='k', length = 6, width = 1.3)
        ax[i].set_xlim(-17.5, 17.5)
        ax[i].set_ylim(-17.5, 17.5)

    # Tight layout and spacing
    plt.tight_layout()

    plt.savefig(os.path.join(paths.figures, 'plot_vel_map_thelma_3.pdf'))
    plt.show()

    return None

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Make a velocity map for the given simulation and snapshot."
        )
    )

    parser.add_argument(
        "-s", "--simname",
        required=True,
        help="Name of the simulation."
    )

    parser.add_argument(
        "-n", "--snap",
        type=str,
        required=True,
        help="Snapshot number in string format."
    )

    parser.add_argument(
        "-g", "--gas_num",
        type=int,
        required=True,
        help=(
            'Number of gas particles that must be in a 2d histogram bin in'
            ' order for the plot to include it.'
        )
    )

    parser.add_argument(
        "-t", "--star_num",
        type=int,
        required=True,
        help=(
            'Number of star particles that must be in a 2d histogram bin in'
            ' order for the plot to include it.'
        )
    )

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    plot(args.simname, args.snap, args.gas_num, args.star_num)
