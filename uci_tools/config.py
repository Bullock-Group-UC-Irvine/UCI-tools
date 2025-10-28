import os
import configparser

env_str = os.getenv('CONDA_DEFAULT_ENV', 'base')
if env_str == 'base':
    env_suffix = ''
else:
    env_suffix = '_' + env_str
home = os.path.expanduser(os.path.join(
    '~/'
))
config_fname = 'config' + env_suffix + '.ini'

def ensure_user_config():
    config = configparser.ConfigParser()
    config_path = os.path.join(
        home, 
        config_fname
    )

    if os.path.isfile(config_path):
        config.read(config_path)
    else:
        print(
            'Before importing {2} for the first time, you were'
            ' supposed to create a'
            ' config file for your environment and/or add a {3}_paths'
            ' section to your existing config file. However, no config file'
            ' existed up to now.'
            ' This code has created one for you in your home directory called' 
            ' {0}.'
            ' If you want to customize anything in the config file,'
            ' you can do so safely.\n'
            .format(config_fname, None, __package__, env_str)
        )

    if not config.has_section(f'{__package__}_paths'):
        config.add_section(f'{__package__}_paths')
        print(
            f'{__package__}_paths section added to {config_fname}'
            f'\n\nNOTE: Anything this code adds to {__package__}_paths assumes'
            f' you are on UC Irvine\'s Greenplanet cluster in the cosmo'
            f' group. If you are not, your code will not work with these'
            f' paths,'
            f' and you must properly configure {config_fname}.\n'
        )

    if not config.has_option(f'{__package__}_paths', 'output_dir'):
        output_dname = 'output' + env_suffix
        output_dir = os.path.join(home, output_dname)
        config.set(f'{__package__}_paths', 'output_dir', output_dir)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
            print(f'{output_dir} created')
        print('output_dir added to gallearn_paths')

    if not config.has_option(f'{__package__}_paths', 'snap_times'):
        snap_times_path = (
            '/DFS-L/DATA/cosmo/grenache/omyrtaj/fofie/snapshot_times.txt'
        )
        config.set(f'{__package__}_paths', 'snap_times', snap_times_path)
        print(f'snap_times added to {__package__}_paths')

    with open(config_path, 'w') as f:
        config.write(f)

    return config_path

def load_config():
    '''
    Check to see if something has created an environment variable specifying
    a config path (GitHub will do this when it runs workflows). If not, look
    for the user defined config. If that doesn't exist, create a default in the
    home dir.
    '''
    # If CI set the path to the config_ci.ini, the following will be something
    # other than `False`.
    config_path = os.environ.get('MYPACKAGE_CONFIG')
    if config_path:
        # In that case, go ahead and get the full path the .ini
        config_path = os.path.expanduser(config_path)
    else:
        # Otherwise, check that the user has a config_<environment_name>.ini.
        # Create one for them if they don't.
        config_path = ensure_user_config()

    config = configparser.ConfigParser()
    config.read(config_path)
    return config

config = load_config()
