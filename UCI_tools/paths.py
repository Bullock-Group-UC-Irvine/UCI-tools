import os
import configparser

env_str = os.getenv('CONDA_DEFAULT_ENV', 'base')
if env_str == 'base':
    env_str = ''
else:
    env_str += '_'
home = os.path.expanduser(os.path.join(
    '~/'
))

def ensure_user_config():
    config_fname = env_str + 'config.ini'
    config_path = os.path.join(
        home, 
        config_fname
    )
    if not os.path.isfile(config_path):
        output_dname = env_str + 'output'
        output_dir = os.path.join(home, output_dname)
        snap_times_path = (
            '/DFS-L/DATA/cosmo/grenache/omyrtaj/fofie/snapshot_times.txt'
        )
        defaults = (
            '[paths]\n'
            'figures = {0}\n'
            'snap_times = {1}'
        ).format(output_dir, snap_times_path)
        with open(config_path, 'w') as f:
            f.write(defaults)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        print(
            'NOTE: Before importing UCI_tools for the first time, you were'
            ' supposed to create a'
            ' config file in your home directory. However, no config file'
            ' existed up to now.'
            ' This code has created one for you in your home directory called' 
            ' {0}. It has also created the directory {1} in your home'
            ' directory. If you want to customize anything in the config file,'
            ' you can do so safely.'
            .format(config_fname, output_dname)
        )
    return config_path

def load_config():
    '''
    Check to see if something has created an environment variable specifying
    a config path (GitHub will do this when it runs workflows). If not, look
    for the user defined config. If that doesn't exist, create a default in the
    home dir.
    '''
    config_path = os.environ.get('MYPACKAGE_CONFIG')
    if config_path:
        config_path = os.path.exapnduser(config_path)
    else:
        config_path = ensure_user_config()

    config = configparser.ConfigParser()
    config.read(config_path)
    return config

config = load_config()
figures = config['paths']['figures']
snap_times = config['paths']['snap_times']
