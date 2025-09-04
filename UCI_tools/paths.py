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
    config_path = os.path.join(
        home, 
        env_str + 'config.ini'
    )
    base = os.path.join(home, 'uci-tools_output')
    if not os.path.isfile(config_path):
        defaults = (
            '[paths]\n'
            'figures = {base}'
        ).format(base=base)
        with open(config_path, 'w') as f:
            f.write(defaults)
    if not os.path.isdir(base):
        os.makedirs(base)
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
