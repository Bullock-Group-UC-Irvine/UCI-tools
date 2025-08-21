import os
import configparser

env_str = os.getenv('CONDA_DEFAULT_ENV', 'base')
if env_str == 'base':
    env_str = ''
else:
    env_str += '_'
config_path = os.path.expanduser(os.path.join(
    '~/',
    env_str + 'config.ini'
))
config = configparser.ConfigParser()
config.read(config_path)

figures = config['paths']['figures']
