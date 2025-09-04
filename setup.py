import setuptools
from UCI_tools import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="UCI_tools",
    version=__version__,
    author="UCI Bullock Galaxy Dark Matter Group",
    author_email="bullock@uci.edu",
    description="Software tools developed by and commonly used by the group.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Bullock-Group-UC-Irvine/shared-analysis-tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'numpy',
        'pandas',
        'pytest',
        'h5py',
        'progressbar',
        'nbformat',
        'astropy',
        'tqdm'
    ],
    entry_points={
        'console_scripts': [
            'UCI-tools.refactor = UCI_tools.refactor:main',
        ],
    }
)
