import setuptools

# Get the version number.
v_dict = {}
exec(open('uci_tools/__version__.py').read(), v_dict)

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="uci_tools",
    version=v_dict['__version__'],
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
            'uci-tools-refactor = uci_tools.refactor:main',
        ],
    }
)
