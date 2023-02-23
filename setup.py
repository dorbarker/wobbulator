from setuptools import setup, find_packages
from wobbulator import __author__, __author_email__, __version__

setup(
    name="wobbulator",
    version=__version__,
    license='GPL-3.0-or-later',
    author=__author__,
    author_email=__author_email__,
    packages=find_packages(),
    install_requires=["biopython>=1.78"],
    entry_points={"console_scripts": ["wobbulator=wobbulator.wobbulator:main"]},
)
