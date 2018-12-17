import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="myvariant",
    version="1.0.0",
    author="Chunlei Wu, Cyrus Afrasiabi, Sebastien Lelong",
    author_email="cwu@scripps.edu",
    description="Python Client for MyVariant.Info services.",
    license="BSD",
    keywords="biology variant annotation web service client api myvariant",
    url="https://github.com/biothings/myvariant.py",
    packages=['myvariant'],
    long_description=read('README.rst'),
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Utilities",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=[
        'biothings_client>=0.2.0'
    ],
)
