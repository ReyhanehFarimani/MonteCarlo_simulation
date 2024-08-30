from setuptools import setup, find_packages

setup(
    name='gcmc_post_processing',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'pandas',
        'scipy',
    ],
)