from setuptools import setup, find_packages

setup(
    name='RFOnM',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'RFOnM=RFOnM.main:main',
        ],
    },
    install_requires=[
        'numpy',
        'scipy',
        'networkx',
    ],
    author='Xu-Wen Wang',
    description='RFOnM for disease module detection',
)
