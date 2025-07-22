from setuptools import setup, find_packages

setup(
    name='rfonm',
    version='0.1',
    packages=['rfonm'],
    entry_points={
        'console_scripts': [
            'rfonm=rfonm.main:main',
        ],
    },
    install_requires=[
    'numpy>=1.26,<2.0',
    'networkx>=3.4,<4.0',
    ],
    author='Xu-Wen Wang',
    description='RFOnM for disease module detection',
)
