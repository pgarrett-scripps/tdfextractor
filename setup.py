from setuptools import setup

setup(
    name='tdfextractor',
    version='0.1.0',
    description='extract mass spec files from bruker d folder',
    url='https://github.com/pgarrett-scripps/tdfextractor.git',
    author='Patrick Garrett',
    author_email='pgarrett@scripps.edu',
    license='MIT',
    packages=['tdfextractor'],
    install_requires=['tdfpy~=0.1.0'],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.10',
    ],
)