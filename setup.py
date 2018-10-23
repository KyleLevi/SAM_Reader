from setuptools import setup

setup(name='sam_reader',
    version='0.1b2',
    python_version='3+',
    author='Kyle Levi',
    description='A package for reading many SAM/BAM files at once using Pysam',
    author_email='pip@kylelevi.com',
    install_requires=['pysam'],
    url='https://github.com/KyleLevi/SAM_Reader')