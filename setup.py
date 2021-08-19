from setuptools import setup

setup(
    name='mdpkg',
    version='0.0.1',
    author='Luís H. Carnevale',
    author_email='lh.carnevale@gmail.com',
    packages=['mdpkg'],
    url='https://github.com/luishcc/mdpkg',
    license='LICENSE.txt',
    description='For Molecular Dynamics analysis',
    long_description=open('README.md').read(),
    install_requires=[
       "numpy",
    ],
    python_requires='>=3.7'
    )
