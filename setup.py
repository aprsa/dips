import setuptools
import dips

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name=dips.__name__,
    version=dips.__version__,
    author='Andrej Prsa',
    author_email='aprsa@villanova.edu',
    description='dips: detrending periodic signals',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/aprsa/dips',
    packages=setuptools.find_packages(),
    scripts=['bin/dips'],
    install_requires=[
        'numpy',
        'scipy >= 0.13.0',
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 3',
        'Environment :: Console',
        'Environment :: X11 Applications',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics'
    ],
)
