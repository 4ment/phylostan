import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="phylostan",
    version="1.0.4",
    author="Mathieu Fourment",
    author_email="mathieu.fourment@uts.edu.au",
    description="Phylogenetic inference with Stan",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/4ment/phylostan",
    packages=['phylostan'],
    license='GPL3',
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=[
          'pystan>=2.19,<3',
          'DendroPy',
          'numpy>=1.7'
    ],
    entry_points = {
        'console_scripts': [ 'phylostan=phylostan.phylostan:main']
    }
)