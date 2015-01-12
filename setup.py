import os
from setuptools import setup, find_packages

setup(name='pyrnapipe',
      version='0.0.1',
      description='pyrnapipe is a pipeline for processing GEO RNA-seq datasets based on the pyrnatools package',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      packages=find_packages(),
      package_data={"pyrnapipe":['data/*']},
      scripts=['scripts/pyrna_pipe.py'],
      install_requires=['pysam', 'pybedtools'],
      license='GPLv3',
      platforms='any',
      classifiers=[
         'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
         'Development Status :: 3 - Alpha',
         'Programming Language :: Python :: 2.7',
         'Environment :: Console',
      ],
      long_description="""

pyrnapipe is a pipeline for processing GEO RNA-seq datasets based on the pyrnatools package

 Contact
=============

If you have any questions or comments about pyrnapipe, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
