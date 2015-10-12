import os
from setuptools import setup, find_packages

setup(name='pyngspipe',
      version='0.0.1',
      description='pyngspipe is a pipeline for processing GEO NGS datasets based on the pyrnatools/pychiptools packages',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      packages=find_packages(),
      scripts=['scripts/pyngs_pipe.py', 'scripts/pyngs_report.py'],
      package_data={"pyngspipe":['data/*']},
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

pyngspipe is a pipeline for processing GEO NGS datasets based on the pyrnatools/pychiptools packages

 Contact
=============

If you have any questions or comments about pyngspipe, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
