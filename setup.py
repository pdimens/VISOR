from setuptools import setup,find_packages

setup(name='VISOR',
  version='1.2.0',
  description='VarIant SimulatOR',
  url='https://github.com/davidebolo1993/VISOR.git',
  requires=['python (>= 3.6)'],
  author='Davide Bolognini',
  author_email='davidebolognini7@gmail.com',
  license='LICENSE.txt',
  dependency_links=['https://github.com/rrwick/Badread/tarball/master#egg=Badread'],
  install_requires=['pyfaidx', 'pysam', 'pywgsim', 'pybedtools', 'mappy', 'plotly==3.10.0', 'numpy'],
  zip_safe=False,
  packages=find_packages(),
  include_package_data=True,
  entry_points={'console_scripts': ['VISOR=VISOR.VISOR:main']}          
)
