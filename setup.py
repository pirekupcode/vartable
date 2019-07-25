from glob import glob
from setuptools import setup, find_packages

# Run setuptools setup
setup( 
    packages = find_packages(),
    scripts = glob('bin/*'),
    entry_points = {
        'console_scripts': [
            'vartable_report = vartable.vartable:main'
          ]
        } )
    
