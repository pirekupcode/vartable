from glob import glob
from setuptools import setup, find_packages
import vartable
# Run setuptools setup
setup(
    name = vartable.__projectname__,
    version = vartable.__version__,
    packages = find_packages(),
    scripts = glob('bin/*'),
    entry_points = {
        'console_scripts': [
            'vartable_report = vartable.vartable:main'
          ]
        } )

