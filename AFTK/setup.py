#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 1/26/23

"""

from setuptools import Extension, setup  # must be on top

import numpy
from Cython.Build import cythonize
from pathlib import Path


cython_compiler_directives = {
    'language_level': '3',
    'embedsignature': True,
    'cdivision': True,
    'boundscheck': False,
    'wraparound': False,
    'profile': True
}
"""Refer to
https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#compiler-directives
"""



# ========== ========== ========== ========== ========== ==========
options = {
    'include_dirs': [numpy.get_include(),
                     '/home/rafael/.pyenv/versions/3.11.1/envs/aftk/lib/python3.11/site-packages/scipy'],
    'define_macros': [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
}

extensions = [
    Extension('__test.cylabs', ['__test/cylabs.pyx'], **options),

    Extension('aftk.mathtools.psdm', ['aftk/mathtools/psdm.pyx'], **options),
    Extension('aftk.mathtools.spharm', ['aftk/mathtools/spharm.pyx'], **options),
    Extension('aftk.mathtools.sptri', ['aftk/mathtools/sptri.pyx'], **options),
    # Extension('aftk.mathtools.linreg', ['aftk/mathtools/linreg.pyx'], **options),

    Extension('aftk.plotting.misc', ['aftk/plotting/misc.pyx'], **options),
    Extension('aftk.plotting.smc', ['aftk/plotting/smc.pyx'], **options),

    Extension('aftk.antenna.eigenantennas', ['aftk/antenna/eigenantennas.pyx'], **options),
    Extension('aftk.antenna.smc', ['aftk/antenna/smc.pyx'], **options),

    Extension('aftk.utils.misc', ['aftk/utils/misc.pyx'], **options),
    Extension('aftk.utils.plotting', ['aftk/utils/plotting.pyx'], **options),
    Extension('aftk.model.eigenfunctions', ['aftk/model/eigenfunctions.pyx'], **options),
    Extension('aftk.measurements.ludwig', ['aftk/measurements/ludwig.pyx'], **options),
    Extension('aftk.measurements.radiation_field', ['aftk/measurements/radiation_field.pyx'], **options),
    # Extension('aftk.identification.ctools', ['aftk/identification/ctools.pyx'], **options)
    # Extension('smet.satellite.propagator.sgp4', ['smet/satellite/propagator/sgp4.pyx'] + sgp4_src_file, **options),
    # Extension('smet.sofa.__init__', ['smet/sofa/__init__.pyx'] + sofa_src_files, **options),
    # Extension('smet.iers.__init__', ['smet/iers/__init__.pyx'] + sofa_src_files, **options),
    # # Extension('smet/**/*', ['smet/**/*.pyx'], include_dirs=[numpy.get_include()]),
    # Extension('smet.utils.custom_exceptions', ['smet/utils/custom_exceptions.pyx'], **options),
    # Extension('smet.utils.convenience', ['smet/utils/convenience.pyx'], **options),
    # Extension('smet.utils.parsers', ['smet/utils/parsers.pyx'], **options),
    # Extension('smet.utils.time', ['smet/utils/time.pyx'] + sofa_src_files, **options),
]

# ========== ========== ========== ========== ========== ==========
setup(ext_modules=cythonize(extensions,
                            annotate=True,
                            compiler_directives=cython_compiler_directives,
                            include_path=[numpy.get_include()]
                            ))



