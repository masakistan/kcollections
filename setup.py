#!/usr/bin/env python

import os, re, sys, sysconfig, platform, subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

__version__ = '1.0.2'

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[sourcedir + 'src/Kset.cc', sourcedir + 'src/Kdict.cc'])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                'CMake is required to install kcollections'
            )

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name))
        )
        cmake_args = ['-DPYTHON=ON',
                      '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        if self.debug:
            cfg = 'Debug'
        else:
            cfg = 'Release'
        #cfg = 'Debug'
        build_args = ['--config', cfg]

        if platform.system() == 'Windows':
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(), extdir
            )]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j4']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env
        )
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp
        )

setup(
    name='kcollections',
    version=__version__,
    author='M. Stanley Fujimoto and Cole A. Lyman',
    author_email='sfujimoto@gmail.com',
    url='https://github.com/masakistan/kcollections',
    license='GPLv3',
    description='A BloomFilterTrie implementation to be generally applicable for genomic applications.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=['kcollections'],
    ext_modules=[CMakeExtension('kcollections._Kdict'), CMakeExtension('kcollections._Kset'), CMakeExtension('kcollections._Kcounter')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False
)
