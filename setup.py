import os
import sys
from glob import glob
import re
import tarfile
import subprocess
import shutil
import zipfile
import tempfile
import re
import time
import configparser

#try:
#    import ConfigParser as configparser
#except ImportError:
#    import configparser

config = configparser.ConfigParser()
install_folder=os.path.dirname(os.path.realpath(__file__))
config_file=os.path.join(install_folder,"configs/baqlava.cfg")
config.read(config_file)
setup_directory = os.path.abspath(os.path.dirname(__file__))

VERSION = config.get("metadata","version",fallback="0.5")
AUTHOR = "BAQLaVa Development Team"
MAINTAINER = "Jordan Jensen"
MAINTAINER_EMAIL = "jordanjensen@g.harvard.edu"


# check for either of the required versions
# required python versions (3.10+)
required_python_version_major = [3]
required_python_version_minor = [10]

pass_check = False
try:
	for major, minor in zip(required_python_version_major, required_python_version_minor):
		if (sys.version_info[0] == major and sys.version_info[1] >= minor):
			pass_check = True
except (AttributeError, IndexError):
	sys.exit("CRITICAL ERROR: The python version found (version 1) " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")

if not pass_check:
	sys.exit("CRITICAL ERROR: The python version found (version " +
	         str(sys.version_info[0]) + "." + str(sys.version_info[1]) + ") " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")

if not pass_check:
	sys.exit("CRITICAL ERROR: The python version found (version " +
	         str(sys.version_info[0]) + "." + str(sys.version_info[1]) + ") " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")


try:
	import setuptools
except ImportError:
	sys.exit("Please install setuptools.")

# check setuptools version
required_setuptools_version_major = 1
try:
	setuptools_version = setuptools.__version__
	setuptools_version_major = int(setuptools_version.split(".")[0])
	if setuptools_version_major < required_setuptools_version_major:
		sys.exit("CRITICAL ERROR: The setuptools version found (version " +
		         setuptools_version + ") does not match the version required " +
		         "(version " + str(required_setuptools_version_major) + "+)." +
					" Please upgrade your setuptools version.")
except (ValueError, IndexError, NameError):
	sys.exit("CRITICAL ERROR: Unable to call setuptools version. Please upgrade setuptools.")

from setuptools.command.install import install as _install

import distutils

# try to import urllib.request.urlretrieve for python3
try:
	from urllib.request import urlretrieve
except ImportError:
	from urllib import urlretrieve

def byte_to_megabyte(byte):
	"""
	Convert byte value to megabyte
	"""

	return byte / (1024.0 ** 2)


class ReportHook():
	def __init__(self):
		self.start_time = time.time()

	def report(self, blocknum, block_size, total_size):
		"""
		Print download progress message
		"""

		if blocknum == 0:
			self.start_time = time.time()
			if total_size > 0:
				print("Downloading file of size: " + "{:.2f}".format(byte_to_megabyte(total_size)) + " MB\n")
		else:
			total_downloaded = blocknum * block_size
			status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

			if total_size > 0:
				percent_downloaded = total_downloaded * 100.0 / total_size
				# use carriage return plus sys.stdout to overwrite stdout
				try:
					download_rate = total_downloaded / (time.time() - self.start_time)
					estimated_time = (total_size - total_downloaded) / download_rate
				except ZeroDivisionError:
					download_rate = 0
					estimated_time = 0
				estimated_minutes = int(estimated_time / 60.0)
				estimated_seconds = estimated_time - estimated_minutes * 60.0
				status += "{:3.2f}".format(percent_downloaded) + " %  " + \
				          "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
				          "{:2.0f}".format(estimated_minutes) + " min " + \
				          "{:2.0f}".format(estimated_seconds) + " sec "
			status += "        \r"
			sys.stdout.write(status)


def download(url, download_file):
	"""
	Download a file from a url
	"""

	try:
		print("Downloading " + url)
		file, headers = urlretrieve(url, download_file, reporthook=ReportHook().report)
		# print final return to start new line of stdout
		print("\n")
	except EnvironmentError:
		print("WARNING: Unable to download " + url)


def download_unpack_tar(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	download(url, download_file)

	error_during_extract = False

	try:
		tarfile_handle = tarfile.open(download_file)
		tarfile_handle.extractall(path=folder)
		tarfile_handle.close()
	except (EnvironmentError, tarfile.ReadError):
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)


def download_unpack_zip(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	download(url, download_file)

	error_during_extract = False

	try:
		zipfile_handle = zipfile.ZipFile(download_file)
		zipfile_handle.extractall(path=folder)
		zipfile_handle.close()
	except EnvironmentError:
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)

def install_databases(final_install_folder, mac_os, replace_install=None):
	url="https://huttenhower.sph.harvard.edu/baqlava-db/databases.V0.5.tar.gz"
	download_file_name="databases.V0.5.tar.gz"
	database_full_path=tempfile.mkdtemp(prefix="baqlava-db-download_",dir=setup_directory)
	print("Downloading baqlava protein and nucleotide databases.")
	error_during_install=False
	download_unpack_tar(url, download_file_name, database_full_path, "database")   
    # copy the installed software to the final bin location
	try:
		shutil.copytree( os.path.join(database_full_path), os.path.join(setup_directory,"baqlava"),dirs_exist_ok=True)
	except (EnvironmentError, shutil.Error):
		error_during_install=True        
	try:
		shutil.rmtree(database_full_path)
	except EnvironmentError:
		print("Warning: Unable to remove temp install folder.")
        
	if error_during_install:
		print("Warning: Unable to install baqlava databases.")
	else:
		print("Installed baqlava database at "+setup_directory)
        

class Install(_install):
    """
    Custom setuptools install command
    """
    
    _install.user_options=_install.user_options+[('bypass-dependencies-install', 
        None, 'bypass install of dependencies')]
    
    def initialize_options(self):
        self.bypass_dependencies_install=False
        _install.initialize_options(self)
    
    def finalize_options(self):
        _install.finalize_options(self)
    
    def run(self):
        # find out the platform
        mac_os=False
        if sys.platform in ["darwin","os2","os2emx"]:
            mac_os=True
        
        # install dependencies if not already installed
        if not self.bypass_dependencies_install:
            # install_databases(self.install_lib,mac_os,replace_install=False)
            database=False
        else:
            print("Bypass install of databases.")
		
        _install.run(self)
        print("----------------------------------------------------------")
        print("Install BAQLaVa database here: "+ os.path.join(self.install_lib,"baqlava"))
        print("USE COMMAND: baqlava_database --download database baqlava-db "+ os.path.join(self.install_lib,"baqlava"))
        print("----------------------------------------------------------")

setuptools.setup(
	name="baqlava",
	author=AUTHOR,
	author_email=MAINTAINER_EMAIL,
	version=VERSION,
	license="MIT",
	description="BAQLaVa: Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa).",
	url="https://github.com/biobakery/baqlava",
	keywords=['microbial', 'microbiome', 'bioinformatics', 'microbiology', 'metagenomic', 'baqlava' ,'anadama2'],
	platforms=['Linux', 'MacOS'],
	classifiers=[
		"Programming Language :: Python",
		"Development Status :: 5 - Production/Stable",
		"Environment :: Console",
		"Operating System :: MacOS",
		"Operating System :: Unix",
		"Programming Language :: Python :: 3.10",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
	#install_requires=['anadama2>=0.7.4'],
	packages=setuptools.find_packages(),
	cmdclass={'install': Install},
	entry_points={
		'console_scripts': [
			'baqlava = baqlava.baqlava:main',
                        'baqlava_join_tables = baqlava.utility_scripts.join_tables:main',
						"baqlava_database = baqlava.download_db:main",
		]},
	package_data={
		'baqlava': [
			'configs/*',
			'utility_files/*',
			'utility_scripts/*',
			'data/BAQLaVa.V0.5.nucleotide/*',
			'data/BAQLaVa.V0.5.protein/*',
		]},
		include_package_data=True, 
	zip_safe=False
)
