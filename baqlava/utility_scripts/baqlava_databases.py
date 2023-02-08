"""
HUMAnN: humann_databases module
Download databases an update config settings
Dependencies: None
To Run: humann_databases --download <database> <build> <install_location>
Copyright (c) 2014 Harvard School of Public Health
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
    
# Try to load one of the humann src modules to check the installation
try:
    from .. import check
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN python package." +
        " Please check your install.") 

# Check the python version
check.python_version()

import argparse
import os

from .. import config
from .. import utilities

# the locations of the current databases to download
current_downloads = [ config.get('hosted_databases','nucleotide_db_full'), 
                     config.get('hosted_databases','uniref_db_full')]

def download_databases(location):
    """
    Download and decompress the selected database
    """
    
    # download the database
    for i in current_downloads:
        os.system("wget -P " + location + " " + i)
        os.system("tar -zxvf " + os.path.join(location,i.split('/')[-1]))
    return None

download_databases(sys.argv[1])
