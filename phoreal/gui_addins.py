# -*- coding: utf-8 -*-
"""
Simple script to load non-Python files to help bundle EXE

Copyright 2019 Applied Research Laboratories, University of Texas at Austin

This package is free software; the copyright holder gives unlimited
permission to copy and/or distribute, with or without modification, as
long as this notice is preserved.

Authors:
    Mike Alonzo
    Eric Guenther
    
Date: September 20, 2019
"""

# Import modules
import os

root = os.path.dirname(__file__)

superFilterFile_windows = os.path.join(root,'closest_x64.dll')

superFilterFile_linux = os.path.join(root,'phorealc.so')

viewerBlank_html = os.path.join(root,'Viewer_blank.html')

viewerBlankOnline_html = os.path.join(root,'Viewer_Online_blank.html')