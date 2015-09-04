Version 0.3.0.beta:  September 4, 2015
==============================================================================

Change-O is a collection of tools for analyzing immunoglobulin sequences. 
Dramatic improvements in high-throughput sequencing technologies now enable 
large-scale characterization of immunoglobulin (Ig) repertoires, defined as the 
collection of trans-membrane antigen-receptor proteins located on the surface 
of T and B lymphocytes. Change-O is a suite of utilities to handle advanced 
analysis of Ig sequences following germline segment assignment. Change-O 
handles output from IMGT/HighV-QUEST and IgBLAST. It includes features for 
creating a personalized genotype, identifying sequences that are from a single 
B cell clone and inferring its lineage tree, analyzing amino acid properties, 
calculating diversity, generating a model of somatic hypermutation, and 
quantifying selection pressure. Record sorting, grouping, and sampling 
operations are also included.


Requirements
-------------------------------------------------------------------------------

Software        | Link
--------------- | -----------------------------
Python 3.4.0    | http://python.org
setuptools 2.0  | http://bitbucket.org/pypa/setuptools
NumPy 1.8       | http://numpy.org
SciPy 0.14      | http://scipy.org
pandas 0.15     | http://pandas.pydata.org
presto 0.5      | http://clip.med.yale.edu/presto
Biopython 1.65  | http://biopython.org
MUSCLE v3.8     | http://www.drive5.com/muscle


Installation - Linux
-------------------------------------------------------------------------------

1. The simplest way to install all Python dependencies is to install the full
   SciPy stack using the instructions at http://scipy.org/install.html, then
   install Biopython according to the
   [instructions](http://biopython.org/DIST/docs/install/Installation.html).
   
2. Install pRESTO 0.5.0+ from http://clip.med.yale.edu/presto.

3. Extract the pRESTO bundle and run `python3 setup.py install --user`.

Installation - Windows
-------------------------------------------------------------------------------

1. Install Python 3.4.0+ from [Python](http://python.org/downloads), selecting
   both the options 'pip' and 'Add python.exe to Path'.

2. Install NumPy, SciPy, pandas and Biopython using the packages available 
   from the [Unofficial Windows](http://www.lfd.uci.edu/~gohlke/pythonlibs)
   binary collection.

4. Unzip the Change-O bundle, open a Command Prompt, and run
   `python setup.py install` from the Change-O folder.
   
5. For a default installation of Python 3.4, the Change-O scripts will be 
   installed into `C:\Python34\Scripts` and should be directly executable from 
   the Command Prompt. If this is not the case, then follow steps 5-6 below.
   
6. Add both the 'C:\Python34' and 'C:\Python34\Scripts' directories to your 
   `%Path%`. On Windows 7 the %Path% setting is located under 'Control Panel' 
    -> 'System and Security' -> 'System' -> 'Advanced System Settings' -> 
    'Environment variables' -> 'System variables' -> 'Path'.

7. Set the file association for Python ('.py') files by right-clicking on a 
   '.py' file, selecting 'Open with' -> 'Choose default program...', choosing the 
   'python.exe' executable from the Python 3.4 folder, and checking 
   'Always use the selected program'.


Installation - Mac OS X
-------------------------------------------------------------------------------

1. Install Xcode 3.2.6
   Available from the Apple store or
   [developer downloads](http://developer.apple.com/downloads).
   If you have a newer version (eg, Xcode 4.6.3) that will work also,
   but Xcode 3 is free of charge.  If Xcode fails to install with an
   "Unknown Error", change the date on your system to some time in 2011,
   install Xcode, and then change the date back to the proper setting.

2. Install XQuartz 2.7.5
   Available from the [XQuartz project](http://xquartz.macosforge.org/landing).

3. Install Homebrew
   Follow the installation and post-installation [instructions](http://brew.sh).

4. Open a terminal and install gfortran (required for SciPy) using Homebrew
   (this can take an hour to install):  
   `> brew install gfortran`  
   If the above fails run this instead:  
   `> brew install --env=std gfortran`  

5. Install Python 3.4.0+ and set it the path to the python executable:  
   `> brew install python3`  
   `> echo 'export PATH=/usr/local/bin:$PATH' >> ~/.profile`  
   Exit and reopen the terminal application so the PATH setting takes effect

6. Install NumPy, SciPy, pandas and Biopyton using the Python package manager:  
   `> pip3 install numpy`  
   `> pip3 install scipy`  
   `> pip3 install pandas`  
   `> pip3 install biopython`  

7. Extract the Change-O bundle, open a terminal window, and run
   `python3 setup.py install` from the Change-O folder.
