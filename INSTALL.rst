Installation
================================================================================

The simplest way to install the latest stable release of Change-O is via pip::

    > pip3 install changeo --user

The current development build can be installed using pip and git in similar fashion::

    > pip3 install git+https://github.com/immcantation/changeo@master --user

If you currently have a development version installed, then you will likely
need to add the arguments ``--upgrade --no-deps --force-reinstall`` to the
pip3 command.

Requirements
--------------------------------------------------------------------------------

The minimum dependencies for installation are:

+ `Python 3.7.0 <http://python.org>`__
+ `setuptools 2.0 <http://bitbucket.org/pypa/setuptools>`__
+ `NumPy 1.8 <http://numpy.org>`__
+ `SciPy 0.14 <http://scipy.org>`__
+ `pandas 0.24 <http://pandas.pydata.org>`__
+ `Biopython 1.81 <http://biopython.org>`__
+ `presto 0.7.0 <http://presto.readthedocs.io>`__
+ `airr 1.3.1 <https://docs.airr-community.org>`__
+ `PyYAML 5.1 <http://pyyaml.org>`__
+ `packaging 23.2 <https://packaging.pypa.io>`__

Some tools wrap external applications that are not required for installation.
Those tools require minimum versions of:

+ AlignRecords requires `MUSCLE 3.8 <http://www.drive5.com/muscle>`__
+ ConvertDb-genbank requires `tbl2asn <https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2>`__
+ AssignGenes requires `IgBLAST 1.6 <https://ncbi.github.io/igblast>`__, but
  version 1.11 or higher is strongly recommended.
+ BuildTrees requires `IgPhyML 1.0.5 <https://github.com/immcantation/igphyml>`_

Linux
--------------------------------------------------------------------------------

1. The simplest way to install all Python dependencies is to install the
   full SciPy stack using the
   `instructions <http://scipy.org/install.html>`__, then install
   Biopython according to its
   `instructions <http://biopython.org/DIST/docs/install/Installation.html>`__.

2. Install `presto 0.7.0 <http://presto.readthedocs.io>`__ or greater.

3. Download the `Change-O bundle <https://github.com/immcantation/changeo/tags>`__
   and run::

   > pip3 install changeo-x.y.z.tar.gz --user

Mac OS X
--------------------------------------------------------------------------------

1. Install Xcode. Available from the Apple store or
   `developer downloads <http://developer.apple.com/downloads>`__.

2. Older versions Mac OS X will require you to install XQuartz 2.7.5. Available
   from the `XQuartz project <http://xquartz.macosforge.org/landing>`__.

3. Install Homebrew following the installation and post-installation
   `instructions <http://brew.sh>`__.

4. Install Python 3.4.0+ and set the path to the python3 executable::

   > brew install python3
   > echo 'export PATH=/usr/local/bin:$PATH' >> ~/.profile

5. Exit and reopen the terminal application so the PATH setting takes effect.

6. You may, or may not, need to install gfortran (required for SciPy). Try
   without first, as this can take an hour to install and is not needed on
   newer releases. If you do need gfortran to install SciPy, you can install it
   using Homebrew::

   > brew install gfortran

   If the above fails run this instead::

   > brew install --env=std gfortran

7. Install NumPy, SciPy, pandas and Biopython using the Python package
   manager::

   > pip3 install numpy scipy pandas biopython

8. Install `presto 0.7.0 <http://presto.readthedocs.io>`__ or greater.

9. Download the `Change-O bundle <https://github.com/immcantation/changeo/tags>`__, 
   open a terminal window, change directories to the download folder, and run::

   > pip3 install changeo-x.y.z.tar.gz

Windows
--------------------------------------------------------------------------------

1. Install Python 3.4.0+ from `Python <http://python.org/downloads>`__,
   selecting both the options 'pip' and 'Add python.exe to Path'.

2. Install NumPy, SciPy, pandas and Biopython using the packages
   available from the
   `Unofficial Windows binary <http://www.lfd.uci.edu/~gohlke/pythonlibs>`__
   collection.

3. Install `presto 0.7.0 <http://presto.readthedocs.io>`__ or greater.

4. Download the `Change-O bundle <https://github.com/immcantation/changeo/tags>`__, 
   open a Command Prompt, change directories to the download folder, and run::

   > pip install changeo-x.y.z.tar.gz

5. For a default installation of Python 3.4, the Change-0 scripts will be
   installed into ``C:\Python34\Scripts`` and should be directly
   executable from the Command Prompt. If this is not the case, then
   follow step 6 below.

6. Add both the ``C:\Python34`` and ``C:\Python34\Scripts`` directories
   to your ``%Path%``. On both Windows 7 and Windows 10, the ``%Path%`` setting is located under Control Panel -> System and Security -> System -> Advanced System Settings -> Environment variables -> System variables -> Path.

7. If you have trouble with the ``.py`` file associations, try adding ``.PY``
   to your ``PATHEXT`` environment variable. Also, try opening a
   Command Prompt as Administrator and run::

    > assoc .py=Python.File
    > ftype Python.File="C:\Python34\python.exe" "%1" %*
