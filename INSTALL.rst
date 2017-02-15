Installation
================================================================================

The simplest way to install the latest stable release of Change-O is via pip::

    > pip3 install changeo --user

The current development build can be installed using pip and mercurial in similar fashion::

    > pip3 install hg+https://bitbucket.org/kleinstein/changeo#default --user

If you currently have a development version installed, then you will likely
need to add the arguments ``--upgrade --no-deps --force-reinstall`` to the
pip3 command.

Requirements
--------------------------------------------------------------------------------

+  `Python 3.4.0 <http://python.org>`__
+  `setuptools 2.0 <http://bitbucket.org/pypa/setuptools>`__
+  `NumPy 1.8 <http://numpy.org>`__
+  `SciPy 0.14 <http://scipy.org>`__
+  `pandas 0.15 <http://pandas.pydata.org>`__
+  `Biopython 1.65 <http://biopython.org>`__
+  `presto 0.5 <http://presto.readthedocs.io>`__

Linux
--------------------------------------------------------------------------------

1. The simplest way to install all Python dependencies is to install the
   full SciPy stack using the
   `instructions <http://scipy.org/install.html>`__, then install
   Biopython according to its
   `instructions <http://biopython.org/DIST/docs/install/Installation.html>`__.

2. Install `presto 0.5.0 <http://presto.readthedocs.io>`__ or greater.

3. Download the Change-O bundle and run::

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

7. Install NumPy, SciPy, pandas and Biopyton using the Python package
   manager::

   > pip3 install numpy scipy pandas biopython

8. Install `presto 0.5.0 <http://presto.readthedocs.io>`__ or greater.

9. Download the Change-O bundle, open a terminal window, change directories
   to the download folder, and run::

   > pip3 install changeo-x.y.z.tar.gz

Windows
--------------------------------------------------------------------------------

1. Install Python 3.4.0+ from `Python <http://python.org/downloads>`__,
   selecting both the options 'pip' and 'Add python.exe to Path'.

2. Install NumPy, SciPy, pandas and Biopython using the packages
   available from the
   `Unofficial Windows binary <http://www.lfd.uci.edu/~gohlke/pythonlibs>`__
   collection.

3. Install `presto 0.5.0 <http://presto.readthedocs.io>`__ or greater.

4. Download the Change-O bundle, open a Command Prompt, change directories to
   the download folder, and run::

   > pip install changeo-x.y.z.tar.gz

5. For a default installation of Python 3.4, the Change-0 scripts will be
   installed into ``C:\Python34\Scripts`` and should be directly
   executable from the Command Prompt. If this is not the case, then
   follow step 5 below.

6. Add both the ``C:\Python34`` and ``C:\Python34\Scripts`` directories
   to your ``%Path%``. On Windows 7 the ``%Path%`` setting is located
   under Control Panel -> System and Security -> System -> Advanced
   System Settings -> Environment variables -> System variables -> Path.

6. If you have trouble with the ``.py`` file associations, try adding ``.PY``
   to your ``PATHEXT`` environment variable. Also, opening a
   command prompt as Administrator and run::

    > assoc .py=Python.File
    > ftype Python.File="C:\Python34\python.exe" "%1" %*
