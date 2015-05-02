Linux
-------------------------------------------------------------------------------

1. The simplest way to install all Python dependencies is to install the full 
   SciPy stack using the instructions at http://scipy.org/install.html, then 
   install Biopython according to the 
   [instructions](http://biopython.org/DIST/docs/install/Installation.html).

2. Unzip the Change-O bundle into a directory of your choice and add that 
   directory to your `$PATH`.


Windows
-------------------------------------------------------------------------------

1. Install Python 2.7.5+ from [Python](http://python.org/download).

2. Install NumPy, SciPy, pandas and Biopython using the packages available 
   from the [Unofficial Windows](http://www.lfd.uci.edu/~gohlke/pythonlibs)
   binary collection.

3. Unzip the Change-O bundle into a directory of your choice and add that 
   directory to your `%Path%`.  On Windows 7 the `%Path%` setting is located under
   'Control Panel' -> 'System and Security' -> 'System' -> 
   'Advanced System Settings' -> 'Environment variables' -> 'System variables' 
   -> 'Path'

4. The Change-O scripts should then be directly executable from the Command Prompt.
   If not, correct the file association for `.py` files by right-clicking on a 
   `.py` file, selecting 'Open with' -> 'Choose default program...', choosing the 
   `python.exe` from Python 2.7, and checking 'Always use the selected program'.
   

Mac OS X
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
   
5. Install Python 2.7.5+ and set it as the default python executable:  
   ```
   > brew install python
   > echo 'export PATH=/usr/local/bin:$PATH' >> ~/.profile
   ```
   Exit and reopen the terminal application so the PATH setting takes effect

6. Install NumPy, SciPy, pandas and Biopyton using the Python package manager:  
   ```
   > pip install numpy
   > pip install scipy
   > pip install pandas
   > pip install biopython
   ```
   
7. Add the Change-O installation to your `PATH` setting. For example,
   if you copy the Change-O scripts into `/Users/Username/changeo`, then set:  
   `> echo 'export PATH=$HOME/changeo:$PATH' >> ~/.profile`  
   Exit and reopen the terminal application so the PATH setting takes effect
