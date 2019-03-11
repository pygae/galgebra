Installation on Linux, Windows, and Mac
=======================================

Install python
--------------

The `galgebra` python module, which is an implementation of geometric algebra in python has two perquisites for a minimal installation, python and `sympy`. For the python language we have the following situation[1].

| os | python installation |
|:-|:-|
| linux | Comes with all versions of linux |
| windows | To install python on windows go to https://www.python.org/downloads/windows/ and install version appropriate for you version of windows. If you wish a more complete/advanced installation go to https://code.google.com/p/pythonxy/. |
| mac | Basic version comes with OSX. For better in-stallation go to http://docs.python-guide.org/en/ latest/starting/install/osx/. |

Install sympy
-------------

For `sympy` there are two alternatives for installation.

| mode | method |
|:-|:-|
| latest release | Go to https://github.com/sympy/sympy/releases and select option appropriate for your system. Note that if you have pip (see https://pip.pypa.io/en/latest/installing.html) installed you can install the latest re-lease by entering the command `pip install sympy`. |
| development version | The method for the development version is preferred since that method always builds `sympy` with the python system you have installed on your system (32-bits verses 64-bits and particular python version you are running). |

Install galgebra
----------------

Since you are reading this document you have already obtained a copy of `galgebra`. If you wish to obtain the very latest version (assuming you have not already done this) go to <https://github.com/brombo/galgebra> and download and extract the zipped archive.

Then with whatever version you are using open a terminal/command line in the `galgebra` directory that is in the top directory of the archive. If you are in the correct the directory it should contain the python program `setgapth.py`. If you are in linux or osx run the program with the command `sudo python setgapth.py`, if in windows use `python setgapth.py`.

This program creates the file `Ga.pth` in the correct directory to simplify importing the `galgebra` modules into your python programs. The modules you will use for programming with geometric algebra/calculus are `ga`, `mv`, `lt`, and `printer`[2]. To import any of these modules into your program, say `mv`, you only have to enter in the program `import mv`. It does not matter where the program file is located.

LaTeX Options
--------------------------

In order to use the latex output of the `galgebra` modules (excluding latex output from *Ipython notebook*) you must install a latex distribution. Directions follow if you do not already have LaTeX installed.[3]

| os | latex installation |
|:-|:-|
| linux | Open a terminal and run `sudo apt-get install texlive-full`. It takes about half an hour to install. |
| windows | Go to http://miktex.org/download (other downloads). Download a net installer. Install a full version of MikTex. |
| mac | Go to http://www.tug.org/mactex/ and follow instructions to install MacTeX. |

“Ipython notebook” Options
--------------------------

To use *ipython notebook* with `galgebra` it must be installed. To install *ipython notebook* do the following.

Google “get-pip.py” and click on the first entry “get-pip.py”. Then follow the instructions to download “get-pip.py”. Open a terminal/command line in the directory of the download and execute `python get-pip.py` for windows or `sudo python get-pip.py` for linux. The reason for install *pip* in this manner is that it insures the correct settings for the version of python you are using. Then run in a terminal/command line `pip install ipython[notebook]`. If you have already installed *ipython notebook* you should enter `pip install ipython[notebook] –upgrade` to make sure you have the latest version. Linux and OSX users will have to use `sudo` with the commands. The version of *ipython notebook* we are using is **jupyter** and that should be shown when the notebook is started.

Note that to correctly print latex from *ipython notebook* one must use the `Format()` function from the *printer* module. Go to the section on latex printing for more information.

The ANSI Console
----------------

The `printer` module of `galgebra` contains the class `Eprint` which is described in section ([stdprint]). This function uses the capabilities of the ansi console (terminal) for enhanced multivector printing where multivector bases, sympy functions and derivatives are printed in different colors. The ansi console is native to Linux and OSX (which is really Unix under the hood), but not windows. The best available free substitute for the ansi console on windows is *ConEmu*. The web page for ConEmu is <http://conemu.github.io/>. In order to install *ConEmu* download the appropriate version of the *ConEmu* installer (exe file) for your system (32 bit or 64 bit) from the website and and execute it. Instructions for using *ConEmu* are given in section ([stdprint]).

Geany Programmers Editor
------------------------

*Geany* is a very nice *free* programmers editor that work well with *python*. From within *geany* you can execute a *python* program. The *galgebra* printing system is setup so that you can display the program output on an ansi terminal or if you are using the LaTeXoptions has the terminal launch a *pdf* browser to view the LaTeXoutput. To install *geany* on Linux use the command line `sudo apt-get install geany`, on Windows go to <http://www.geany.org/Download/Releases> or to install *geany* in OSX go to <http://wiki.geany.org/howtos/osx/running>.

[1] Currently `galgebra` supports python versions 2.7+, but not versions 3.0+ of python.

[2] All these modules are in the same directory as `setgapth.py`.

[3] In order for `galgebra` to output latex formatted pdf files your distribution of latex must have `pdflatex` installed.