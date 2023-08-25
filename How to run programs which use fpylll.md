Information collected by Vivian. This is also largely for me in the future, in case I need to do something like this again.
# Introduction
[fpylll](https://github.com/fplll/fpylll) is a Python wrapper for the C library [fplll](https://github.com/fplll/fplll): it provides an interface in Python to use the fplll library.

Essentially, it uses another library, Cython, to be able to call C code from Python. Thus the numerical algorithms are all implemented in C, and the fplll library is a prerequisite for fpylll.
# Install fpylll
## Note on Operating Systems
I believe fpylll should be able to run on both Linux and Mac. However, all three of us use Windows, so we were not able to confirm whether Mac works or not.
## Try first (Mac or Linux)
If you're on Mac or Linux, you can first try installing using PyPI, through the command line:
```
pip install fpylll
```

It's probably worth trying to make this option work, updating/installing pip if necessary. This option will not work on Windows (I think).
## Linux (Ubuntu/Debian)
If installing from pip didn't work, fpylll can be installed from the terminal using Advanced Package Tool (APT).

The code to install fpylll is:
```
sudo apt-get install python3-fpylll
```

This is how us three got fpylll working on our laptops, but I think pip would have been easier.
## Mac
I'm not sure of another simple way to install fpylll. You can try downloading the code from the [fpylll GitHub](https://github.com/fplll/fpylll) and running `bootstrap.sh` using the terminal (this is from the GitHub):
```sh
$ ./bootstrap.sh
$ source ./activate
```
## Windows
One option is to install "Windows Subsystem for Linux" (WSL), to avoid using a virtual machine. Once WSL is installed (see a [tutorial](https://learn.microsoft.com/en-us/windows/wsl/install) for more details), you can connect to WSL through the terminal and install fpylll onto your new Linux subsystem.

To manage our code in WSL, we used Visual Studio Code, which I would recommend.
# Run the code
You should now be able to run fpylll code. If it still doesn't work, your Python might not be able to "see" the fpylll library - it may have installed to the wrong version of Python, if you have multiple on your computer.