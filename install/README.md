# BEMRosetta Install

## Windows
No install is required. Just open BEMRosetta.exe and dodge Windows warnings. 

## U++ POSIX/X11 installation

Standard POSIX/X11 distribution of BEMRosetta comes as a source tarball.

Before compiling source code, you must install a few developpement packages. Many POSIX/X11 distributions provides developpement packages with the same names. Sometimes, development package names do not match. You will have to find the corresponding names for your distribution.


### Debian/apt-get based distributions

Build requires: g++  make  libgtk2.0-dev  libnotify-dev  libbz2-dev  sox libgtkglext1-dev

How to install them:

if sudo is available and enabled on your distribution, copy/paste this in a terminal:
```
sudo apt-get install  g++  make  libgtk2.0-dev  libnotify-dev  libbz2-dev  sox libgtkglext1-dev
```

if sudo is not available:
```
su -c 'apt-get install  g++  make  libgtk2.0-dev  libnotify-dev  libbz2-dev  sox libgtkglext1-dev'
```


### Compile source code

First, uncompress source tarball and change dir to the new created directory.

Example: for bemrosetta_linux.tar.tar.gz
```
tar zxvf bemrosetta_linux.tar.tar.gz
cd bemrosetta_linux.tar
```

Use 'make' to compile and generate BEMRosetta.out:
```
make
```

Now you can start playing with BEMRosetta by invoking ./BEMRosetta/BEMRosetta.out
You might want to put BEMRosetta elsewhere later, e.g. inside ~/bin/ for example

