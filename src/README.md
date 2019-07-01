# BEMRosetta Install

## Windows
No install is required. Just copy bin/BEMRosetta.exe and bin/BEMRosetta_cl.exe anywhere and dodge Windows warnings. 

## Linux POSIX 

Standard POSIX distribution of BEMRosetta comes as a source tarball.

Before compiling source code, you must install a few developpement packages. Many POSIX distributions provides developpement packages with the same names. However, sometimes development package names do not match. In this case, you will have to find the corresponding names for your distribution.


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

Change dir to the new created directory.

```
cd bemrosetta_files
```

Open a terminal and run 'make_gui.sh' to compile and generate BEMRosetta, and 'make_cli.sh' to compile and generate BEMRosetta_cl:

You might want to put BEMRosetta binaries elsewhere, e.g. inside ~/bin/ for example.
