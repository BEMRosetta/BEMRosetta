# BEMRosetta Install

## Binaries

### Windows
Just copy bin/BEMRosetta.exe or bin/BEMRosetta_cl.exe anywhere in your computer and dodge Windows warnings. 

## Compiling

BEMRosetta is proudly developed with the [U++ framework](https://www.ultimatepp.org/) using the C++ language. U++ makes C++ so simple and efficient that it appeals even to python users.

U++ is simple to install. Uninstalling it is just a matter of deleting the install folder. You can download it from the [U++ nightly builds](https://www.ultimatepp.org/www$uppweb$download$en-us.html), choosing the latest U++ for Windows (with CLANG) or for Linux/BSD/Solaris. You can find a repository in [Github](https://github.com/ultimatepp/ultimatepp) too.

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Download.png" width="600" title="U++ download"></p>

### Windows

In Windows, U++ binaries are precompiled and U++ comes with all dependencies needed, including CLANG compiler. Simply:
* Unpack the archive
* Start theide.exe to finish the install
* Copy src/BEMRosetta and src/BEMRosetta_cl folders to upp/bazaar folder
* Open a terminal, cd to upp folder, and compile BEMRosetta:
```
umk BEMRosetta BEMRosetta    CLANGX64 -r +GUI .\BEMRosetta.exe
umk BEMRosetta BEMRosetta_cl CLANGX64 -r      .\BEMRosetta_cl.exe
```
Put BEMRosetta binaries where you want them. No install required.

### Linux POSIX 

Standard POSIX distribution of BEMRosetta comes as a source tarball.

Before compiling source code, you must install a few developpement packages. Many POSIX distributions provides developpement packages with the same names. However, sometimes development package names do not match. In this case, you will have to find the corresponding names for your distribution.


### Debian/apt-get based distributions

Build requires: g++  make  libgtk2.0-dev  libnotify-dev  libbz2-dev  sox  libgtkglext1-dev  libxtst-dev

How to install them:

if sudo is available and enabled on your distribution, copy/paste this in a terminal:
```
sudo apt-get install g++ make libgtk2.0-dev libnotify-dev libbz2-dev sox libgtkglext1-dev libxtst-dev
```

if sudo is not available:
```
su -c 'apt-get install g++ make libgtk2.0-dev libnotify-dev libbz2-dev sox libgtkglext1-dev libxtst-dev'
```


### Compile source code

First, uncompress source tarball bemrosetta_files.tar.gz and change dir to the new created directory.

```
tar zxvf bemrosetta_files.tar.gz
cd bemrosetta_files
```

Open a terminal and run 'make_gui.sh' to compile and generate BEMRosetta, and 'make_cli.sh' to compile and generate BEMRosetta_cl:

You might want to put BEMRosetta binaries elsewhere, e.g. inside ~/bin/ for example.
