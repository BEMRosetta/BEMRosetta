# BEMRosetta Install

## Binaries

### Windows
No compiling is required. Just go to the [last release](https://github.com/BEMRosetta/BEMRosetta/releases), copy bin/BEMRosetta.exe and bin/BEMRosetta_cl.exe anywhere in your computer and dodge Windows warnings. 

## Compiling

BEMRosetta is proudly developed with the [U++ framework](https://www.ultimatepp.org/) using the C++ language. 
U++ makes C++ so simple and efficient that it appeals even to scripting language programmers.

U++ is simple to [install](https://www.ultimatepp.org/www$uppweb$download$en-us.html). Uninstalling it is just a matter of deleting the install folder. 

To compile BEMRosetta, open theide (the U++ IDE), choose in the right the "MyApps" assembly, and tap on UppHub:

![Untitled](https://github.com/BEMRosetta/BEMRosetta/assets/38589221/0d0fceba-9bb7-4311-88db-b568a1c07ad7)

In the left grid, choose BEMRosetta and click Install:

![Untitled2](https://github.com/BEMRosetta/BEMRosetta/assets/38589221/e7fe4d46-4fb5-4d2b-8d37-13eb8274acb9)

UppHub will use Git to install all necessary packages. Wait until "Done" appears, and tap on Close ...

![Untitled3](https://github.com/BEMRosetta/BEMRosetta/assets/38589221/2bf64a0d-5e64-4512-b9c5-b03e461f7949)

and next in Exit:

![Untitled4](https://github.com/BEMRosetta/BEMRosetta/assets/38589221/0297e39d-73de-41a5-8baf-93512e63a6d7)

The sources are downloaded. To compile them, choose MyApps in the in the right, and BEMROsetta in the lower list:

![Untitled5](https://github.com/BEMRosetta/BEMRosetta/assets/38589221/f0365218-5eb0-40a4-abaa-e504ff2d2819)

Then BEMRosetta and BEMRosetta_cl will appear. Choose which one do you want to compile.

![Untitled6](https://github.com/BEMRosetta/BEMRosetta/assets/38589221/d4f487eb-60d3-42f2-b7d1-3af787138fd7)

### Windows

In Windows, U++ binaries are precompiled and U++ comes with all dependencies needed, including [CLANG](https://clang.llvm.org/) compiler. Simply:
* Open a terminal, cd to the upp folder
* Compile BEMRosetta with umk.exe. the command line U++ make tool:
```
umk BEMRosetta BEMRosetta    CLANGX64 -r +GUI .\BEMRosetta.exe
umk BEMRosetta BEMRosetta_cl CLANGX64 -r      .\BEMRosetta_cl.exe
```
Put BEMRosetta binaries where you want them. No install is required.

### Linux POSIX 

In Linux, U++ has to be compiled so that the dependencies of the target distribution are complied. However the process is quick and mostly automated for majority of popular Linux flavors. Simply:
* Open a terminal and cd to the upp folder
* Run the ./install script. It should detect your Linux distro and suggest the command to install the required dependencies. It will also install them, or else you can copy the command to another terminal and run it yourself.
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Install.png" width="600" title="U++ install"></p>
&nbsp;&nbsp;&nbsp;&nbsp;In case of any problem the normal dependencies to install are g++, make, libgtk2.0-dev, libnotify-dev, libbz2-dev, sox, libgtkglext1-dev. To install them:

&nbsp;&nbsp;&nbsp;&nbsp;- if sudo is available and enabled on your distribution, copy/paste this in a terminal:
```
sudo apt-get install  g++  make  libgtk2.0-dev  libnotify-dev  libbz2-dev  sox libgtkglext1-dev
```

&nbsp;&nbsp;&nbsp;&nbsp;- if sudo is not available:
```
su -c 'apt-get install  g++  make  libgtk2.0-dev  libnotify-dev  libbz2-dev  sox libgtkglext1-dev'
```

* Compile BEMRosetta with umk.exe. umk (U++ make) is a command line utility to build U++ programs:
```
umk BEMRosetta BEMRosetta    CLANG -r +GUI,SHARED ~/bemrosetta.exe
umk BEMRosetta BEMRosetta_cl CLANG -r +SHARED     ~/bemrosetta_cl.exe
```
Put BEMRosetta binaries where you want them, e.g. inside ~/bin . No install is required.
