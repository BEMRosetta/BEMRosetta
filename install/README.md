# BEMRosetta Install

## Binaries

### Windows
No install is required. Just copy bin/BEMRosetta.exe and bin/BEMRosetta_cl.exe anywhere in your computer and dodge Windows warnings. 

## Compiling

BEMRosetta is proudly developed with the [U++ framework](https://www.ultimatepp.org/) using the C++ language. U++ makes C++ so simple and efficient that it appeals even to scripting language programmers.

U++ is simple to install. Uninstalling it is just a matter of deleting the install folder. 
1. Download it from the [U++ nightly builds](https://www.ultimatepp.org/www$uppweb$download$en-us.html), choosing the latest U++ for Windows (with CLANG) or for Linux/BSD/Solaris. It includes some additional tools and precompiled binaries.
2. Download or clone it in the same folder from [Github](https://github.com/ultimatepp/ultimatepp) frequently to get the most updated sources.
3. Start theide to finish the install
4. Install [Anboto libraries](https://github.com/anboto/Anboto) in the same folder as U++.
5. Clone or copy BEMRosetta to your computer. The best place to copy BEMRosetta is the same folder as U++. 

A normal folder structure could be:
* U++ParentFolder/upp
* U++ParentFolder/Anboto
* U++ParentFolder/BEMRosetta

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Download.png" width="600" title="U++ download"></p>

* Go to U++ folder (normally 'upp') and open uppsrc.var with a text editor like Notepad.
  * In Windows, it will be in 'upp' folder
  * In Linux, it will be in 'upp/.config/u++/theide' folder 
* uppsrc.var will contain something like:
<pre>
UPP    = "U++ParentFolder/upp/uppsrc";
OUTPUT = "U++ParentFolder/upp/out";
</pre>
* Replace the content of "UPP = ..." line with:
<pre>
UPP    = "U++ParentFolder/BEMRosetta;U++ParentFolder/Anboto;U++ParentFolder/upp/uppsrc;U++ParentFolder/upp/bazaar";
</pre>
* Save the file as 'BEMRosetta.var'

From now on you will have 'BEMRosetta' ready to be selected in the 'Select main package' dialog of TheIDE, and called from umk, the command line U++ make tool.

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
* Unpack the archive
* Open a terminal and cd to the upp folder
* Run the ./install script. It should detect your Linux distro and suggest the command to install the required dependencies. It will also install them, or else you can copy the command to another terminal and run it yourself.
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Install.png" width="600" title="U++ install"></p>
&nbsp;&nbsp;&nbsp;&nbsp;In case of any problem the normal dependencies to install are g++, make, libgtk2.0-dev, libnotify-dev, libbz2-dev, sox, libgtkglext1-dev. To install them:

  * if sudo is available and enabled on your distribution, copy/paste this in a terminal:
```
sudo apt-get install  g++  make  libgtk2.0-dev  libnotify-dev  libbz2-dev  sox libgtkglext1-dev
```

  * if sudo is not available:
```
su -c 'apt-get install  g++  make  libgtk2.0-dev  libnotify-dev  libbz2-dev  sox libgtkglext1-dev'
```

* Compile BEMRosetta with umk.exe. umk (U++ make) is a command line utility to build U++ programs:
```
umk BEMRosetta BEMRosetta    CLANG -r +GUI,SHARED ~/bemrosetta.exe
umk BEMRosetta BEMRosetta_cl CLANG -r +SHARED     ~/bemrosetta_cl.exe
```
Put BEMRosetta binaries where you want them, e.g. inside ~/bin . No install is required.
