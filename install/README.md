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
* Open a terminal, cd to the upp folder, and compile BEMRosetta with umk.exe. umk (U++ make) is a command line utility to build U++ programs:
```
umk BEMRosetta BEMRosetta    CLANGX64 -r +GUI .\BEMRosetta.exe
umk BEMRosetta BEMRosetta_cl CLANGX64 -r      .\BEMRosetta_cl.exe
```
Put BEMRosetta binaries where you want them. No install required.

### Linux POSIX 

In Linux, U++ has to be compiled so that the dependencies of the target distribution are complied. However the process is quick and mostly automated for majority of popular Linux flavors. Simply:
* Open a terminal
* Unpack the archive
* cd to the upp folder
* Run the ./install script:
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Install.png" width="600" title="U++ download"></p>
The script should detect your Linux distribution and suggest the command to install required dependencies. It can even install these if you do not mind sharing your sudo password, or you can just copy the command to separate terminal and run it yourself.
* Compile BEMRosetta with umk.exe. umk (U++ make) is a command line utility to build U++ programs:
```
umk BEMRosetta BEMRosetta    CLANG -r +GUI,SHARED ~/bemrosetta.exe
umk BEMRosetta BEMRosetta_cl CLANG -r +SHARED     ~/bemrosetta_cl.exe
```
Put BEMRosetta binaries where you want them, e.g. inside ~/bin/ . No install is required
