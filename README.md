[![License](https://img.shields.io/badge/license-GPL(>=2)-C11B17.svg)](http://www.gnu.org/licenses/gpl-2.0.html)
[![Library](https://img.shields.io/badge/library-GSL(>=2.3)-C11B17.svg)](https://www.gnu.org/software/gsl/)

# rtmpt
*R* package for fitting Response-Time extended Multinomial Processing Tree models by Klauer and Kellen (2018)

## Installation

### CRAN
Unfortunately, an installation via [CRAN](https://cran.r-project.org/) is not yet possible, but will become possible in the future.

### Windows
The installation from source requires the installation of [GSL](https://www.gnu.org/software/gsl/) (>=2.3) and is rather complicated.

An easier way is to install the `rtmpt`package via binaries:
1. Download the binaries for your corresponding R version (choose `R(3.6.X)` if you have an R subversion of 3.6)
* [![Binary](https://img.shields.io/badge/binaries-R(3.5.X)-blue.svg)](https://github.com/RaphaelHartmann/rtmpt-files/blob/master/binaries/3.5/rtmpt_0.1-14.zip) - R old release
* [![Binary](https://img.shields.io/badge/binaries-R(3.6.X)-blue.svg)](https://github.com/RaphaelHartmann/rtmpt-files/blob/master/binaries/3.6/rtmpt_0.1-14.zip) - R release
* [![Binary](https://img.shields.io/badge/binaries-R(3.7.X)-blue.svg)](https://github.com/RaphaelHartmann/rtmpt-files/blob/master/binaries/3.7/rtmpt_0.1-14.zip) - R devel
2. Open Rstudio and install via *Tools -> Install Packages... ->* select *Install from: Package Archive File (.zip)* and *browse* the .zip file *->* click *Install*

### Linux
Before installing the *R* package `rtmpt` you need to install the latest version of the GNU Scientific Library (GSL):
1. Download the latest `*tar.gz` file of [GSL](http://ftpmirror.gnu.org/gsl/)
2. Unzip the file
3. Assuming you downloaded version 2.5 and unziped the file to ~/Downloads/gsl-2.5, open terminal and change directory via ```$ cd ~/Downloads/gsl-2.5```
4. Execute the following lines:
    ```
    ./configure
    make
    sudo make install
    ```
For some Linux distributions it is necessary to run `sudo ldconfig -v` in the terminal before installing `rtmpt`.

Installing `rtmpt` is possible either 

- with the *R* package [devtools](https://cran.r-project.org/web/packages/devtools/index.html) and the command `devtools::install_github("RaphaelHartmann/rtmpt@master")` in *R* **or** 
- with the following steps:
  
    1. Download all the source files from the `rtmpt` package above
  
    2. In the terminal change directory via `$ cd ~/Downloads` (assuming the *rtmpt* source folder is in /Downloads)
  
    3. Execute:
    ```
    R CMD build rtmpt
    R CMD INSTALL rtmpt_<XYZ>.tar.gz
    ```
    where XYZ stands for the version number.

### MacOS
Short instruction for installing rtmpt on MacOS:
1. If not already installed, install [Homebrew](https://brew.sh/)
2. Install GSL via Homebrew using terminal and command `brew install gsl`
3. Install Xcode from the App Store and then install the Command Line Tools using `xcode-select --install`
4. If you use MacOS Mojave, make sure to also run the following command:
```
open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg
```
5. Also execute `brew install gcc` in your terminal.
6. Check via `ls /usr/local/bin/` what endings the gcc/g++ executables have (in my case -9)
7. If `~/.R/Makevars` does not already exist, create it with `mkdir ~/.R/` and `touch ~/.R/Makevars`
8. Use `open ~/.R/Makevars` to write the following in the Makevars file (assuming your gcc/g++ executables have the ending -9):
    ```
    CC=/usr/local/bin/gcc-9
    CXX=/usr/local/bin/g++-9
    CXX11=/usr/local/bin/g++-9
    ```
    and save.

Installing the `rtmpt` package:

This is possible either

- with the *R* package [devtools](https://cran.r-project.org/web/packages/devtools/index.html) and the command `devtools::install_github("RaphaelHartmann/rtmpt@master")` in *R* **or** 
- with the following steps:
  
    1. Download all the source files from the `rtmpt` package above
  
    2. In the terminal change directory via `$ cd ~/Downloads` (assuming the *rtmpt* source folder is in /Downloads)
  
    3. Execute:
    ```
    R CMD build rtmpt
    R CMD INSTALL rtmpt_<XYZ>.tar.gz
    ```
    where XYZ stands for the version number.

## Example
For an example download the [example](https://github.com/RaphaelHartmann/rtmpt-files/blob/master/an_example.zip) and run the R script. Note that this will take some time. You can also check out the help files for the `rtmpt` package in R via the commands
```
?to_rtmpt_model
?to_rtmpt_data
?set_params
?set_resps
?fit_rtmpt
?SimData
```
or the vignette via `vignette("rtmpt_intro")`.

## References
Klauer, K. C., & Kellen, D. (2018). RT-MPTs: Process models for response-time distributions based on multinomial processing trees with applications to recognition memory. *Journal of Mathematical Psychology, 82*, 111-130.
