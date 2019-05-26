# rtmpt
R package for fitting Response-Time extended Multinomial Processing Tree models

## Installation
### Windows

### Linux
Before installing the R package `rtmpt` you need to install the latest version of the GNU Scientific Library (GSL):
1. Download the latest version of [GSL](https://www.gnu.org/software/gsl/)
2. Unzip the file
3. Assuming you downloaded version 2.5 and unziped the file to ~/Downloads/gsl-2.5, open terminal and change directory via ```$ cd ~/Downloads/gsl-2.5```
4. Execute the following lines:
```
./configure
make
sudo make install
```
After the installation of GSL is complete you can install `rtmpt`:
1. Download all the source files from the `rtmpt` package above
2. In the terminal change directory via `$ cd ~/Downloads` (assuming the *rtmpt* source folder is in /Downloads)
3. Execute:
```
R CMD build rtmpt
R CMD INSTALL rtmpt_XYZ.tar.gz
```
where XYZ stands for the version number.

If the error "libgsl.so.23: cannot open shared object file: No such file or directory" occurs while testing the `rtmpt` package, run also in the terminal `$ sudo ldconfig -v`

## References
Klauer, K. C., & Kellen, D. (2018). RT-MPTs: Process models for response-time distributions based on multinomial processing trees with applications to recognition memory. *Journal of Mathematical Psychology, 82*, 111-130.
