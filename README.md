# Fourier ptychographic microscopy (FPM)
FPM is a computational microscopy technique that allows obtaining microscopic image of a measured sample with much higher resolution than it would follow from the numerical aperture (NA) of used microscope objective. Such high resolution is obtained by combining in the Fourier domain information about the measured specimen from different illumination directions. This advantage combined with the fact that low NA objectives are characterized by low magnification, allows obtaining a high-resolution image of a specimen with a large field of view. Moreover, iterative algorithms used in the process of reconstruction, enable obtaining information not only about the amplitude but also about the phase of the sample The phase distribution is especially significant in the case of bioimaging, as living cells are known to be semi-transparent objects with low amplitude, hence difficult to see under the classical brightfield microscope).

# _FPM app_
_FPM app_ is first, to the best of our knowledge: simple, intuitive, universal, semi-automatic, open to modify, GUI open-source app that allows users to perform straightforward reconstruction of original FPM datasets, without requiring the user to have specialized knowledge in the field of imaging or programming. _FPM app_ was created in MATLAB (based on Lei Tian algorithm - http://sites.bu.edu/tianlab/open-source/) and all MATLAB codes (that are free to modify for personal use) along with standalone executable version and full _FPM app_ documentation can be found in https://github.com/MRogalski96/FPM-app/releases.

# Software
_FPM app_ is released in 2 versions:
- MATLAB version - it contains a pack of MATLAB codes, which are used to open FPM app through MATLAB. These codes are open to be modified, to adjust _FPM app_ to a given set of preferences or to further improve it.
- Executable version – it contains FPMAppInstaller_web.exe that installs FPMapp.exe along with all the necessary files and MATLAB Runtime that is required to run MATLAB standalone applications.

# Documentation
To _FPM app_ is attached documentation that consists: 
- Description of _FPM app_ working overflow and all _FPM app_ windows
- Methods used for processing the FPM data
- Description of our innovations into the FPM processing path
- A tips for modifying _FPM app_
- Exemplary description of collecting and reconstructing FPM data process

# Author information
_FPM app_ was developed by **Mikolaj Rogalski** as a part of master's thesis at Faculty of Mechatronics, Warsaw University of Technology, Warsaw, Poland

# How to cite the work
Mikołaj Rogalski, Piotr Zdańkowski, Maciej Trusiak, FPM app: an open-source MATLAB application for simple and intuitive Fourier ptychographic reconstruction, Bioinformatics, 2021;, btab237, https://doi.org/10.1093/bioinformatics/btab237

# Contact
In case of any problem with _FPM app_ please contact the author:
mikolaj.rogalski.dokt@pw.edu.pl
