SPECTRA
=======

ExoCross can be used with the SPECTRA line list format developed by iao.ru  by adding the `SPECTRA` keyword anywhere in the input file (outside any section). 
This will also require the definition of the (i) reference temperature, (ii)  partition function for the target temperature, (iii)
partition function for the reference temperature, and (iv) the molecule/isotope pair  (`iso`). 

SPECTRA broadening parameters (as part of the format) will be used unless the species-section is given, which specifies the broadening.  
The .states file is not required and ignored if given. The SPECTRA reference intensity is used directly to compute the intensity for the target temperature.

For example: 
::

    Temperature  500.0  ref 296.0
    Range 0.0  10000.0
    
    Npoints 10001
    
    absorption
    Voigt
    
    spectra
    iso 26 1
    pf 1000.0 ref 500
    
    output CH4_voigt_T500K

    Transitions  SpectraMol_CH4_296K.txt
    
    

