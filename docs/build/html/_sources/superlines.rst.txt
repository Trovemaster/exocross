Super-lines
===========

ExoCross can use the precomputed _binned_ cross-sections (i.e. with keyword `bin`) in order to put a line-profile on such a _binned_ transition treating it as a single line with a combined intensity. It is therefore assumed that the corresponding binned-cross-sections are pre-computed for a given T. The advantage of this approach is that calculations of the cross-sections for line-profiles are much faster due to the smaller number of bins relative to the number of lines. 

The approximation used is that the line-centres for all transitions within a given bin sifted to the centre of the bin. 

This option can be invoked by giving a keyword `histogram` (or `histogram-J`, see below) anywhere in the body of the inout file. The states file in this case is used only for the partition function. The histograms (binned-cross-sections) are then given as usual transition files. Otherwise the input file is as usual. 


For example: 
::
    
    Temperature  296.0
    Range 0.0  10000.0
        
    Npoints 10001
    
    absorption
    gauss
    hwhm 0.5 (cm-1)
    
    histrogram
    
    output C2H2_ab_g0.5
    States       H2O.states
    Transitions  H2O_T1000K_bin.xsec
        


In order to be able using J-dependent Voigt line-profiles, a J-value can be given together with the transition filename. In this case the histrogramJ keyword has to be used instead of histogram.   


For example: 
::
    
    histogram-J
    
    Temperature  1000.0
    pressure 1
    Range 0.0  10000.0
    
    Npoints 1000001
    
    species.
      He  gamma 0.0043 n 0.02  t0 298.0 file  1H2-16O__Nina__He__a1.broad model JJ ratio 0.16
      H2  gamma 0.0207 n 0.027 t0 298.0 file  1H2-16O__Nina__H2__a1.broad model JJ ratio 0.84
    end
    
    absorption
    voi-quad
    mass 18
    
    output H2O_voi-q_1000K_1bar
    offset 25.0
    nquad  40
    
    States  1H2-16O__BT2.states
    
    Transitions
      bin-J0_T3500_grid1_1e-4.xsec     0
      bin-J1_T3500_grid1_1e-4.xsec     1
      bin-J2_T3500_grid1_1e-4.xsec     2
      bin-J3_T3500_grid1_1e-4.xsec     3
      bin-J4_T3500_grid1_1e-4.xsec     4
      bin-J5_T3500_grid1_1e-4.xsec     5
    end

            
A histogram can be produced using BIN as the "line-profile" type. This can be combined with the wavelength-range used instead of wavenumber-range (default) by providing `um` or `micron` as units in RANGE. `Bin' is used to bin all intensities into different wavenumber/wavelength grid points.  
::

    Temperature  296.0
    Range 1.0  1000.0  um 
    
    Npoints 10001
    
    absorption
    bin
    


 
