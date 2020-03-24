Using filters
=============

Cross-sections using Gaussian profile where only the states labelled A2Pi (upper) and X2Pi (lower) are selected. 

Example::

     
    (PS cross-sections)
    
    Temperature 2000.0
    Range 0.  45000.0
    
    Npoints 45001
    
    filter
     upper 7 A2Pi
     lower 7 X2Pi
    end
    
    emission
    gaussian
    hwhm 1.0
    
    output PS_emis_2000_AX
    
    States      PS.states
    Transitions  PS.trans
    


The filter command has the following structure:

    upper col Label 
    lower col Label 


where `upper` and `lower` refer to the upper and lower states, respectively; `label` is the reference label (number or string) which appear in the `col` -column  (integer, 1,2,3...) . 


An example of a multiple QNs-filter::

    filter
     upper 6 0  
     upper 7 0  
     upper 8 0
     upper 9 1 
     upper 10 1 
     upper 11 1 
     upper 12 0 
     upper 13 0 
     upper 14 0
     lower 6 0  
     upper 7 1  
     upper 8 1  
     upper 9 0 
     upper 10 0 
     upper 11 0 
     upper 12 0 
     upper 13 0 
     upper 14 0
    end


Here the first number in each pair indicates the column number in the .states file (counting from the 1st column with the state index), while the second is the value used in the filter of the state. In this case this is to select the vibrational band nu3 (0 00 111 000 <- 0 00 000 000) of methane. 

>Filters are case-sensitive.




