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
    ......
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


Uncertainty filter 
------------------

For the line lsits with uncertainties (column 5), the uncertainty filter (threshold) 
can be used to select accurate transitions only. Both states (upper and lower) must satisfy the 
uncertainty criteria. The uncertainly keywords are `Unc` or `uncertainty` with the uncertainty value as
the parameter. This should be specified in the `filter` section. 


Example:

    filter
     unc 0.001
    end


The default column with uncertainty is number 5, which can be changed in the QN section using the keyword `unc`:


    QN
     unc 6
    end



Using QN 
--------


A QN section has been added which explicitly tells ExoCross program which columns in the .states file it should read for QNs.  
The keywords for the QN section include `density` and `vib`, followed by one and three numbers respectively.  
The one number following the `density` keyword is the column which contains the density information, 
usually the last column unless edited by hand.  The `vib` keyword is followed by two numbers denoting the first and last column number 
for the last vibrational quantum number used.  For example, if the density column is the 18th column in the .states file 
and the vibrational normal mode quantum numbers v1, v2, v3 are in columns 7, 8 and 9 then the QN section would appear as::



    QN
     density 18
     vib 7 10 
    end

