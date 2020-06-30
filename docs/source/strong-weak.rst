Strong-weak partitioning 
========================

The line list (.trans) can be partition into strong/weak lines parts using a reference absorption intensity threshold. 
This option is invoked using a combination of two keywords: `TRANS` as a type of `profile` (main body keyword) and 
a `cutoff` type `exp`. The reference temperature is defined using `Temperature`. 
The dynamic partitioning cutoff is given by 

:math:`I_{\rm cut} = I0 \exp(-\nu/\alpha)` 

The format of the `cutoff` line is as follows
::
    
    cutoff  Exp 1e-25  alpha  2000   
        

where :math`I0=1e-25` and :math:`\alpha=2000`. 


For example: 
::
    
    Temperature  3000.0 
    Range 0.0  10000.0
    
    absorption
    
    cutoff  Exp 1e-25  alpha  2000
    
    trans
    
    output part

    States       H2O.states
    Transitions  H2O.trans
        


This input will creates two files H2O.trans.part.w and H2O.trans.part.s, with the weak and strong lines in 
the standard .trans format. In case of multiple .trans files in `Transitions` section, each of .trans  
generates two .w and .s files. 


 
