Computing Lifetimes 
===================

The lifetimes are computed by summing Einstein A coefficients as

:math:`\tau_{i'}= \sum_{i''< i'} \frac{1}{A_{i',i''}}`

For states with no channels to lower states, the lifetimes are considered undefined and assign NaN. This feature is based on the external 
library IEEE. In the older version of ExoCross, such states were assigned -1.0000. 



Example  of the input file for lifetimes:: 

    
    lifetime
    
    output 28SiO-16
    
    States 28Si-16O__EBJT.states
    Transitions 28Si-16O__EBJT.trans 
    
