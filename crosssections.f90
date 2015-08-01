   program crosssections

    use timer
    use accuracy
    !
    use spectrum,only : Readinput,verbose,intensity
     !
     if (verbose>=1) call TimerStart('Spectrum')
     !
     ! initializing some constants 
     !
     call accuracyInitialize
     !
     ! read input data
     !
     call ReadInput
     !
     call intensity
     !
     if (verbose>=1) call MemoryReport
     !
     if (verbose>=1) call TimerStop('Spectrum')
     !
     if (verbose>=1) call TimerReport
     !
    end program crosssections
 
