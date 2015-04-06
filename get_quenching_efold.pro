function get_quenching_efold, mstar, instant=instant, constant=constant, linear=linear
; calculates the characteristic e-fold time for quenching
; as a function of mass using different models 
    if keyword_set(instant) then begin 
        ; near zero. super fast quenching 
        tau = 0.0000001 
    endif 
    if keyword_set(constant) then begin
        tau = replicate(0.8, n_elements(mstar))   ; Gyr
    endif 
    if keyword_set(linear) then begin 
        tau = -(0.9/1.5)*(mstar-9.5)+1.0
        if (min(tau) LT 0.1) then begin 
            for i=0L,n_elements(tau)-1L do begin 
                if (tau[i] LT 0.1) then tau[i]=0.1
            endfor 
        endif 
    endif 
    return, tau
end
