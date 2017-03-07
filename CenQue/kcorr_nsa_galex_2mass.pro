; IDL code to combine NSA SDSS, GALEX UV, and 2MASS IR photometry
; and then run k-correct on them. 
pro kcorr_nsa_galex_2mass, sdss_flux

    ; import NSA DRP catalog (contains elliptical petrosian fluxes for 
    ; SDSS and also GALEX FNUV) 
    ;nsa = mrdfits('dat/observations/drpall-v2_0_1.fits', 1) ; NSA DRP catalog for MaNGA 
    
    ; import NSA catalogs (contains sersic fluxes for SDSS and also GALEX FN UV) 
    nsa = mrdfits('dat/observations/nsa_v0_1_2.fits', 1)

    ; import AB fluxes for 2MASS from 
    ; SDSS DR7 VAGC kcorr file 
    vagc = mrdfits('dat/observations/kcorrect.none.petro.z0.00.fits', 1) 

    ; match the catalogs together using spherematch and construct combined 
    ; AB fluxes 
    
    ; run spherematch on them with match radius of 3'' (in degrees)
    spherematch, vagc.ra, vagc.dec, nsa.ra, nsa.dec, 0.0003, m_vagc_nsa, m_nsa_vagc, d_match
    
    print, n_elements(m_vagc_nsa)
    print, n_elements(m_nsa_vagc)
    
    ; maggies
    mgys = fltarr(10, n_elements(m_vagc_nsa)) 
        
    ; UV bands
    mgys[0,*] = nsa.nmgy[0,*][m_nsa_vagc]
    mgys[1,*] = nsa.nmgy[1,*][m_nsa_vagc]

    ; sdss bands
    if (sdss_flux eq 'sersic') then begin 
        for i_b=2L, 6L  do begin 
            mgys[i_b,*] = nsa.nmgy[i_b,*][m_nsa_vagc]
        endfor 
    endif 
    if (sdss_flux eq 'petro') then begin
        for i_b=2L, 6L  do begin 
            mgys[i_b,*] = nsa.petroflux[i_b,*][m_nsa_vagc] * nsa.extinction[i_b,*][m_nsa_vagc]
        endfor 
    endif 

    ; 2MASS bands
    mgys[7,*] = vagc.abmaggies[5,*][m_vagc_nsa]
    mgys[8,*] = vagc.abmaggies[6,*][m_vagc_nsa]
    mgys[9,*] = vagc.abmaggies[7,*][m_vagc_nsa]
end 
