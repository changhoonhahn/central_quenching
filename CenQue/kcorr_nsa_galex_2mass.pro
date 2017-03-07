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
    
    ; store matched AB fluxes to 'mgys' 
    mgys = fltarr(10, n_elements(m_vagc_nsa)) 
    mgys_ivar = fltarr(10, n_elements(m_vagc_nsa)) 
        
    ; UV bands
    fnuv = nsa.nmgy[0:1,*]
    fnuv_ivar = nsa.nmgy_ivar[0:1,*]
    mgys[0:1,*] = fnuv[*,m_nsa_vagc]*1.e-9
    mgys_ivar[0:1,*] = fnuv_ivar[*,m_nsa_vagc]*1.e+18

    ; SDSS bands
    if (sdss_flux eq 'sersic') then begin 
        ugriz = nsa.nmgy[2:6,*]
        ugriz_ivar = nsa.nmgy_ivar[2:6,*]
        mgys[2:6,*] = ugriz[*,m_nsa_vagc]*1.e-9
        mgys_ivar[2:6,*] = ugriz_ivar[*,m_nsa_vagc]*1.e+18
    endif 
    if (sdss_flux eq 'petro') then begin
        ugriz = nsa.petroflux[2:6,*]
        ugriz_ivar = nsa.petroflux_ivar[2:6,*]
        extinction = nsa.extinction[2:6,*]
        mgys[2:6,*] = ugriz[*,m_nsa_vagc]*10.^(0.4*extinction[*,m_nsa_vagc])*1.e-9
        mgys[2:6,*] = ugriz_ivar[*,m_nsa_vagc]*10.^(-0.8*extinction[*,m_nsa_vagc])*1.e+18
    endif 

    ; 2MASS bands
    jhk = vagc.abmaggies[5:7,*]
    jhk_ivar = vagc.abmaggies_ivar[5:7,*]
    mgys[7:9,*] = jhk[*,m_vagc_nsa]
    mgys_ivar[7:9,*] = jhk_ivar[*,m_vagc_nsa]

    ; redshifts
    zz = nsa.z
    z_red = zz[m_nsa_vagc]
    
    ; filter list
    filterlist=['galex_FUV.par', 'galex_NUV.par', 'sdss_u0.par', $
        'sdss_g0.par', 'sdss_r0.par', 'sdss_i0.par', 'sdss_z0.par', $
        'twomass_J.par', 'twomass_H.par', 'twomass_Ks.par']

    ; call kcorrect 
    kcorrect, mgys, mgys_ivar, z_red, kcorr, band_shift=0.0, $
        coeffs=k_coeff, rmaggies=k_maggies, mass=k_mass, $
        mtol=k_mtol, absmag=k_absmag, amivar=k_absmag_ivar, $
        filterlist=filterlist

    out_struct = replicate({ $
        ra:                 -999.0, $
        dec:                -999.0, $
        z:                  -999.0, $
        maggies:        fltarr(10), $
        maggies_ivar:   fltarr(10), $
        k_mass:             -999.0, $
        k_mtol:         fltarr(10), $
        k_coeff:         fltarr(5), $
        k_maggies:      fltarr(10), $
        absmag:         fltarr(10), $
        absmag_ivar:    fltarr(10), $ 
        kcorrect:       fltarr(10)},$
        n_elements(m_vagc_nsa))

    ra = nsa.ra
    ra = ra[m_nsa_vagc]
    dec = nsa.dec
    dec = dec[m_nsa_vagc]

    out_struct.ra              = ra
    out_struct.dec             = dec
    out_struct.z               = z_red
    out_struct.maggies         = mgys
    out_struct.maggies_ivar    = mgys_ivar
    out_struct.k_mass          = k_mass
    out_struct.k_mtol          = k_mtol
    out_struct.k_coeff         = k_coeff
    out_struct.k_maggies       = k_maggies
    out_struct.absmag          = k_absmag
    out_struct.absmag_ivar     = k_absmag_ivar
    out_struct.kcorrect        = kcorr
    
    ;out_file = '/mount/sirocco1/hahn/cenque/observations/nsa_galex_2mass_kcorrect_'+sdss_flux+'.fits'
    ;print, out_file

    mwrfits, out_struct, 'catalog.fits', /CREATE
end 
