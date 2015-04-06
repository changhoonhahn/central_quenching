pro build_treepm_sfr_dist, silent=silent
; Builds a sample of galaxies with SFR distribution using SFR main sequence
; from PRIMUS SFR-M* relation and mass distribution from Andrew's merger tree
    min_mass = 9.5
    max_mass = 11.5
    massbin = mass_bins(nmassbin, minmass=min_mass, maxmass=max_mass, delta_mass=0.01)
    if (keyword_set(silent) EQ 0 )  then struct_print, massbin
    z = 0.9     ; this is currently hardcoded and represent the starting redshift THINK TWICE BEFORE CHANGING
    mbin = 0.01

    ; import Andrew's merger tree snapshot on harmattan containing only central galaxies
    ; these subhalo merger tree files have M* derived from SHAM into fits table
    nsnap = get_nsnap_z(z)   ; n_snapshot of the snapshot closests in redshift to our redshift
    print, 'SNAPSHOT Number', nsnap

    treepm_dir = '/data1/hahn/wetzel_tree/'
    treepm_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(nsnap),2)+'.fits'  
    treepm_snap = mrdfits(treepm_dir+treepm_file,1)

    ; impose mass range
    mass_range = where((treepm_snap.mass GT min_mass) AND (treepm_snap.mass LT max_mass), n_massrange)  
    treepm_snap = treepm_snap[mass_range]
    has_child = where(treepm_snap.child GE 0, n_has_child)
    treepm_snap = treepm_snap[has_child]
    print, n_has_child, 'has children'

    tag_additions = replicate({SFR:0.0, SSFR:0.0, sfq:''},n_has_child)              ; add SFR, SSFR, and sfq tags
    treepm_snap = struct_addtags(treepm_snap, tag_additions)

    for i=0L,nmassbin-1L do begin 
        bin_qf = get_qf(massbin[i].massbin, z, /cosmosinterp)      ; Quiescent fraction of bin as a function of Mstar and z
        bin_mstar_index = where((treepm_snap.mass GE massbin[i].masslo) AND (treepm_snap.mass LT massbin[i].massup),ngal_bin) ; M* in mass bin 
        if (ngal_bin GT 0L) then begin 
            ngal_quiescent = round(float(ngal_bin)*bin_qf)
            ngal_active = ngal_bin-ngal_quiescent
            print, 'Ngal_bin =',ngal_bin ,'  Ngal_Q = ', ngal_quiescent,'    Ngal_SF = ', ngal_active
            
            ; Active galaxies
            if (ngal_active GT 0L) then begin 
                if (ngal_active LT ngal_bin) then active_index = cgrandomindices(ngal_bin,ngal_active) $     ; randomly selected bin_mstar indices that will be designated 'active'
                    else active_index = range(0,ngal_active-1L) 
                bin_active_index = bin_mstar_index[active_index]
                treepm_snap[bin_active_index].sfq = replicate('active',ngal_active)
                bin_active_mean_SFR = get_primus_sfr_mstar(massbin[i].massbin, z, sigma=bin_sigma_SFR, $
                    ssfr=bin_ssfr, minmass=min_mass, maxmass=max_mass, silent=silent, /harmattan)    ; hardcoded to run on harmattan
                treepm_snap[bin_active_index].SFR = bin_sigma_SFR*randomn(rseed,ngal_active,/normal)+bin_active_mean_SFR   ; log-normal sampling of SFR based on observed PRIMUS mean SFR and sigma SFR
                treepm_snap[bin_active_index].SSFR = treepm_snap[bin_active_index].SFR-treepm_snap[bin_active_index].mass
            endif 

            ; Quiescent galaxies   
            if (ngal_quiescent GT 0L) then begin
                bin_quiescent_index = bin_mstar_index[cmset_op(range(0,ngal_bin-1L), 'and', /not2, active_index)]    ; bin_mstar indices that are not bin_active_index
                treepm_snap[bin_quiescent_index].SSFR = 0.2*randomn(rseed,ngal_quiescent,/normal)-12.0     ; log-normal distribution of quiescent SSFR with mean ~10^-12 and 0.2 dex scatter
                treepm_snap[bin_quiescent_index].SFR  = treepm_snap[bin_quiescent_index].SSFR+treepm_snap[bin_quiescent_index].mass
                treepm_snap[bin_quiescent_index].sfq  = replicate('quiescent',ngal_quiescent)
            endif 
        endif 
    endfor
    output_dir  = '/data1/hahn/wetzel_tree/'
    output_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(nsnap),2)+'_ssfr_mbin'+$
        strmid(strtrim(string(mbin),2),0,4)+'.fits'  
    mwrfits, treepm_snap, output_dir+output_file, /create
return 
end
