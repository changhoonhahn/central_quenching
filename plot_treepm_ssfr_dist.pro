pro plot_treepm_ssfr_dist, nsnap=nsnap, constant=constant, linear=linear
; plots the SSFR distribution for Andrew's TreePM+SHAM with SFR assignment and evolution
; select mass bin
    massbin = mass_bins(nmassbin,/simple)
    mbin = 0.01
    nsnap_i = 13    ; original snapshot, hardcoded
    ; set quenching mechanism flag 
    if keyword_set(constant) then quench_str = '_constanttau'
    if keyword_set(linear) then quench_str = '_lineartau'

    ; one snapshot, or multiple
    if keyword_set(nsnap) then begin 
        nsnap_arr = [nsnap]
    endif else begin 
        nsnap_arr = range(1,13)
    endelse 
    
    ; set figure name 
    psdir = get_figure_path(/tinker) 
    if keyword_set(nsnap) then begin 
        if (nsnap EQ nsnap_i) then begin 
            psfile = psdir+'subhalo_sham_centrals_ssfr_dist_snapshot'+$
                strtrim(string(nsnap),2)+'_schechter'+quench_str+'.eps'
        endif else begin 
            psfile = psdir+'subhalo_sham_centrals_ssfr_dist_snapshot'+$
                strtrim(string(nsnap),2)+'_evolfrom_snapshot'+strtrim(string(nsnap_i),2)+$
                '_schechter'+quench_str+'.eps'
        endelse
    endif else begin 
        psfile = psdir+'subhalo_sham_centrals_ssfr_dist_allsnapshots'+$
            '_evolfrom_snapshot'+strtrim(string(nsnap_i),2)+'_schechter'+quench_str+'.eps'
    endelse
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.
    !P.Multi = [0, 2, 2] 

    for i_snap=0L,n_elements(nsnap_arr)-1L do begin 
        ; what redshift does the snapshot correspond to? 
        zsnap = get_z_nsnap(nsnap_arr[i_snap])

        ; import ssfr distribution data
        data_dir = '/data1/hahn/wetzel_tree/'
        if (nsnap_arr[i_snap] EQ nsnap_i) then begin 
            ssfr_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(nsnap_i),2)+'_ssfr_mbin'+$
                strmid(strtrim(string(mbin),2),0,4)+'.fits' 
        endif else begin
            ssfr_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(nsnap_arr[i_nsnap]),2)+'_ssfr_mbin'+$
                strmid(strtrim(string(mbin),2),0,4)+'_evolfrom_snapshot'+strtrim(string(nsnap_i),2)+$
                quench_str+'.fits'
        endelse 
        ssfrdata = mrdfits(data_dir+ssfr_file,1)        ; read treepm SSFR data

        for i=0L,nmassbin-1L do begin 
            ; divide SSFR distributions into 4 mass bins  
            mass_range = where((ssfrdata.mass GE massbin[i].masslo) and (ssfrdata.mass LT massbin[i].massup), n_mass_range)
            ssfr_bin = ssfrdata[mass_range] 
            masslabel ='['+strmid(strtrim(string(massbin[i].masslo),2),0,4)+$
                ','+strmid(strtrim(string(massbin[i].massup),2),0,4)+']'
            cgHistoplot, ssfr_bin.ssfr, BINSIZE=0.2, /frequency, /Fill, $
                xtitle=textoidl('log SSFR [yr^{-1}]'), ytitle=textoidl('P(log SSFR)'),$
                polycolor='white', datacolorname='black', xrange=[-13.0, -7.0], yrange=[0.0, 0.12]
            
            if (i EQ 0L) then begin 
                im_legend, 'Snapshot '+strtrim(string(nsnap_arr[i_snap]),2)+' '+$
                    textoidl('z='+strmid(strtrim(string(zsnap),2),0,4)), box=0,charsize=1,$
                    /bottom, /left,/normal
            endif 
            
            im_legend, textoidl(masslabel), box=0,charsize=2,$
                /top, /left, /normal
        endfor
    endfor 
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end
