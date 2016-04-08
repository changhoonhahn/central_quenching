pro plot_groupcatalog_ssfr_dist, Mrcut=Mrcut
; plots the SSFR distribution for Jeremy's SDSS group catalog 
    if (Mrcut eq 18) then masscut='9.4'            ; Mr cuts and mass cuts from Tinker et al. 2011
    if (Mrcut eq 19) then masscut='9.8' 
    if (Mrcut eq 20) then masscut='10.2' 

    psdir = get_figure_path(/tinker)            
    psfile = psdir+'massSFR_clf_groups_M'+strtrim(string(Mrcut),2)+'_'+masscut+'_D360.galdata_corr.central.eps'
    im_plotconfig, 5, pos, psfile=psfile, charsize=1.               ; set .eps file 
    !P.Multi = [0, 2, 2] 
    
    ; import ssfr distribution data (currently hardcoded) 
    ssfr_file_dir = '/data1/hahn/group_catalog/'
    ssfr_file = 'massSFR_clf_groups_M'+strtrim(string(Mrcut),2)+'_'+masscut+'_D360.galdata_corr.central.fits'
    print, ssfr_file_dir+ssfr_file
    ssfrdata = mrdfits(ssfr_file_dir+ssfr_file,1)

    massbin = mass_bins(nmassbin, /simple)              ; divide into four mass bins
    for i=0L,nmassbin-1L do begin 
        ; impose mass range on the ssfr data
        mass_range = where((ssfrdata.mass GE massbin[i].masslo) and (ssfrdata.mass LT massbin[i].massup), n_mass_range)
        print, massbin[i].masslo, '     to ', massbin[i].massup, '      ', n_mass_range 

        masslabel ='['+strmid(strtrim(string(massbin[i].masslo),2),0,4)+$
            ','+strmid(strtrim(string(massbin[i].massup),2),0,4)+']'
        if (n_mass_range GT 0L) then begin 
            ssfr_bin = ssfrdata[mass_range] 
            cgHistoplot, ssfr_bin.ssfr, BINSIZE=0.1, /frequency, /Fill, $
                xtitle=textoidl('log SSFR [yr^{-1}]'), ytitle=textoidl('P(log SSFR)'),$
                polycolor='white', datacolorname='black', xrange=[-13.0, -7.0], yrange=[0.0, 0.12]
        endif else begin 
            plot, [0], [0], /nodata, position=pos[*,0], xsty=0, ysty=0, $
                xrange=[-13.0, -7.0], yrange=[0.0, 0.12], $
                ;xtitle=textoidl('log SSFR [yr^{-1}]'), ytitle=textoidl('P(log SSFR)'),$
                color=im_color('black',255), charsize=2.1
            ;cgHistoplot, [], BINSIZE=0.1, /frequency, /Fill, $
            ;    xtitle=textoidl('log SSFR [yr^{-1}]'), ytitle=textoidl('P(log SSFR)'),$
            ;    polycolor='white', datacolorname='black', xrange=[-13.0, -7.0], yrange=[0.0, 0.12]
        endelse 
        if (i EQ 0L) then begin 
            im_legend, 'SDSS Group Catalog', box=0,charsize=1,$
                /bottom, /left,/normal
        endif 
        im_legend, textoidl(masslabel), box=0,charsize=2,$
            /top, /left, /normal
    endfor
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end
