pro plot_sfr_activeseq_powerlaw_comp 
; Plots the evolution of SFR as a function of cosmic time 
; in order to compare two methods of SFR evolution: 
; a simple power law evolution or an SFR obtained by a
; parameterization of the active sequence where M* is varied 
; by SFR integration (essentially SFR decreases as M* accretes mass
; and over cosmic time)
; configure initial conditions
    mass = [9.5, 10.0, 10.5, 11.0] 
    zform = 1.0
    zbin = 0.9
    SFR0 = fltarr(n_elements(mass))
    for i=0L,n_elements(mass)-1L do SFR0[i] = get_envcount_sfr_mstar_fit(mass[i], zbin, /fixslope, /silent, /harmattan) 
    print, mass
    print, SFR0
; configure plot
    psdir = get_figure_path(/tinker) 
    psfile = psdir+'sfr_activeseq_powerlaw_comp.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8 
; load colors
    loadct,24, /silent, file='/home/users/hahn/research/pro/im_pro/etc/fsc_brewer.tbl'
    mass_color = reverse([250, 223, 197, 180])
    mass_label = strarr(n_elements(mass)) 
    for i=0L,n_elements(mass)-1L do mass_label[i] = textoidl('M_{*} \sim')+' '+strmid(strtrim(string(mass[i]),2),0,4)
    xrange = [0.3*10.0^9, 7.0*10.0^9]
    yrange = [-0.2, 1.2]
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        xrange=xrange, yrange=yrange, xtickinterval=2*10.0^9.0, ytickinterval=0.25,$ 
        xtitle='Cosmic Time', ytitle=textoidl('log(SFR)'), color=im_color('black',255)
    im_legend, mass_label, /top, /right, box=0,$
        color=mass_color, psym=replicate(16,n_elements(mass_label)), $
        symsize=replicate(1.5,n_elements(mass_label))*1.7, $
        symthick=8, spacing=2.7, charsize=2.1
    tbin = galage(zbin, zform)
    for i=0L,n_elements(mass)-1L do begin 
        djs_oplot, [tbin], [SFR0[i]], psym=symcat(16, thick=5), $
            symsize=2, color=im_color(mass_color[i], 255) 
    endfor 
    while (zbin GT 0.101) do begin 
        print, zbin 
        zbin = zbin-0.2
        tbin = galage(zbin, zform)
        for j=0L,n_elements(mass)-1L do begin 
            delta_mass = get_envcount_mass_active_evol(mass[j], SFR0[j], zbin+0.2, zbin, SFRf=SFRf, /harmattan) 
            mass[j] = delta_mass 
            SFR0[j] = SFRf
            djs_oplot, [tbin], [SFR0[j]], psym=symcat(16, thick=5), $
                symsize=2, color=im_color(mass_color[j], 255) 
        endfor 
        print, mass
        print, SFR0
    endwhile 
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end
