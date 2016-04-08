pro plot_wetzel_cosmos_qf_comp 
    z = [0.36, 0.66, 0.88]
;    types = 'all' ;'cen';['all']; , 'cen', 'sat']
    types = 'cen'
    exp_sigma = [1.1972271, 1.05830526, 0.9182575]
    gauss_sigma = [0.89393588, 0.81232316, 0.70241647]
    for i=0L,n_elements(z)-1L do begin 
        psdir = get_figure_path(/tinker) 
        psfile = psdir+'plot_qf_z'+strmid(strtrim(string(z[i]),2),0,4)+types+'.eps'

        im_plotconfig, 0, pos, psfile=psfile, charsize=1.8
    ; import cosmos data
        cosmos_dir = '/data1/hahn/wetzel_tree/'
        cosmos_file = 'qf_z'+strmid(strtrim(string(z[i]),2),0,4)+types+'.dat'
        readcol, cosmos_dir+cosmos_file, mass, qf, qflow, qfhigh

        xrange = [min(mass), max(mass)] 
        yrange = [0.0, 1.0]
        
        plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
            xrange=xrange, yrange=yrange, xtickinterval=1.0, ytickinterval=0.25,$
            xtitle=textoidl('log M_{*}'), ytitle=textoidl('f_{Q}'), color=im_color('black',255)
        im_legend, ['COSMOS', 'Wetzel', 'Exp fit'], /top, /left, box=0,$
            color=['dodger blue', 'red', 'black'], psym=replicate(16,3), $
            symsize=replicate(1.5,3)*1.7, $
            symthick=8, spacing=2.7, charsize=2.1
        im_legend, textoidl('z \sim'+strmid(strtrim(string(z[i]),2),0,4)), $
            /bottom, /left, box=0, charsize=2.1
        djs_oplot, mass, qf, psym=symcat(16, thick=5), color=im_color('dodger blue', 255) 
        djs_oplot, mass, exp((mass-12.0)/exp_sigma[i]), line=0, linewidth=5, $
            color=im_color('black', 255)
;        djs_oplot, mass, exp(-0.5*(mass-12.0)^2/gauss_sigma[i]), line=2, linewidth=5,$
;            color=im_color('black', 255)
        
        massrange = where(mass GT 9.5 and mass LE 11.5) 
        mass = mass[massrange] 
        qf = qf[massrange] 
        wetzel = fltarr(n_elements(qf))
        for ii=0L,n_elements(qf)-1L do begin 
            wetzel[ii] = get_qf(mass[ii], z[i])
        endfor 
        djs_oplot, mass, wetzel, psym=symcat(16, thick=5), color=im_color('red', 255) 
        im_plotconfig, /psclose, psfile=psfile
    endfor
return 
end 
