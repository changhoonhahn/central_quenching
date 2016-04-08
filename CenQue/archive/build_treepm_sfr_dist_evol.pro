pro build_treepm_sfr_dist_evol, nsnap_i, nsnap_f, silent=silent, $
    instant=instant, constant=constant, linear=linear, mass_resid=mass_resid
; build SFR distribution for Andrew's merger tree given snapshot zi 
; delta redshift and cosmic time between the two snapshots
    ; define mass bins
    min_mass = 9.5 
    max_mass = 11.5
    mbin = 0.01
    mbins = mass_bins(nbin, minmass=min_mass, maxmass=max_mass, delta_mass=mbin)
    
    ; loop through snapshots 
    for n=0L,nsnap_i-nsnap_f-1L do begin 
        zi = nsnap_i-n 
        zf = zi-1L 
        zinitial = get_z_nsnap(zi,tcosmic=tinitial) 
        zfinal = get_z_nsnap(zf, tcosmic=tfinal) 
        print, 'Evolving from snapshot ', zi,'   to ', zf 
        print, 'Evolving from redhsift ', zinitial, '   to ', zfinal 

        ; set quenching efold flag
        if keyword_set(instant) then quench_str = '_instanttau'
        if keyword_set(constant) then quench_str = '_constanttau'
        if keyword_set(linear) then quench_str = '_lineartau'

        ; import mass, SFR, and SSFR data from snapshot zi (currently hardcoded) (parent)
        massSFR_dir  = '/data1/hahn/wetzel_tree/'
        if (zi EQ nsnap_i) then begin 
            massSFR_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(zi),2)+'_ssfr_mbin'+$
           strmid(strtrim(string(mbin),2),0,4)+'.fits' 
        endif else begin 
           massSFR_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(zi),2)+'_ssfr_mbin'+$
           strmid(strtrim(string(mbin),2),0,4)+'_evolfrom_snapshot'+strtrim(string(nsnap_i),2)+quench_str+'.fits'
        endelse 

        ; test snapshot file
        if file_test(massSFR_dir+massSFR_file) then begin
            massSFR_data = mrdfits(massSFR_dir+massSFR_file,1,silent=silent)
        endif else begin 
            stop 
        endelse 

        ; import SHAM mass from snapshot zf (child)
        massSFR_evol_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(zf),2)+'.fits'
        if file_test(massSFR_dir+massSFR_evol_file) then begin 
            massSFR_evol = mrdfits(massSFR_dir+massSFR_evol_file,1,silent=silent) 
        endif else begin 
            STOP
        endelse 
        ; add tags to snapshot zf data 
        tag_additions = replicate({SFR:0.0, SSFR:0.0, sfq:'', mass_parent:0.0, SFR_parent:0.0, tau:0.0, qSSFR:0.0}, $
            n_elements(massSFR_evol))
        massSFR_evol = struct_addtags(massSFR_evol, tag_additions)

        ; --------------------------------parent and child match--------------------------------------------------
        ; match child's index to parent's child and parent's index to child's parent
        parentmatch = cmset_op(massSFR_data.index,'and', massSFR_evol.parent, /index)   ; matched with index of snapshot zi centrals
        parent_cond = lonarr(n_elements(massSFR_data))
        parent_cond[parentmatch] = 1L

        childmatch = cmset_op(massSFR_evol.parent,'and',massSFR_data.index, /index)
        child_cond = lonarr(n_elements(massSFR_evol))
        child_cond[childmatch] = 1L

        if (NOT keyword_set(silent)) then begin
            print, n_elements(childmatch), n_elements(parentmatch), n_elements(massSFR_evol)
            print, 100.*float(n_elements(massSFR_evol)-n_elements(childmatch))/float(n_elements(massSFR_evol)), '% of galaxies are orphans'
            print, 1.0-float(total(child_cond))/float(n_elements(child_cond))
            massrange = where((massSFR_evol.mass GT min_mass) and (massSFR_evol.mass LE max_mass), n_massrange)
            print, 100.*float(n_massrange)/float(n_elements(massSFR_evol)), '% of galaxies are within mass range'
            print, 100.*(1.0-float(total(child_cond))/float(n_massrange)), '% of galaxies within mass range are orphans'
        endif 

        ; match active child to active parent and quiescent child to quiescent parent  
        ; match SF/Q classification 
        massSFR_evol[childmatch].sfq = massSFR_data[parentmatch].sfq
        matchcheck = where(massSFR_evol[childmatch].index NE massSFR_data[parentmatch].child, n_mismatch)

        ; child/parent active/quiescent indices 
        child_active_gal        = where(child_cond AND (strtrim(massSFR_evol.sfq,2) EQ 'active'), ngal_child_active)
        child_quiescent_gal     = where(child_cond AND (strtrim(massSFR_evol.sfq,2) EQ 'quiescent'), ngal_child_quiescent)
        parent_active_gal       = where(parent_cond AND (strtrim(massSFR_data.sfq,2) EQ 'active'), ngal_parent_active)
        parent_quiescent_gal    = where(parent_cond AND (strtrim(massSFR_data.sfq,2) EQ 'quiescent'), ngal_parent_quiescent)

        child_active_sort       = child_active_gal[sort(massSFR_evol[child_active_gal].index)]
        child_quiescent_sort    = child_quiescent_gal[sort(massSFR_evol[child_quiescent_gal].index)]
        parent_active_sort      = parent_active_gal[sort(massSFR_data[parent_active_gal].child)]
        parent_quiescent_sort   = parent_quiescent_gal[sort(massSFR_data[parent_quiescent_gal].child)]
        
        ; save parent SFR and mass to child data structure for convenience
        massSFR_evol[childmatch].SFR_parent = massSFR_data[parentmatch].SFR
        massSFR_evol[childmatch].mass_parent = massSFR_data[parentmatch].mass
   
        ; investigate mass residual between SHAM and SFR integrated mass  
        if keyword_set(mass_resid) then begin 
            ; integrated SFR mass using RK4 
            integrated_mass = get_envcount_mass_active_evol(massSFR_data[parent_active_sort].mass, $
                massSFR_data[parent_active_sort].SFR, zinitial, zfinal, SFRf=evolvedSFR, /harmattan)
            sham_mass = massSFR_evol[child_active_sort].mass                ; SHAM mass
            delta_mass = sham_mass-integrated_mass                          ; SHAM mass - integrated SFR mass
            ; record the mass residuals between SHAM and integrated SFR mass
            openw, lun, massSFR_dir+'delta_mass_evolved_sham_snapshot'+strtrim(string(zf),2)+'.dat', /get_lun 
            for i=0L,n_elements(delta_mass)-1L do begin 
                printf, lun, delta_mass[i], sham_mass[i], integrated_mass[i]
            endfor 
            free_lun, lun 
        endif 
       
        ; ASSUME active galaxy SFR evolves by the same amount for all SFR, retains same scatter in SFR as z_initial
        ; delta SFR = SFR(M_child, z_final) - SFR(M_parent, z_initial)      (for active galaxies only)
        dSFR = get_envcount_sfr_mstar_fit(massSFR_evol[child_active_sort].mass, zfinal, $
            /fixslope, silent=silent, /harmattan)-$
            get_envcount_sfr_mstar_fit(massSFR_data[parent_active_sort].mass, zinitial, $
            /fixslope, silent=silent, /harmattan)
        massSFR_evol[child_active_sort].SFR  = massSFR_data[parent_active_sort].SFR+dSFR        ; "evolve" SFR
        massSFR_evol[child_active_sort].SSFR = massSFR_evol[child_active_sort].SFR-$            ; compute SSFR
            massSFR_evol[child_active_sort].mass

        if (zi EQ nsnap_i) then begin   ; does not have tau assignment yet 
            massSFR_evol[child_quiescent_sort].SSFR = massSFR_data[parent_quiescent_sort].SSFR 
            massSFR_evol[child_quiescent_sort].SFR = massSFR_evol[child_quiescent_sort].SSFR+$
                massSFR_evol[child_quiescent_sort].mass
        endif else begin                ; use tau to quench galaxies appropriately 
            ; keep quenching time scale tau
            massSFR_evol[childmatch].tau = massSFR_data[parentmatch].tau 
            ; keep pre-determined post-quenched SSFR
            massSFR_evol[childmatch].qSSFR = massSFR_data[parentmatch].qSSFR            

            ; galaxies that are still quenching (tau > 0.01 and post-quenched qSSFR < SSFR) 
            stillquenching = where((massSFR_evol[child_quiescent_sort].tau GT 0.001) and $
                (massSFR_evol[child_quiescent_sort].qSSFR LT massSFR_data[parent_quiescent_sort].SSFR), n_stillquenching)

            if (n_stillquenching GT 0L) then begin 
                ; continue quenching 
                ; log(e^(-tfinal-tinial)/tau)
                massSFR_evol[child_quiescent_sort[stillquenching]].SFR = massSFR_data[parent_quiescent_sort[stillquenching]].SFR+$
                    alog10(exp(-(tfinal-tinitial)/massSFR_data[parent_quiescent_sort[stillquenching]].tau)) 
                massSFR_evol[child_quiescent_sort[stillquenching]].SSFR = $
                    massSFR_evol[child_quiescent_sort[stillquenching]].SFR-massSFR_evol[child_quiescent_sort[stillquenching]].mass

                ; galaxies that would quench too much. Galaxies that would overquench within this evolution step  
                overquenched = where(massSFR_evol[child_quiescent_sort[stillquenching]].SSFR LT $
                    massSFR_evol[child_quiescent_sort[stillquenching]].qSSFR, n_overquenched) 
                if (n_overquenched GT 0L) then begin 
                    ; if quenching the galaxy from t_initial to t_final 'over-quenches' then we simply set the child SSFR to the 
                    ; pre-determined post-quenching SSFR
                    massSFR_evol[child_quiescent_sort[stillquenching[overquenched]]].SSFR = $
                        massSFR_evol[child_quiescent_sort[stillquenching[overquenched]]].qSSFR
                    massSFR_evol[child_quiescent_sort[stillquenching[overquenched]]].tau = 0.0 
                endif  
            endif 

            ; fully quenched quiescent galaxies 
            ; these galaxy remain at the same SSFR and their SFR only change because SHAM masses change
            quenched = where(massSFR_evol[child_quiescent_sort].tau LE 0.001)    ; (quiescent galaxies from first snapshot) 
            massSFR_evol[child_quiescent_sort[quenched]].SSFR = massSFR_data[parent_quiescent_sort[quenched]].SSFR
            ; SFR needs to be updated since mass is upgraded by SHAM still
            massSFR_evol[child_quiescent_sort[quenched]].SFR = massSFR_data[parent_quiescent_sort[quenched]].SSFR+ $
                massSFR_data[parent_quiescent_sort[quenched]].mass
        endelse 

; --------------------- Quench active galaxies in order to match the quiescent fraction at z_final ---------------------
        ngal_append = 0L
        for i=0L,nbin-1L do begin 
            bin_qf = get_qf(mbins[i].massbin, zfinal, /cosmosinterp)  ; quiescent fraction at mass bin and zfinal (from COSMOS)
            bin_mass = where((massSFR_evol.mass GE mbins[i].masslo) AND (massSFR_evol.mass LT mbins[i].massup), ngal_bin) 
            massSFR_evol_bin = massSFR_evol[bin_mass] 

            ; active and quiescent based on parent classification
            bin_active = where(strtrim(massSFR_evol_bin.sfq,2) EQ 'active', ngal_bin_active0)
            bin_quiescent = where(strtrim(massSFR_evol_bin.sfq,2) EQ 'quiescent', ngal_bin_quiescent0)
            
            ; expected number of quiescent galaxies
            nexp_bin_quiescent = round(bin_qf*float(ngal_bin_active0+ngal_bin_quiescent0))      

            if (nexp_bin_quiescent GT 0L) then begin 
                n_quenched = nexp_bin_quiescent-ngal_bin_quiescent0         ; number of active galaxies to be quenched
                if (n_quenched GT 0L) then begin 
                    ; randomly select n_quenched active galaxies with parents to be quenched
                    ; and then quench them using quenching e-fold time tau 
                    rand_index = cgrandomindices(ngal_bin_active0, n_quenched) 
                    massSFR_evol_bin[bin_active[rand_index]].sfq = replicate('quiescent', n_quenched)   ; change to 'quiescent'

                    ; assign quenching timescale to active galaxies that are starting to quench 
                    massSFR_evol_bin[bin_active[rand_index]].tau = $
                        get_quenching_efold(massSFR_evol_bin[bin_active[rand_index]].mass_parent, $
                        instant=instant, constant=constant, linear=linear) 

                    ; start quenching (these are galaxies that started quenching at t_initial/z_initial 
                    massSFR_evol_bin[bin_active[rand_index]].SFR = massSFR_evol_bin[bin_active[rand_index]].SFR_parent+$
                        alog10(exp(-(tfinal-tinitial)/massSFR_evol_bin[bin_active[rand_index]].tau))
                    massSFR_evol_bin[bin_active[rand_index]].SSFR = massSFR_evol_bin[bin_active[rand_index]].SFR-$
                        massSFR_evol_bin[bin_active[rand_index]].mass

                    ; determine the final post-quenching SSFR of the quenched galaxy
                    massSFR_evol_bin[bin_active[rand_index]].qSSFR = 0.2*randomn(rseed, n_quenched, /normal)-12.0   

                    ; check to see if it overquenched somehow during first quenching time step  
                    active_overquenched = where(massSFR_evol_bin[bin_active[rand_index]].SSFR LT $
                        massSFR_evol_bin[bin_active[rand_index]].qSSFR, n_active_overquenched)  
                    if (n_active_overquenched GT 0L) then begin
                        print, 'some galaxies were overquenched right away at mass ', mbins[i].massbin
                        massSFR_evol_bin[bin_active[rand_index[active_overquenched]]].SSFR = $
                            massSFR_evol_bin[bin_active[rand_index[active_overquenched]]].qSSFR
                        massSFR_evol_bin[bin_active[rand_index[active_overquenched]]].tau = 0.0
                    endif 

                    ; incase code breaks due to alog10
                    if (max(massSFR_evol_bin[bin_active[rand_index]].SFR NE massSFR_evol_bin[bin_active[rand_index]].SFR) EQ 1) then begin 
                        struct_print, massSFR_evol_bin[bin_active[rand_index]]
                    endif 
                endif 

                if (NOT keyword_set(silent)) then begin 
                    print, 'mass ', mbins[i].massbin
                    print, 'nbin_gal', ngal_bin, '  from parent: ', ngal_bin_active0+ngal_bin_quiescent0, '     extra: ', ngal_bin-(ngal_bin_active0+ngal_bin_quiescent0)
                    print, 'from parent: SF', ngal_bin_active0, '   Q', ngal_bin_quiescent0
                    print, mbins[i].massbin, 'extra percentage: ', float(ngal_bin-(ngal_bin_active0+ngal_bin_quiescent0))/float(ngal_bin)
                endif 
            endif 

; ----------------------------------------- Galaxy accounting (fudging numbers) ---------------------------------------------------
; due to failures in the child-parent tracking there are extra galaxies.
; these extra galaxies are assigned active/quiescent based on the QF in 
; the mass bin. These extra galaxies only constitute a small fraction
            bin_extra = where((strtrim(massSFR_evol_bin.sfq,2) NE 'active') AND (strtrim(massSFR_evol_bin.sfq,2) NE 'quiescent'), $
                ngal_extra)            ; neither active or quiescent because we couldn't find the parent

            if (ngal_extra GT 0) then begin 
                ngal_extra_quiescent = round(bin_qf*float(ngal_extra)) 
                ngal_extra_active = ngal_extra - ngal_extra_quiescent
                if (NOT keyword_set(silent)) then begin  
                    print, 'extra galaxies', ngal_extra, '   extra_active', ngal_extra_active, $
                        '   extra_quiescent', ngal_extra_quiescent
                endif 

                ; extra active galaxies
                if (ngal_extra_active GT 0) then begin 
                    ; randomly selected bin_mstar indices that will be designated 'active'
                    if (ngal_extra_quiescent GT 0L) then begin  
                        extra_active_index = cgrandomindices(ngal_extra, ngal_extra_active)
                    endif else begin 
                        extra_active_index = range(0, ngal_extra_active-1L)     ; all quiescent 
                    endelse 
                    bin_extra_active_index = bin_extra[extra_active_index]
                    massSFR_evol_bin[bin_extra_active_index].sfq = replicate('active', ngal_extra_active)
                    bin_extra_active_mean_SFR = get_sfr_mstar(mbins[i].massbin, zfinal, $          ; mean SFR(M_extra,active, zfinal)
                        sigma=bin_extra_sigma_SFR, ssfr=bin_extra_ssfr, $
                        /envcount, /silent, /harmattan)    ; hardcoded to run on harmattan
                    massSFR_evol_bin[bin_extra_active_index].SFR = $                                ; scatter SFR 
                        bin_extra_sigma_SFR*randomn(rseed,ngal_extra_active,/normal)+$
                        bin_extra_active_mean_SFR   
                    massSFR_evol_bin[bin_extra_active_index].SSFR = $
                        massSFR_evol_bin[bin_extra_active_index].SFR-$
                        massSFR_evol_bin[bin_extra_active_index].mass
                endif else begin 
                    extra_active_index = []
                endelse 
                ; extra quiescent galaxies
                if (ngal_extra_quiescent GT 0L) then begin
                    bin_extra_quiescent_index = bin_extra[cmset_op(range(0,ngal_extra-1L), 'and', /not2, extra_active_index)]    
                    massSFR_evol_bin[bin_extra_quiescent_index].sfq  = replicate('quiescent',ngal_extra_quiescent)
                    ; SSFR for quiescent is a gaussian around -12.0
                    massSFR_evol_bin[bin_extra_quiescent_index].SSFR = 0.2*randomn(rseed,ngal_extra_quiescent,/normal)-12.0     
                    massSFR_evol_bin[bin_extra_quiescent_index].SFR  = massSFR_evol_bin[bin_extra_quiescent_index].SSFR+$
                        massSFR_evol_bin[bin_extra_quiescent_index].mass
                endif
            endif 
            massSFR_evol[bin_mass].sfq = massSFR_evol_bin.sfq 
            massSFR_evol[bin_mass].SFR = massSFR_evol_bin.SFR 
            massSFR_evol[bin_mass].SSFR = massSFR_evol_bin.SSFR 
            massSFR_evol[bin_mass].tau = massSFR_evol_bin.tau 
            massSFR_evol[bin_mass].qSSFR = massSFR_evol_bin.qSSFR 
        endfor
    ; remove galaxies outside of min and max mass range 
        mass_range = where((massSFR_evol.mass GE min_mass) AND $
            (massSFR_evol.mass LE max_mass), n_output)   
        output = massSFR_evol[mass_range]
        output_dir = '/data1/hahn/wetzel_tree/'
        output_file = $
            'subhalo_sham_centrals_snapshot'+strtrim(string(zf),2)+'_ssfr_mbin'+$
            strmid(strtrim(string(mbin),2),0,4)+'_evolfrom_snapshot'+$
            strtrim(string(nsnap_i),2)+quench_str+'.fits'
        mwrfits, output, output_dir+output_file, /create
    endfor 
return 
end
