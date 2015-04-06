pro test_treepm_child_parent, zi
; build SFR distribution for Andrew's merger tree given snapshot zi 
; delta redshift and cosmic time between the two snapshots
    zf = zi-1L 
    zinitial = get_z_nsnap(zi,tcosmic=tinitial) 
    zfinal = get_z_nsnap(zf, tcosmic=tfinal) 
    delz = abs(zinitial-zfinal)
    delt = abs(tinitial-tfinal)
    mbin = 0.01
; import mass, SFR, and SSFR data from snapshot zi
    massSFR_dir  = '/data1/hahn/wetzel_tree/'
    massSFR_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(zi),2)+'_ssfr_mbin'+$
       strmid(strtrim(string(mbin),2),0,4)+'.fits' 
    if file_test(massSFR_dir+massSFR_file) then massSFR_data = mrdfits(massSFR_dir+massSFR_file,1) $
        else STOP

; import mass, SFR, and SSFR data from snapshot zf
    massSFR_evol_file = 'subhalo_sham_centrals_snapshot'+strtrim(string(zf),2)+'.fits'
    if file_test(massSFR_dir+massSFR_evol_file) then massSFR_evol = mrdfits(massSFR_dir+massSFR_evol_file,1) $
        else STOP
    tag_additions = replicate({SFR:0.0, SSFR:0.0, sfq:''},n_elements(massSFR_evol))
    massSFR_evol = struct_addtags(massSFR_evol, tag_additions)
; assign SFR and sfq to child galaxy 
    parentchild = cmset_op(massSFR_evol.index,'and',massSFR_data.child, /index)
    parentmatch = cmset_op(massSFR_data.index,'and', massSFR_evol.parent, /index)   ; matched with index of snapshot zi centrals
    parent_cond = lonarr(n_elements(massSFR_data))
    parent_cond[parentmatch] = 1L 
    childmatch = cmset_op(massSFR_evol.parent,'and',massSFR_data.index, /index)
    child_cond = lonarr(n_elements(massSFR_evol))
    child_cond[childmatch] = 1L 
    print, parentchild[0:10]
    print, childmatch[0:10]
    notorphan = where(massSFR_evol.parent GE 0, n_notorphan)
    print, n_elements(childmatch), n_elements(parentmatch), n_notorphan, n_elements(massSFR_evol)
    massSFR_evol[childmatch].sfq = massSFR_data[parentmatch].sfq

; evolve the mass, SFR, and SSFR of the galaxies from z0 to zf
    parent_active_gal = where((strtrim(massSFR_data.sfq,2) EQ 'active') AND parent_cond, ngal_parent_active)
    parent_quiescent_gal = where((strtrim(massSFR_data.sfq,2) EQ 'quiescent') AND parent_cond, ngal_parent_quiescent)
    child_active_gal = where((strtrim(massSFR_evol.sfq,2) EQ 'active') AND child_cond, ngal_child_active)
    child_quiescent_gal = where((strtrim(massSFR_evol.sfq,2) EQ 'quiescent') AND child_cond, ngal_child_quiescent)
    
    print, ngal_parent_active, ngal_parent_quiescent
    print, ngal_child_active, ngal_child_quiescent
return 
end
