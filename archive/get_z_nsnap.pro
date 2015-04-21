function get_z_nsnap, nsnap, tcosmic=tcosmic
; returns redshift given snapshot number using Andrew's table 
    readcol, '/home/users/hahn/research/pro/tinker/central_quenching/snapshot_table.dat', $
        zi, aexp, redshift, t, t_wid
    tcosmic = t[closest(zi,nsnap)]
    return, redshift[closest(zi,nsnap)]
end
