function get_nsnap_z, z, tcosmic=tcosmic
; returns n_snapshot given the redshift using Andrew's table 
    readcol, '/home/users/hahn/research/pro/tinker/wetzel_tree/snapshot_table.dat', $
        zi, aexp, redshift, t, t_wid
    tcosmic = t[closest(redshift,z)]
    return, long(zi[closest(redshift, z)])
end
