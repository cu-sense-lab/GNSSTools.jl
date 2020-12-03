"""
    acquire(data::GNSSSignal, T, type, prn; T_fine=10*T, sig_freq=missing,
            fd_center=0, fd_range=5000, fd_rate=0, threads=nthreads(), σω=10,
            iterations=1)
"""
function acquire(data::GNSSSignal, T, type, prn; T_fine=10*T, sig_freq=missing,
                 fd_center=0, fd_range=5000, fd_rate=0, threads=nthreads(),
                 σω=10, iterations=1)

end
