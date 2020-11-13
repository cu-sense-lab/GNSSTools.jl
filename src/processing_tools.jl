"""
    process(data, signal_type, prn)
"""
function process(signal, signal_type, prn; σω=1000., fd_center=0., fd_range=5000.,
                 RLM=10, replica_t_length=1e-3, cov_mult=1, q_a=1, q_mult=1,
                 dynamickf=true, dll_b=2, state_num=3, fd_rate=0.)
    f_s = signal.f_s;
    replica = definesignal(signal_type, f_s, replica_t_length);
    replicalong = definesignal(signal_type, f_s, RLM*replica_t_length);
    Δfd = 1/replica.t_length;  # Hz
    corr_result = gencorrresult(fd_range, Δfd, replica.sample_num);
    courseacquisition!(corr_result, signal, replica, prn;
                       fd_center=fd_center, fd_range=fd_range,
                       fd_rate=fd_rate, Δfd=Δfd);
    n0_est, fd_est, SNR_est = course_acq_est(corr_result, fd_center, fd_range,
                                             Δfd)
    results = fineacquisition(signal, replicalong, prn, fd_est,
                              n0_est, Val(:fft); σω=σω)
    trackresults = trackprn(signal, replica, prn, results.phi_init,
                        results.fd_est, results.n0_idx_course,
                        results.P, results.R; DLL_B=dll_b,
                        state_num=state_num, dynamickf=dynamickf,
                        cov_mult=cov_mult, qₐ=q_a, q_mult=q_mult);
    plotresults(trackresults)
end
