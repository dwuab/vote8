function Pt=three_state_driver( p, eta, N, P0 )

filename = sprintf('data_trans_prob/N%d_p_%.3f_eta_%d_W12.dat',N,p,eta);
W12 = importdata(filename);
Pt=three_state_evolve( W12,p, eta, N, P0);
end