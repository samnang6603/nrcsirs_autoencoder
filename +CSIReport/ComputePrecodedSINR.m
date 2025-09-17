function sinr = ComputePrecodedSINR(H,W,nVar)

R = pagemtimes(H,W);
[~,S,V] = pagesvd(R,'econ','vector');

end