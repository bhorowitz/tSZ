USE linearpk
double precision :: linear_pk, k_ov_h
character(len=128) :: filename
integer :: n
external linear_pk
filename='wmap5baosn_max_likelihood_matterpower.dat'
n=896 ! # of lines in the file
CALL open_linearpk(filename,n)
k_ov_h=1d0
print*,'P(k) at k=',k_ov_h,' h Mpc^-1 is',linear_pk(k_ov_h),' h^3 Mpc^-3'
end
