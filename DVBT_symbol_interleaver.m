mmax8k = 8192;
nmax8k = 6048;
mmax2k = 2048;
nmax2k = 1512;
rip8k = zeros(mmax8k,13);
rip8k(3,1)=1;
rip8k(2,13) = 1;
rip2k = zeros(mmax2k,11);
rip2k(3,1)=1;
rip2k(2,11) = 1;
for i = 3:mmax8k-1
    rip8k(i+1,1:11) = rip8k(i,2:12);
    rip8k(i+1,12) = mod(rip8k(i,1)+rip8k(i,2)+rip8k(i,5)+rip8k(i,7),2);
    rip8k(i+1,13) = mod(i,2);
end
for i = 3:mmax2k-1
    rip2k(i+1,1:9) = rip2k(i,2:10);
    rip2k(i+1,10) = mod(rip2k(i,1)+rip2k(i,4),2);
    rip2k(i+1,11) = mod(i,2);
end
ri8k(:,[8,2,5,3,10,7,9,11,1,4,12,6,13]) = rip8k;
ri2k(:,[5,4,10,7,3,9,2,6,8,1,11]) = rip2k;
H8k = bi2de(ri8k,'right-msb');
H8k = H8k(H8k<nmax8k);
H2k = bi2de(ri2k,'right-msb');
H2k = H2k(H2k<nmax2k);

HB8k = reshape(reshape(mod((0:125)+[0 63 105 42 21 84].',126)*6+[1,4,2,5,3,6].',[],1)+(0:47)*756,[],1);