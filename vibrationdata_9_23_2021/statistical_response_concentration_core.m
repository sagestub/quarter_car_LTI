
% statistical_response_concentration_core.m  ver 1.0  by Tom Irvine

function[pure_tone,broadband]=statistical_response_concentration_core(fc,delta_f,nnet,mdens,D)


C=(pi*nnet*fc/2);

q=(C*mdens)*(2^D);
w=10*log10(q);

ratio= 1+(q-1)*(C/delta_f);

ba=10*log10(ratio);

pure_tone=w;
broadband=ba;
