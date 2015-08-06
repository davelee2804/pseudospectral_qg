#!/usr/bin/env octave

E = load( 'phi.50000.en' );
K = E(:,1);
X = E(:,2);
Y = E(:,3);
A = E(:,4);
B = E(:,5);

loglog( K, X, '-o', K, Y, '-+' );
print -dpng energySpectra_50000_k12_beta_3.png;
loglog( K, A, '-o', K, B, '-+', K(16:129), 10.0*K(16:129).^(-3.0), '-' );
print -dpng enstrophySpectra_50000_k12_beta_3.png;
