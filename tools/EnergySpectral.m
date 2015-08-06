#!/usr/bin/env octave

E = load( 'phi.40000.en' );
K = E(:,1);
X = E(:,2);
Y = E(:,3);

loglog( K, X, '-o', K, Y, '-+', K, K.^(-3.0), '-' );
print -dpng energySpectra_40000.png;
