function x = color( n, alpha )

  hfa = zeros ( 2 * n, 1 );
  hfa(1) = 1.0; 
  for i = 2 : n
    hfa(i) = hfa(i-1) * ( 0.5 * alpha + ( i - 2 ) ) / ( i - 1 );
  end
  hfa(n+1:2*n) = 0.0; 

  wfa = rand ( n, 1 );

  wfa = 1 * ( 2.0 * wfa - 1.0 );

  z = zeros ( n, 1 );
  wfa = [ wfa; z ];

  fh = fft ( hfa );
  fw = fft ( wfa );

  fh = fh(1:n+1);
  fw = fw(1:n+1);

  fw = fh .* fw;    

  fw(1)   = fw(1)   / 2.0;
  fw(end) = fw(end) / 2.0;

  z = zeros ( n - 1, 1 );
  fw = [ fw; z ];

  x = ifft ( fw );

  x = 2.0 * real ( x(1:n) );

end
